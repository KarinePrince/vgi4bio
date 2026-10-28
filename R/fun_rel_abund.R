#' Extract relative species abundances from a fitted combination model
#'
#' @description
#' Extracts from a model fitted by [fun_abund_estim()] the relative abundance
#' of each species on each site, i.e. the \eqn{\tilde N_{ij}} of Giraud et al.
#' (2016), with confidence intervals.
#'
#' The abundance of a species on a site is only identifiable up to a
#' species-specific constant. Two scales are therefore proposed:
#' * `scale = "sum"` (default): for each species, abundances are rescaled to
#'   sum to one across sites, \eqn{\tilde N_{ij} / \sum_{j'} \tilde N_{ij'}}
#'   (the standardisation used by Giraud et al. 2016 to compare methods).
#'   Available for all species, including those recorded only in the
#'   opportunistic dataset(s).
#' * `scale = "none"`: expected count per unit of effort of the standardized
#'   dataset (e.g. per point-count-year of ACT). Only defined for species
#'   recorded in the standardized dataset; `NA` for the others.
#'
#' Confidence intervals are computed on the log scale (delta method for
#' `scale = "sum"`) using the covariance matrix of the GLM, which includes
#' the overdispersion when `family = "quasipoisson"`, then back-transformed.
#'
#' For a Bayesian model fitted by [fun_abund_hab_estim()], the relative
#' abundance of species \eqn{i} on site \eqn{j} is computed for each
#' posterior sample as \eqn{N_{ij} \sum_h S_{ih} V_{hj}} (Coron et al. 2018),
#' where \eqn{V_{hj}} is the area of habitat \eqn{h} in site \eqn{j} given
#' in `sites` (\eqn{N_{ij} V_j} for the `"no_habitat"` model), then
#' rescaled to sum to one across sites. Only `scale = "sum"` is available.
#' The posterior mean, standard deviation and equal-tailed credible interval
#' are returned.
#'
#' @param fit The output of [fun_abund_estim()] or of [fun_abund_hab_estim()].
#' @param scale `"sum"` (default) or `"none"`, see Description.
#' @param level Confidence level of the intervals. Default 0.95.
#' @param species Optional character vector: species for which abundances are
#'   extracted. Default: all species of the model.
#' @param sites For a model fitted by [fun_abund_hab_estim()] only: a
#'   data.frame with a `site` column and one column per habitat giving the
#'   area of each habitat in each site (all the sites of the model). For the
#'   `"no_habitat"` model, if `NULL`, all sites are assumed to have the same
#'   area.
#'
#' @return A data.frame with one row per species x site and the columns:
#'   * `species`, `site`;
#'   * `rel_abund`: relative abundance on the chosen scale (posterior mean
#'     for a Bayesian model);
#'   * `se_log`: standard error of `log(rel_abund)` (GLM), or `sd`: posterior
#'     standard deviation of `rel_abund` (Bayesian model);
#'   * `lower`, `upper`: confidence (GLM) or credible (Bayesian) interval;
#'   * `in_standardized`: whether the species is recorded in the
#'     standardized dataset.
#'
#' @references
#' Giraud C., Calenge C., Coron C. & Julliard R. (2016). Capitalizing on
#' opportunistic data for monitoring relative abundances of species.
#' *Biometrics*, 72(2), 649-658. \doi{10.1111/biom.12431}
#'
#' Coron C., Calenge C., Giraud C. & Julliard R. (2018). Bayesian estimation
#' of species relative abundances and habitat preferences using opportunistic
#' data. *Environmental and Ecological Statistics*, 25(1), 71-93.
#' \doi{10.1007/s10651-018-0398-2}
#'
#' @seealso [fun_abund_estim()], [fun_abund_hab_estim()], [fun_hab_pref()],
#'   [fun_pearsons()], [fun_braycurt()]
#' @importFrom stats coef delete.response model.frame model.matrix qnorm terms vcov
#' @export
#'
#' @examples
#' data(Aquitaine)
#' sp <- c("ALAARV", "COLPAL", "TURMER", "STRDEC")
#' st <- levels(Aquitaine$Site)[1:15]
#' aq <- droplevels(Aquitaine[Aquitaine$Species %in% sp & Aquitaine$Site %in% st, ])
#' df <- rbind(
#'   data.frame(species = aq$Species, site = aq$Site, N_obs = aq$X_ACT,
#'              effort = aq$Effort_ACT, dataset = "ACT"),
#'   data.frame(species = aq$Species, site = aq$Site, N_obs = aq$X_LPO,
#'              effort = 1, dataset = "LPO"))
#' res <- fun_abund_estim(df, standardized = "ACT")
#' ra <- fun_rel_abund(res)
#' head(ra)
#' # relative abundances sum to one for each species
#' tapply(ra$rel_abund, ra$species, sum)
fun_rel_abund <- function(fit, scale = c("sum", "none"), level = 0.95,
                          species = NULL, sites = NULL) {
  scale <- match.arg(scale)
  if (inherits(fit, "vgi4bio_bayes")) {
    if (scale != "sum")
      stop("Only scale = \"sum\" is available for Bayesian models.", call. = FALSE)
    return(.rel_abund_bayes(fit, level = level, species = species, sites = sites))
  }
  if (!is.list(fit) || !inherits(fit$mod, "glm") || is.null(fit$standardized))
    stop("'fit' must be the output of fun_abund_estim() or fun_abund_hab_estim().",
         call. = FALSE)
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1)
    stop("'level' must be a number between 0 and 1.", call. = FALSE)

  mod <- fit$mod
  dat <- mod$data
  std <- fit$standardized
  sp <- levels(dat$species)
  if (!is.null(species)) {
    unknown <- setdiff(species, sp)
    if (length(unknown) > 0)
      stop("Unknown species: ", paste(unknown, collapse = ", "), call. = FALSE)
    sp <- species
  }
  st <- levels(dat$site)
  in_std <- sp %in% unique(as.character(dat$species[dat$dataset == std]))

  ## design matrix of the abundance part: standardized dataset, unit effort
  ## (its site:dataset terms are 0 by construction)
  nd <- expand.grid(site = st, species = sp, KEEP.OUT.ATTRS = FALSE)
  nd <- nd[, c("species", "site")]
  nd$dataset <- factor(std, levels = levels(dat$dataset))
  nd$effort <- 1
  tt <- stats::delete.response(stats::terms(mod))
  mf <- stats::model.frame(tt, nd, xlev = mod$xlevels)
  X <- stats::model.matrix(tt, mf, contrasts.arg = mod$contrasts)

  b <- stats::coef(mod)
  ok <- !is.na(b)                      # aliased coefficients are dropped
  X <- X[, names(b)[ok], drop = FALSE]
  V <- stats::vcov(mod)[ok, ok, drop = FALSE]
  eta <- drop(X %*% b[ok])

  rel <- se <- numeric(nrow(nd))
  for (s in seq_along(sp)) {
    i <- which(nd$species == sp[s])
    Xi <- X[i, , drop = FALSE]
    Sigma <- Xi %*% V %*% t(Xi)        # covariance of log-abundances
    if (scale == "sum") {
      e <- exp(eta[i] - max(eta[i]))
      p <- e / sum(e)
      G <- diag(length(i)) - matrix(p, length(i), length(i), byrow = TRUE)
      rel[i] <- p
      se[i] <- sqrt(pmax(rowSums((G %*% Sigma) * G), 0))
    } else {
      rel[i] <- if (in_std[s]) exp(eta[i]) else NA_real_
      se[i] <- if (in_std[s]) sqrt(pmax(diag(Sigma), 0)) else NA_real_
    }
  }
  if (scale == "none" && any(!in_std))
    message("Species not recorded in the standardized dataset get NA with scale = \"none\": ",
            paste(sp[!in_std], collapse = ", "))

  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(species = nd$species, site = nd$site,
                    rel_abund = rel, se_log = se,
                    lower = rel * exp(-z * se), upper = rel * exp(z * se),
                    in_standardized = rep(in_std, each = length(st)))
  rownames(out) <- NULL
  out
}
