#' Habitat preferences of species and observers (Bayesian model)
#'
#' @description
#' Extracts from a model fitted by [fun_abund_hab_estim()] the habitat
#' selection parameters of the species (\eqn{S_{ih}}) or the habitat
#' preferences of the observers of the opportunistic dataset
#' (\eqn{q_{hk}}), as in Coron et al. (2018).
#'
#' These parameters are only identifiable relative to the first habitat
#' (\eqn{S_{i1} = 1}, \eqn{q_{1k} = 1}). Two scales are proposed:
#' * `scale = "ratio"` (default): ratio to the first (reference) habitat,
#'   \eqn{S_{ih} / S_{i1}}: a value of 2 means that, per unit area, the
#'   species is twice as dense in habitat \eqn{h} as in the reference habitat;
#' * `scale = "sum"`: rescaled to sum to one across habitats for each species,
#'   \eqn{S_{ih} / \sum_{h'} S_{ih'}}.
#'
#' @param fit A model fitted by [fun_abund_hab_estim()] with
#'   `model = "habitat"` or `"standardized_only"`.
#' @param who `"species"` (default) for the habitat selection of the
#'   species, or `"observers"` for the habitat preferences of the observers of
#'   the opportunistic dataset (`model = "habitat"` only).
#' @param scale `"ratio"` (default) or `"sum"`, see Description.
#' @param level Level of the credible intervals. Default 0.95.
#'
#' @return A data.frame with the columns `species` (only for
#'   `who = "species"`), `habitat`, `pref` (posterior mean), `sd`, `lower`
#'   and `upper` (equal-tailed credible interval).
#'
#' @references
#' Coron C., Calenge C., Giraud C. & Julliard R. (2018). Bayesian estimation
#' of species relative abundances and habitat preferences using opportunistic
#' data. *Environmental and Ecological Statistics*, 25(1), 71-93.
#' \doi{10.1007/s10651-018-0398-2}
#'
#' @seealso [fun_abund_hab_estim()], [fun_rel_abund()]
#' @export
#'
#' @examples
#' # see ?fun_abund_hab_estim
fun_hab_pref <- function(fit, who = c("species", "observers"),
                         scale = c("ratio", "sum"), level = 0.95) {
  who <- match.arg(who)
  scale <- match.arg(scale)
  if (!inherits(fit, "vgi4bio_bayes"))
    stop("'fit' must be the output of fun_abund_hab_estim().", call. = FALSE)
  if (fit$model == "no_habitat")
    stop("The 'no_habitat' model has no habitat parameter.", call. = FALSE)
  .check_level(level)
  H <- length(fit$habitats)
  M <- .draws(fit$samples)

  if (who == "observers") {
    if (fit$model != "habitat")
      stop("Observer preferences are only estimated by the 'habitat' model.", call. = FALSE)
    Q <- cbind(1, M[, .node_cols(M, "q2", H - 1), drop = FALSE])
    if (scale == "sum") Q <- Q / rowSums(Q)
    out <- data.frame(habitat = fit$habitats, .summ_draws(Q, level))
    names(out)[2] <- "pref"
    return(out)
  }

  I <- length(fit$species)
  out <- do.call(rbind, lapply(seq_len(I), function(i) {
    Si <- M[, paste0("S[", i, ",", seq_len(H), "]"), drop = FALSE]
    if (scale == "sum") Si <- Si / rowSums(Si)
    data.frame(species = fit$species[i], habitat = fit$habitats,
               .summ_draws(Si, level))
  }))
  names(out)[3] <- "pref"
  rownames(out) <- NULL
  out
}


## ---------------------------------------------------------------------------
## internal helpers
## ---------------------------------------------------------------------------

.check_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1)
    stop("'level' must be a number between 0 and 1.", call. = FALSE)
  invisible(TRUE)
}

## matrix of posterior draws (all chains) from a coda::mcmc.list, without
## relying on coda's as.matrix() method being registered (coda not loaded)
.draws <- function(samples) {
  do.call(rbind, lapply(samples, function(x) {
    x <- unclass(x)
    attr(x, "mcpar") <- NULL
    x
  }))
}

## column names of a vector node in a matrix of draws ("x" if length 1)
.node_cols <- function(M, node, n) {
  nm <- paste0(node, "[", seq_len(n), "]")
  if (n == 1 && !nm %in% colnames(M) && node %in% colnames(M)) nm <- node
  nm
}

## posterior summaries of the columns of a matrix of draws
.summ_draws <- function(D, level) {
  a <- (1 - level) / 2
  data.frame(mean = colMeans(D), sd = apply(D, 2, stats::sd),
             lower = apply(D, 2, stats::quantile, probs = a, names = FALSE),
             upper = apply(D, 2, stats::quantile, probs = 1 - a, names = FALSE),
             row.names = NULL)
}

## relative abundances from a Bayesian fit
.rel_abund_bayes <- function(fit, level, species, sites) {
  .check_level(level)
  sp_all <- fit$species
  st <- fit$sites
  I <- length(sp_all); J <- length(st)
  if (!is.null(species)) {
    unknown <- setdiff(species, sp_all)
    if (length(unknown))
      stop("Unknown species: ", paste(unknown, collapse = ", "), call. = FALSE)
  } else species <- sp_all

  ## site areas
  Vs <- NULL
  if (!is.null(sites)) {
    if (!is.data.frame(sites) || !"site" %in% names(sites))
      stop("'sites' must be a data.frame with a 'site' column.", call. = FALSE)
    m <- match(st, as.character(sites$site))
    if (anyNA(m))
      stop("Missing site(s) in 'sites': ", paste(st[is.na(m)], collapse = ", "),
           call. = FALSE)
    miss <- setdiff(fit$habitats, names(sites))
    if (fit$model != "no_habitat" && length(miss))
      stop("Missing habitat column(s) in 'sites': ", paste(miss, collapse = ", "),
           call. = FALSE)
    Vs <- as.matrix(sites[m, intersect(fit$habitats, names(sites)), drop = FALSE])
  } else if (fit$model != "no_habitat") {
    stop("'sites' (area of each habitat in each site) is required for the '",
         fit$model, "' model.", call. = FALSE)
  }

  M <- .draws(fit$samples)
  H <- length(fit$habitats)
  out <- do.call(rbind, lapply(match(species, sp_all), function(i) {
    Ni <- M[, paste0("N[", i, ",", seq_len(J), "]"), drop = FALSE]
    if (fit$model == "no_habitat") {
      w <- if (is.null(Vs)) rep(1, J) else rowSums(Vs)
      A <- sweep(Ni, 2, w, `*`)
    } else {
      Si <- M[, paste0("S[", i, ",", seq_len(H), "]"), drop = FALSE]
      A <- Ni * (Si %*% t(Vs))               # draws x sites
    }
    A <- A / rowSums(A)
    data.frame(species = sp_all[i], site = st, .summ_draws(A, level))
  }))
  names(out)[3] <- "rel_abund"
  out$in_standardized <- fit$monitored[as.character(out$species), fit$standardized]
  rownames(out) <- NULL
  out
}
