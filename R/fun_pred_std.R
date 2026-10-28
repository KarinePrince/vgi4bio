#' Predict counts of an independent standardized survey (Bayesian model)
#'
#' @description
#' From a model fitted by [fun_abund_hab_estim()], predicts for each species
#' and each site the count per unit of effort expected for an independent
#' standardized survey (e.g. the STOC programme in Coron et al. 2018), given
#' the location, effort and habitat composition of its sampling cells
#' (Appendix A.3 of Coron et al. 2018):
#' \deqn{\hat X_{ij} = N_{ij} \frac{\sum_{c \in j} e_c \sum_h S_{ih}
#'   V_{hc} / V_c}{\sum_{c \in j} e_c}}
#' where \eqn{e_c} is the effort of the cell \eqn{c}. For the
#' `"no_habitat"` model, \eqn{\hat X_{ij} = N_{ij}}.
#'
#' The detection probability of the species in the new survey is unknown:
#' predictions are only comparable between sites for a given species (e.g.
#' through per-species correlations with [fun_pearsons()]).
#'
#' @param fit A model fitted by [fun_abund_hab_estim()].
#' @param cells data.frame describing the cells of the new survey, with
#'   columns `cell`, `site`, `effort` (strictly positive) and, except for the
#'   `"no_habitat"` model, one column per habitat of the model giving the
#'   area of the habitat in the cell. Cells of sites absent from the model
#'   are ignored (with a message).
#' @param level Level of the credible intervals. Default 0.95.
#'
#' @return A data.frame with the columns `species`, `site`, `pred`
#'   (posterior mean of the expected count per unit of effort, up to a
#'   species-specific constant), `sd`, `lower`, `upper` and
#'   `in_standardized`.
#'
#' @references
#' Coron C., Calenge C., Giraud C. & Julliard R. (2018). Bayesian estimation
#' of species relative abundances and habitat preferences using opportunistic
#' data. *Environmental and Ecological Statistics*, 25(1), 71-93.
#' \doi{10.1007/s10651-018-0398-2}
#'
#' @seealso [fun_abund_hab_estim()], [fun_pearsons()]
#' @export
#'
#' @examples
#' # see vignette("habitat", package = "vgi4bio")
fun_pred_std <- function(fit, cells, level = 0.95) {
  if (!inherits(fit, "vgi4bio_bayes"))
    stop("'fit' must be the output of fun_abund_hab_estim().", call. = FALSE)
  .check_level(level)
  hab <- if (fit$model == "no_habitat") character(0) else fit$habitats
  miss <- setdiff(c("cell", "site", "effort", hab), names(cells))
  if (length(miss))
    stop("Missing column(s) in 'cells': ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyNA(cells$effort) || any(cells$effort <= 0))
    stop("'effort' must be strictly positive.", call. = FALSE)
  cells$site <- as.character(cells$site)
  out_sites <- setdiff(unique(cells$site), fit$sites)
  if (length(out_sites)) {
    message("Cells of sites absent from the model are ignored: ",
            paste(out_sites, collapse = ", "))
    cells <- cells[cells$site %in% fit$sites, ]
  }
  if (!nrow(cells)) stop("No cell in the sites of the model.", call. = FALSE)
  if (length(hab)) {
    V <- as.matrix(cells[, hab])
    if (anyNA(V) || any(V < 0) || any(rowSums(V) <= 0))
      stop("Habitat areas must be non-negative, without NA, with a positive total.",
           call. = FALSE)
    prop <- V / rowSums(V)
  }

  M <- .draws(fit$samples)
  sts <- sort(unique(cells$site))
  jj <- match(sts, fit$sites)
  H <- length(fit$habitats)
  out <- do.call(rbind, lapply(seq_along(fit$species), function(i) {
    Ni <- M[, paste0("N[", i, ",", jj, "]"), drop = FALSE]
    if (length(hab)) {
      Si <- M[, paste0("S[", i, ",", seq_len(H), "]"), drop = FALSE]
      w <- vapply(sts, function(st) {
        k <- cells$site == st
        colSums(cells$effort[k] * prop[k, , drop = FALSE]) / sum(cells$effort[k])
      }, numeric(H))                                 # H x sites
      A <- Ni * (Si %*% w)
    } else {
      A <- Ni
    }
    data.frame(species = fit$species[i], site = sts, .summ_draws(A, level))
  }))
  names(out)[3] <- "pred"
  out$in_standardized <- fit$monitored[as.character(out$species), fit$standardized]
  rownames(out) <- NULL
  out
}
