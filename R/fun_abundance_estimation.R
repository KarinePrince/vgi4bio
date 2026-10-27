#' Estimate relative species abundances by combining standardized and
#' opportunistic data
#'
#' @description
#' Fits the data-combination model of Giraud et al. (2016): a Poisson (or
#' quasi-Poisson) GLM in which the expected count of species \eqn{i} at site
#' \eqn{j} in dataset \eqn{k} is
#' \deqn{\log E(N_{ijk}) = \log(\mathrm{effort}_{jk}) + (species \times site)_{ij}
#'   + (species \times dataset)_{ik} + (site \times dataset)_{jk}.}
#' The `species:site` terms carry the relative abundances; the
#' `species:dataset` terms absorb the species detection probabilities in each
#' dataset, and the `site:dataset` terms the unknown sampling effort of the
#' opportunistic dataset(s).
#'
#' Identifiability requires the site effects of the standardized dataset to
#' be known (through `effort`): the `site:dataset` terms are therefore fixed
#' to zero for the dataset given in `standardized`, which is set as the
#' reference level of `dataset`. Several opportunistic datasets can be
#' combined; only one standardized dataset is currently supported.
#'
#' Use [fun_rel_abund()] on the result to extract the relative abundances of
#' each species on each site.
#'
#' @param df A data.frame of observed counts with (at least) the columns:
#'   * `species`: name or code of the species;
#'   * `site`: identifier of the observation site;
#'   * `dataset`: dataset (e.g. monitoring programme) the count comes from;
#'   * `effort`: sampling effort, strictly positive (e.g. number of
#'     point-count-years of the standardized programme at the site; set to 1
#'     for opportunistic data whose effort is unknown);
#'   * `N_obs`: observed count of `species` at `site` in `dataset`.
#'
#'   Only include, for each dataset, the species actually monitored by it:
#'   zeros for a species a programme does not record are structural zeros,
#'   not absences (a warning is issued if a species has only zeros in a
#'   dataset).
#' @param standardized Name of the standardized dataset (the value of
#'   `dataset` whose sampling effort is known), e.g. `"ACT"`. If `NULL`
#'   (default), the first level of `dataset` is used, with a warning: always
#'   set it explicitly.
#' @param family Error distribution of the GLM: `"quasipoisson"` (default) or
#'   `"poisson"`.
#' @param newdata Optional data.frame with columns `species`, `site`,
#'   `dataset` and `effort` for which predictions are computed. If `NULL`
#'   (default), no prediction is made.
#' @param type Type of prediction, passed to [stats::predict.glm()]: `"link"`
#'   (default, log scale), `"response"` (count scale) or `"terms"`.
#' @param verbose If `TRUE`, prints progress messages and the model summary.
#'   Default `FALSE`.
#'
#' @return A list with three elements:
#'   * `mod`: the fitted `glm` object;
#'   * `N_pred`: the predictions for `newdata`, or `NULL` if `newdata` is
#'     not provided;
#'   * `standardized`: the name of the standardized dataset.
#'
#' @references
#' Giraud C., Calenge C., Coron C. & Julliard R. (2016). Capitalizing on
#' opportunistic data for monitoring relative abundances of species.
#' *Biometrics*, 72(2), 649-658. \doi{10.1111/biom.12431}
#'
#' @author Karine Princé
#' @seealso [fun_rel_abund()], [fun_pearsons()], [fun_braycurt()], `vignette("tutorial", package = "vgi4bio")`
#' @importFrom stats glm predict relevel
#' @export
#'
#' @examples
#' data(Aquitaine)
#' # small subset to keep the example fast
#' # 4 species monitored by both programmes, 15 sites
#' sp <- c("ALAARV", "COLPAL", "TURMER", "STRDEC")
#' st <- levels(Aquitaine$Site)[1:15]
#' aq <- droplevels(Aquitaine[Aquitaine$Species %in% sp & Aquitaine$Site %in% st, ])
#' df <- rbind(
#'   data.frame(species = aq$Species, site = aq$Site, N_obs = aq$X_ACT,
#'              effort = aq$Effort_ACT, dataset = "ACT"),
#'   data.frame(species = aq$Species, site = aq$Site, N_obs = aq$X_LPO,
#'              effort = 1, dataset = "LPO"))
#' res <- fun_abund_estim(df, standardized = "ACT", family = "quasipoisson",
#'                        newdata = transform(df[df$dataset == "ACT", ], effort = 1),
#'                        type = "response")
#' head(res$N_pred)
#' head(fun_rel_abund(res))
fun_abund_estim <- function(df, standardized = NULL,
                            family = c("quasipoisson", "poisson"),
                            newdata = NULL,
                            type = c("link", "response", "terms"),
                            verbose = FALSE) {

  ## Check input data and parameters prior modelling
  family <- match.arg(family)
  type <- match.arg(type)

  if (!is.data.frame(df))
    stop("'df' must be a data.frame.", call. = FALSE)
  list_var <- c("species", "site", "N_obs", "dataset", "effort")
  missing_var <- setdiff(list_var, colnames(df))
  if (length(missing_var) > 0)
    stop("Missing column(s) in 'df': ", paste(missing_var, collapse = ", "),
         call. = FALSE)
  if (anyNA(df$N_obs))
    stop("'N_obs' contains NA values.", call. = FALSE)
  if (any(df$N_obs < 0))
    stop("'N_obs' contains negative values.", call. = FALSE)
  if (anyNA(df$effort) || any(df$effort <= 0))
    stop("'effort' must be strictly positive and without NA (it is used as log(effort) offset).",
         call. = FALSE)
  if (length(unique(df$dataset)) < 2)
    warning("'df' contains a single dataset: the model does not combine data sources.",
            call. = FALSE)

  df$species <- as.factor(df$species)
  df$site <- as.factor(df$site)
  df$dataset <- as.factor(df$dataset)

  ## Standardized dataset = reference level of 'dataset', so that its
  ## site:dataset terms are fixed to 0 (known effort, through the offset)
  if (is.null(standardized)) {
    standardized <- levels(df$dataset)[1]
    if (nlevels(df$dataset) > 1)
      warning("'standardized' not given: assuming '", standardized,
              "' is the standardized dataset (known effort). Set 'standardized' explicitly.",
              call. = FALSE)
  }
  if (length(standardized) != 1)
    stop("Only one standardized dataset is currently supported.", call. = FALSE)
  if (!standardized %in% levels(df$dataset))
    stop("'standardized' (", standardized, ") is not a value of 'dataset'.",
         call. = FALSE)
  df$dataset <- stats::relevel(df$dataset, ref = standardized)

  ## Species never observed in a dataset: probably not monitored by it
  ## (structural zeros), which makes the species:dataset estimates diverge.
  tot <- tapply(df$N_obs, list(df$species, df$dataset), sum)
  zero <- which(!is.na(tot) & tot == 0, arr.ind = TRUE)
  if (nrow(zero) > 0)
    warning("Some species have only zero counts in a dataset (",
            paste0(rownames(tot)[zero[, 1]], " in ", colnames(tot)[zero[, 2]],
                   collapse = ", "),
            "). If these species are not monitored by that dataset, remove the ",
            "corresponding rows from 'df' (structural zeros).", call. = FALSE)

  ## Fit model
  if (verbose) message("Fitting abundance model...")
  mod <- stats::glm(N_obs ~ species * dataset + species * site + site * dataset +
                      offset(log(effort)),
                    data = df, family = family)
  if (!mod$converged)
    warning("The GLM did not converge.", call. = FALSE)

  ## Predictions
  N_pred <- NULL
  if (!is.null(newdata)) {
    list_predvar <- c("species", "site", "dataset", "effort")
    missing_pred <- setdiff(list_predvar, colnames(newdata))
    if (length(missing_pred) > 0)
      stop("Missing column(s) in 'newdata': ", paste(missing_pred, collapse = ", "),
           call. = FALSE)
    if (verbose) message("Computing predictions...")
    newdata$dataset <- factor(newdata$dataset, levels = levels(df$dataset))
    N_pred <- stats::predict(mod, newdata = newdata, type = type)
  }

  if (verbose) print(summary(mod))
  list(mod = mod, N_pred = N_pred, standardized = standardized)
}


#' Pearson correlation between predicted and reference abundances
#'
#' For each species, computes the Pearson correlation coefficient across sites
#' between the abundances predicted by a given approach (e.g.
#' [fun_abund_estim()]) and "reference" abundances (e.g. from an independent
#' standardized survey), then returns the quartiles of these coefficients
#' across species. Species whose predictions sum to zero are ignored.
#'
#' @param N_pred Numeric vector of predicted abundances.
#' @param N_ref Numeric vector of reference abundances, in the same order as
#'   `N_pred`.
#' @param newdata The data.frame used to compute `N_pred` (must contain a
#'   `species` column and have one row per element of `N_pred`).
#'
#' @return A named numeric vector: first quartile, median and third quartile
#'   of the per-species Pearson correlation coefficients.
#' @seealso [fun_braycurt()]
#' @importFrom stats cor na.omit quantile
#' @export
#'
#' @examples
#' nd <- data.frame(species = rep(c("a", "b"), each = 4))
#' fun_pearsons(N_pred = c(1, 2, 3, 4, 2, 2, 5, 1),
#'              N_ref  = c(1, 3, 3, 5, 1, 2, 4, 1), newdata = nd)
fun_pearsons <- function(N_pred, N_ref, newdata) {
  .check_comparison_input(N_pred, N_ref, newdata)
  ns <- split(as.numeric(N_pred), newdata$species)
  n <- split(as.numeric(N_ref), newdata$species)
  r <- vapply(seq_along(ns), function(i) {
    if (sum(ns[[i]], na.rm = TRUE) > 0)
      suppressWarnings(stats::cor(ns[[i]], n[[i]], use = "complete.obs"))
    else NA_real_
  }, numeric(1))
  stats::quantile(stats::na.omit(r), c(0.25, 0.5, 0.75))
}


#' Bray-Curtis dissimilarity between predicted and reference abundances
#'
#' For each species, rescales predicted and reference abundances so that they
#' sum to one across sites, computes the Bray-Curtis dissimilarity between the
#' two distributions, then returns the quartiles of these dissimilarities
#' across species. Species whose predictions sum to zero are ignored.
#'
#' @inheritParams fun_pearsons
#'
#' @return A named numeric vector: first quartile, median and third quartile
#'   of the per-species Bray-Curtis dissimilarities.
#' @seealso [fun_pearsons()]
#' @importFrom stats na.omit quantile
#' @export
#'
#' @examples
#' nd <- data.frame(species = rep(c("a", "b"), each = 4))
#' fun_braycurt(N_pred = c(1, 2, 3, 4, 2, 2, 5, 1),
#'              N_ref  = c(1, 3, 3, 5, 1, 2, 4, 1), newdata = nd)
fun_braycurt <- function(N_pred, N_ref, newdata) {
  .check_comparison_input(N_pred, N_ref, newdata)
  ns <- lapply(split(as.numeric(N_pred), newdata$species), function(x) x / sum(x))
  n <- lapply(split(as.numeric(N_ref), newdata$species), function(x) x / sum(x))
  d <- vapply(seq_along(ns), function(i) {
    if (isTRUE(sum(ns[[i]]) > 0))
      sum(abs(ns[[i]] - n[[i]])) / sum(ns[[i]] + n[[i]])
    else NA_real_
  }, numeric(1))
  stats::quantile(stats::na.omit(d), c(0.25, 0.5, 0.75))
}


## internal: checks shared by fun_pearsons() and fun_braycurt()
.check_comparison_input <- function(N_pred, N_ref, newdata) {
  if (missing(newdata) || is.null(newdata))
    stop("'newdata' must be provided (the same data.frame as used for the predictions).",
         call. = FALSE)
  if (!"species" %in% colnames(newdata))
    stop("'newdata' must contain a 'species' column.", call. = FALSE)
  if (length(N_pred) != length(N_ref) || length(N_pred) != nrow(newdata))
    stop("'N_pred', 'N_ref' and 'newdata' must have the same length / number of rows.",
         call. = FALSE)
  invisible(TRUE)
}
