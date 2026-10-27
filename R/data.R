#' Bird counts in Aquitaine from standardized and opportunistic surveys
#'
#' Counts of 34 bird species on 63 sites of the Aquitaine region (France),
#' from the standardized ACT and STOC monitoring programmes and from
#' opportunistic records of the LPO. These are the data used in Giraud et al.
#' (2016) and Coron et al. (2018).
#'
#' @format A data.frame with 2142 rows (34 species x 63 sites) and 12 variables:
#' \describe{
#'   \item{Species}{factor, species code (first three letters of the genus
#'     followed by the first three letters of the species).}
#'   \item{Site}{factor, site code.}
#'   \item{Site_in_STOC}{integer (0/1), whether the site is also surveyed by
#'     the STOC programme.}
#'   \item{Species_in_ACT}{integer (0/1), whether the species is monitored
#'     by the ACT programme. For species with 0, `X_ACT` is a structural zero
#'     and must not be used as an observation.}
#'   \item{Effort_ACT}{numeric, number of point-count-years of the ACT
#'     programme on the site.}
#'   \item{Effort_STOC}{integer, sampling effort of the STOC programme on the
#'     site.}
#'   \item{Effort_LPO1, Effort_LPO2}{numeric, two alternative proxies of the
#'     opportunistic (LPO) sampling effort on the site (see Giraud et al.
#'     2016).}
#'   \item{X_ACT}{integer, count of the species on the site in the ACT
#'     programme.}
#'   \item{X_LPO}{integer, count of the species on the site in the
#'     opportunistic LPO data.}
#'   \item{X_ACTred}{integer, count of the species on the site in a reduced
#'     ACT dataset (subsample of `X_ACT`, see Giraud et al. 2016).}
#'   \item{X_STOC}{integer, count of the species on the site in the STOC
#'     programme (`NA` for sites not surveyed by STOC).}
#' }
#'
#' @source Giraud C., Calenge C., Coron C. & Julliard R. (2016). Capitalizing
#'   on opportunistic data for monitoring relative abundances of species.
#'   *Biometrics*, 72(2), 649-658. \doi{10.1111/biom.12431}
#'
#'   Coron C., Calenge C., Giraud C. & Julliard R. (2018). Bayesian estimation
#'   of species relative abundances and habitat preferences using
#'   opportunistic data. *Environmental and Ecological Statistics*, 25(1),
#'   71-93. \doi{10.1007/s10651-018-0398-2}
"Aquitaine"
