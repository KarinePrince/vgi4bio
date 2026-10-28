#' Simulate standardized and opportunistic counts with habitat selection
#'
#' @description
#' Simulates a data set following the model of Coron et al. (2018), with the
#' design of their simulation study (Section 3.1): `n_species` species,
#' `n_sites` sites, `n_hab` habitat types, a standardized dataset whose
#' cells each contain a single habitat and have a known effort of 1, and an
#' opportunistic dataset whose cells contain a mixture of habitats and have
#' an unknown effort. The count of species \eqn{i} in cell \eqn{c} (of site
#' \eqn{j}) for dataset \eqn{k} is drawn from
#' \deqn{\mathcal{P}\left(N_{ij} E_{ck} P_{ik} \sum_h
#'   \frac{q_{hk} V_{hc}}{\sum_{h'} q_{h'k} V_{h'c}}
#'   \frac{S_{ih}}{\sum_{h'} S_{ih'} V_{h'j}}\right).}
#'
#' The output has the format expected by [fun_abund_hab_estim()]: useful to
#' test the package, to illustrate it, or to check that a given design allows
#' the parameters to be recovered.
#'
#' @param n_species Number of species. Default 20.
#' @param n_sites Number of sites. Default 30.
#' @param n_hab Number of habitat types. Default 2.
#' @param n_std Number of standardized cells per site. Default 10.
#' @param n_opp Number of opportunistic cells per site. Default 30.
#' @param seed Optional seed for reproducibility.
#'
#' @return A list with:
#'   * `obs`: data.frame of counts (`species`, `cell`, `dataset`, `count`),
#'     datasets being `"STD"` (standardized) and `"OPP"` (opportunistic);
#'   * `cells`: data.frame describing the cells (`cell`, `dataset`, `site`,
#'     `effort`, and one column of area per habitat `hab1`, `hab2`, ...);
#'   * `sites`: data.frame with the area of each habitat in each site;
#'   * `truth`: list of the simulated parameters (`N`, `S`, `q`, `P`, `E`)
#'     and `rel_abund`, the true relative abundances rescaled to sum to one
#'     across sites for each species.
#'
#' @references
#' Coron C., Calenge C., Giraud C. & Julliard R. (2018). Bayesian estimation
#' of species relative abundances and habitat preferences using opportunistic
#' data. *Environmental and Ecological Statistics*, 25(1), 71-93.
#' \doi{10.1007/s10651-018-0398-2}
#'
#' @seealso [fun_abund_hab_estim()]
#' @importFrom stats rpois runif
#' @export
#'
#' @examples
#' sim <- fun_simul_hab(n_species = 5, n_sites = 8, seed = 1)
#' head(sim$obs)
#' head(sim$cells)
#' sim$sites
fun_simul_hab <- function(n_species = 20, n_sites = 30, n_hab = 2,
                          n_std = 10, n_opp = 30, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  I <- n_species; J <- n_sites; H <- n_hab
  if (H < 2) stop("'n_hab' must be at least 2.", call. = FALSE)
  hab <- paste0("hab", seq_len(H))
  sp <- sprintf("sp%02d", seq_len(I))
  st <- sprintf("site%02d", seq_len(J))

  ## cells: standardized cells contain one habitat, opportunistic cells a mixture
  n_c_std <- J * n_std
  n_c_opp <- J * n_opp
  U <- ceiling(stats::runif(n_c_std, 0, H))
  V <- rbind(t(vapply(U, function(u) as.numeric(seq_len(H) == u), numeric(H))),
             10 * matrix(stats::runif(n_c_opp * H), n_c_opp, H))
  site_c <- c(rep(seq_len(J), n_std), rep(seq_len(J), n_opp))
  ds_c <- rep(c(1L, 2L), c(n_c_std, n_c_opp))
  E <- c(rep(1, n_c_std), 100 * stats::runif(n_c_opp))
  Vsite <- 1000 * matrix(stats::runif(J * H), J, H)

  ## parameters
  S <- matrix(stats::runif(I * H), I, H)
  q <- cbind(1, stats::runif(H))
  N <- 10000 * matrix(stats::runif(I * J), I, J)
  P <- matrix(stats::runif(I * 2), I, 2)

  ## counts
  cells_id <- c(sprintf("std%04d", seq_len(n_c_std)), sprintf("opp%04d", seq_len(n_c_opp)))
  obs <- do.call(rbind, lapply(seq_along(cells_id), function(c) {
    k <- ds_c[c]; j <- site_c[c]
    w <- q[, k] * V[c, ] / sum(q[, k] * V[c, ])
    lam <- N[, j] * E[c] * P[, k] * drop((S %*% diag(w, H)) %*% rep(1, H)) /
      drop(S %*% Vsite[j, ])
    data.frame(species = sp, cell = cells_id[c],
               dataset = c("STD", "OPP")[k], count = stats::rpois(I, lam))
  }))

  cells <- data.frame(cell = cells_id, dataset = c("STD", "OPP")[ds_c],
                      site = st[site_c], effort = ifelse(ds_c == 1, 1, NA_real_))
  cells[hab] <- as.data.frame(V)
  sites <- data.frame(site = st)
  sites[hab] <- as.data.frame(Vsite)

  rel <- N / rowSums(N)
  truth_rel <- data.frame(species = rep(sp, each = J), site = rep(st, I),
                          rel_abund = as.vector(t(rel)))
  dimnames(N) <- list(sp, st); dimnames(S) <- list(sp, hab)
  dimnames(q) <- list(hab, c("STD", "OPP")); dimnames(P) <- list(sp, c("STD", "OPP"))

  list(obs = obs, cells = cells, sites = sites,
       truth = list(N = N, S = S, q = q, P = P, E = E, rel_abund = truth_rel))
}
