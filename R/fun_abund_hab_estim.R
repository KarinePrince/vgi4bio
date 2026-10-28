#' Estimate relative abundances and habitat preferences (Bayesian model)
#'
#' @description
#' Fits with JAGS the model of Coron et al. (2018), which combines a
#' standardized dataset (known sampling effort) and an opportunistic dataset
#' (unknown sampling effort) and accounts for the habitat preferences of both
#' the species and the observers. The count of species \eqn{i} in cell
#' \eqn{c} of site \eqn{j} for dataset \eqn{k} follows
#' \deqn{X_{ick} \sim \mathcal{P}\left(N_{ij} E_{ck} P_{ik}
#'   \frac{\sum_h q_{hk} S_{ih} V_{hc}}{\sum_h q_{hk} V_{hc}}\right)}
#' where \eqn{V_{hc}} is the area of habitat \eqn{h} in cell \eqn{c},
#' \eqn{E_{ck}} the sampling effort (known for the standardized dataset),
#' \eqn{P_{ik}} the detection/reporting probability, \eqn{S_{ih}} the
#' selection of habitat \eqn{h} by species \eqn{i} and \eqn{q_{hk}} the
#' preference of the observers of dataset \eqn{k} for habitat \eqn{h}
#' (equation 1 of Coron et al. 2018).
#'
#' Three models are available:
#' * `"habitat"`: the model above (Coron et al. 2018);
#' * `"no_habitat"`: the same model without habitat preferences
#'   (\eqn{S_{ih} = q_{hk} = 1}), i.e. the model of Giraud et al. (2016)
#'   at the cell level;
#' * `"standardized_only"`: the habitat model fitted on the standardized
#'   data only.
#'
#' @details
#' **Unknown opportunistic efforts.** With `effort = "marginalized"`
#' (default), the efforts \eqn{E_{ck}} of the opportunistic cells are
#' integrated out analytically: conditionally on the total count
#' \eqn{T_c = \sum_i X_{ic}} of an opportunistic cell, the counts of the
#' species follow a multinomial distribution with probabilities
#' \eqn{\mu_{ic} / \sum_{i'} \mu_{i'c}}, where \eqn{\mu_{ic}} is the
#' intensity above without \eqn{E_{ck}}, and \eqn{T_c} only carries
#' information on \eqn{E_{ck}}. Under the authors' prior on \eqn{E_{ck}}
#' (Gamma(1e-11, 1e-11), i.e. a scale-invariant prior \eqn{\propto 1/E}),
#' the posterior distribution of all the other parameters is the same as
#' with the Poisson model, but the MCMC no longer has to explore one effort
#' per opportunistic cell, which strongly correlate with the abundances and
#' the observer preferences: convergence is much faster. Cells with a total
#' count of zero carry no information and are dropped.
#' With `effort = "estimated"`, the Poisson model is fitted with one effort
#' parameter per opportunistic cell, as in the authors' code (but with the
#' normalisation by \eqn{\sum_h q_{hk} V_{hc}} of equation 1, which mixes
#' slightly better than the authors' form \eqn{\sum_h q_{hk} S_{ih}
#' V_{hc} / V_c}; both forms give the same relative abundances and habitat
#' preferences). This option converges much more slowly and is mainly kept
#' for comparison.
#'
#' **MCMC engine.** With `engine = "jags"`, the model is fitted by JAGS, as
#' by the authors. With `engine = "nimble_hmc"`, it is fitted by the
#' No-U-Turn sampler (NUTS) of nimble/nimbleHMC (de Valpine et al. 2017): the
#' posterior distribution is the same (positive parameters are sampled on the
#' log scale with the same priors), but all parameters are updated jointly
#' using the gradient of the log-posterior. The abundances and the habitat
#' selection parameters compensate each other (\eqn{N_{ij} S_{ih}}, and
#' globally the observer preferences \eqn{q_{hk}} against all the
#' \eqn{S_{ih}}), which makes the univariate updates of JAGS mix very
#' slowly: on simulated data (10 species, 15 sites), JAGS needed about 80000
#' iterations per chain to converge, NUTS a few hundred, for a wall-clock
#' time about 4 times shorter. `"nimble_hmc"` is therefore recommended,
#' with much smaller `n.burnin` and `n.iter` (e.g. 500 and 1000, `thin = 1`).
#' It requires a C++ compiler; compilation takes one to a few minutes.
#' Whatever the engine, always check the convergence (`rhat` < 1.1).
#'
#' **Computing time.** On the data of Coron et al. (2018) (34 species, 65
#' sites, about 1500 opportunistic cells), one JAGS iteration takes about 4
#' seconds and one NUTS iteration (after warm-up) at least 20 seconds: a
#' converged fit takes many hours and should be run on a dedicated machine.
#'
#' **Identifiability constraints** (as in the code of Coron et al.):
#' \eqn{P_{i,opp} = 1} for all species, \eqn{P_{i,std} = 1} for the first
#' species monitored by the standardized dataset, \eqn{S_{i1} = 1} (first
#' habitat as reference), \eqn{q_{h,std} = 1} for all habitats (no observer
#' preference in the standardized protocol) and \eqn{q_{1,opp} = 1}.
#'
#' **Priors** are those of the authors' JAGS code (`ModelWithHabitat.txt`):
#' Gamma(1e-11, 1e-11) on \eqn{N} and on the opportunistic efforts,
#' Gamma(0.01, 0.01) on \eqn{P}, Gamma(1e-4, 1e-4) on \eqn{S} and \eqn{q}.
#'
#' @param obs data.frame of counts with columns `species`, `cell`, `dataset`
#'   and `count`. A species is considered as monitored by a dataset if it has
#'   at least one row (possibly with a zero count) for this dataset; for these
#'   species, missing cells of the dataset are filled with zero counts. Every
#'   species must be monitored by the opportunistic dataset.
#' @param cells data.frame describing the cells, with columns `cell`,
#'   `dataset`, `site`, `effort` (sampling effort, required and strictly
#'   positive for the standardized dataset, ignored otherwise) and one column
#'   per habitat giving the area of the habitat in the cell.
#' @param standardized Name of the standardized dataset.
#' @param habitats Names of the habitat columns of `cells`. The first habitat
#'   is the reference for the habitat selection parameters.
#' @param model `"habitat"` (default), `"no_habitat"` or
#'   `"standardized_only"`, see Description.
#' @param effort `"marginalized"` (default) or `"estimated"`: how the
#'   unknown efforts of the opportunistic cells are handled, see Details.
#' @param engine `"jags"` (default) or `"nimble_hmc"`: MCMC software used,
#'   see Details. `"nimble_hmc"` requires `effort = "marginalized"`.
#' @param n.chains Number of MCMC chains. Default 3 (the authors used 1; at
#'   least 2 are needed for the Gelman-Rubin convergence diagnostic).
#' @param n.adapt Number of adaptation iterations (JAGS only; NUTS adapts
#'   during the burn-in, which is its warm-up). Default 1000.
#' @param n.burnin Number of burn-in iterations, discarded. Default 1000.
#' @param n.iter Number of iterations kept (before thinning). Default 10000.
#' @param thin Thinning interval. Default 10.
#' @param seed Optional seed, for reproducible results.
#' @param quiet If `TRUE` (default), JAGS progress is not displayed.
#'
#' @return An object of class `"vgi4bio_bayes"`: a list with the posterior
#'   samples (`samples`, a `coda::mcmc.list`), the model name, the
#'   standardized and opportunistic dataset names, the levels of species,
#'   sites and habitats, the monitoring table, the convergence diagnostics
#'   (`rhat`, `n_eff`) and the data passed to JAGS. Use [fun_rel_abund()]
#'   and [fun_hab_pref()] to extract the estimates.
#'
#' @references
#' Coron C., Calenge C., Giraud C. & Julliard R. (2018). Bayesian estimation
#' of species relative abundances and habitat preferences using opportunistic
#' data. *Environmental and Ecological Statistics*, 25(1), 71-93.
#' \doi{10.1007/s10651-018-0398-2}
#'
#' Giraud C., Calenge C., Coron C. & Julliard R. (2016). Capitalizing on
#' opportunistic data for monitoring relative abundances of species.
#' *Biometrics*, 72(2), 649-658. \doi{10.1111/biom.12431}
#'
#' de Valpine P., Turek D., Paciorek C.J., Anderson-Bergman C., Temple Lang D.
#' & Bodik R. (2017). Programming with models: writing statistical
#' algorithms for general model structures with NIMBLE. *Journal of
#' Computational and Graphical Statistics*, 26(2), 403-413.
#' \doi{10.1080/10618600.2016.1172487}
#'
#' @seealso [fun_simul_hab()], [fun_rel_abund()], [fun_hab_pref()]
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("rjags", quietly = TRUE)) {
#'   sim <- fun_simul_hab(n_species = 4, n_sites = 6, n_std = 5, n_opp = 8, seed = 1)
#'   fit <- fun_abund_hab_estim(sim$obs, sim$cells, standardized = "STD",
#'                              habitats = c("hab1", "hab2"),
#'                              n.adapt = 200, n.burnin = 200, n.iter = 1000,
#'                              thin = 2, seed = 1)
#'   fit
#'   head(fun_rel_abund(fit, sites = sim$sites))
#'   fun_hab_pref(fit)
#' }
#' }
fun_abund_hab_estim <- function(obs, cells, standardized, habitats,
                                model = c("habitat", "no_habitat", "standardized_only"),
                                effort = c("marginalized", "estimated"),
                                engine = c("jags", "nimble_hmc"),
                                n.chains = 3, n.adapt = 1000, n.burnin = 1000,
                                n.iter = 10000, thin = 10, seed = NULL,
                                quiet = TRUE) {
  model <- match.arg(model)
  effort <- match.arg(effort)
  engine <- match.arg(engine)
  if (engine == "jags" && !requireNamespace("rjags", quietly = TRUE))
    stop("Package 'rjags' (and the JAGS software, https://mcmc-jags.sourceforge.io) ",
         "is required for engine = \"jags\".", call. = FALSE)
  if (engine != "jags" && effort == "estimated")
    stop("engine = \"nimble_hmc\" is only available with effort = \"marginalized\".",
         call. = FALSE)

  dat <- .prep_hab_data(obs, cells, standardized, habitats, model, effort)
  jd <- dat$jags_data
  inits <- .jags_inits(jd, n.chains, seed)

  if (engine == "jags") {
    code <- .jags_code(model, effort)
    monitor <- c("N", "P1", "S", "q2", "E2")
    monitor <- monitor[vapply(monitor, function(m) .jags_has(m, jd, model, effort),
                              logical(1))]
    progress <- if (quiet) "none" else "text"
    jm <- rjags::jags.model(textConnection(code), data = jd, inits = inits,
                            n.chains = n.chains, n.adapt = n.adapt, quiet = quiet)
    if (n.burnin > 0) stats::update(jm, n.iter = n.burnin, progress.bar = progress)
    samples <- rjags::coda.samples(jm, variable.names = monitor, n.iter = n.iter,
                                   thin = thin, progress.bar = progress)
  } else {
    samples <- .fit_nimble(jd, model, n.chains, n.burnin, n.iter, thin, seed,
                           quiet, inits)
  }

  diag <- .mcmc_diag(samples)
  if (is.finite(diag$max_rhat) && diag$max_rhat > 1.1)
    warning("Some parameters have not converged (max Rhat = ",
            round(diag$max_rhat, 2), "): increase 'n.burnin' and 'n.iter'.",
            call. = FALSE)

  structure(list(samples = samples, model = model, effort = effort, engine = engine,
                 standardized = standardized, opportunistic = dat$opportunistic,
                 species = dat$species, sites = dat$sites, habitats = habitats,
                 monitored = dat$monitored, rhat = diag$rhat, n_eff = diag$n_eff,
                 jags_data = jd, call = match.call()),
            class = "vgi4bio_bayes")
}


#' @export
print.vgi4bio_bayes <- function(x, ...) {
  cat("Bayesian model of relative abundances (vgi4bio)\n")
  cat("  model:        ", x$model,
      if (x$model != "standardized_only") paste0("(opportunistic effort ", x$effort, ")"),
      "\n")
  cat("  datasets:     ", x$standardized, "(standardized)",
      if (!is.null(x$opportunistic)) paste0("+ ", x$opportunistic, " (opportunistic)"), "\n")
  cat("  species:      ", length(x$species), "\n")
  cat("  sites:        ", length(x$sites), "\n")
  if (x$model != "no_habitat")
    cat("  habitats:     ", paste(x$habitats, collapse = ", "), "\n")
  n_ch <- coda::nchain(x$samples)
  cat("  engine:       ", if (is.null(x$engine)) "jags" else x$engine, "\n")
  cat("  MCMC:         ", n_ch, "chain(s) x", coda::niter(x$samples), "samples\n")
  if (n_ch > 1)
    cat("  max Rhat:     ", round(max(x$rhat, na.rm = TRUE), 3),
        "(should be < 1.1)\n")
  cat("  min eff. size:", round(min(x$n_eff, na.rm = TRUE)), "\n")
  invisible(x)
}


## ---------------------------------------------------------------------------
## internal helpers
## ---------------------------------------------------------------------------

.prep_hab_data <- function(obs, cells, standardized, habitats, model, effort) {
  if (!is.data.frame(obs) || !is.data.frame(cells))
    stop("'obs' and 'cells' must be data.frames.", call. = FALSE)
  miss <- setdiff(c("species", "cell", "dataset", "count"), names(obs))
  if (length(miss)) stop("Missing column(s) in 'obs': ", paste(miss, collapse = ", "), call. = FALSE)
  hab_cols <- if (model == "no_habitat") character(0) else habitats
  miss <- setdiff(c("cell", "dataset", "site", "effort", hab_cols), names(cells))
  if (length(miss)) stop("Missing column(s) in 'cells': ", paste(miss, collapse = ", "), call. = FALSE)
  if (model != "no_habitat" && length(habitats) < 2)
    stop("At least two habitats are needed.", call. = FALSE)
  if (anyNA(obs$count) || any(obs$count < 0))
    stop("'count' must be non-negative, without NA.", call. = FALSE)

  obs$species <- as.character(obs$species)
  obs$cell <- as.character(obs$cell)
  obs$dataset <- as.character(obs$dataset)
  cells$cell <- as.character(cells$cell)
  cells$dataset <- as.character(cells$dataset)
  cells$site <- as.character(cells$site)

  ds <- unique(obs$dataset)
  if (length(standardized) != 1 || !standardized %in% ds)
    stop("'standardized' must be one value of obs$dataset.", call. = FALSE)
  opp <- setdiff(ds, standardized)
  if (model == "standardized_only") {
    obs <- obs[obs$dataset == standardized, ]
    opp <- character(0)
  } else if (length(opp) != 1) {
    stop("The model currently handles exactly one standardized and one ",
         "opportunistic dataset.", call. = FALSE)
  }

  key_c <- paste(cells$cell, cells$dataset, sep = "\r")
  if (anyDuplicated(key_c))
    stop("Duplicated cell x dataset combinations in 'cells'.", call. = FALSE)
  key_o <- paste(obs$cell, obs$dataset, sep = "\r")
  if (!all(key_o %in% key_c))
    stop("Some cells of 'obs' are not described in 'cells'.", call. = FALSE)
  if (anyDuplicated(paste(obs$species, key_o)))
    stop("Duplicated species x cell x dataset combinations in 'obs'.", call. = FALSE)

  cells <- cells[cells$dataset %in% c(standardized, opp), ]
  if (model != "no_habitat") {
    Vc <- as.matrix(cells[, habitats])
    if (anyNA(Vc) || any(Vc < 0)) stop("Habitat areas must be non-negative, without NA.", call. = FALSE)
    if (any(rowSums(Vc) <= 0)) stop("Some cells have a null total habitat area.", call. = FALSE)
  }
  is_std <- cells$dataset == standardized
  if (anyNA(cells$effort[is_std]) || any(cells$effort[is_std] <= 0))
    stop("'effort' must be strictly positive for the cells of the standardized dataset.",
         call. = FALSE)

  species <- sort(unique(obs$species))
  sites <- sort(unique(cells$site))
  dsets <- c(standardized, opp)
  monitored <- matrix(vapply(dsets, function(d) species %in% obs$species[obs$dataset == d],
                             logical(length(species))),
                      ncol = length(dsets), dimnames = list(species, dsets))
  if (length(opp) && !all(monitored[, opp]))
    stop("All species must be monitored by the opportunistic dataset (",
         paste(species[!monitored[, opp]], collapse = ", "), " are not).", call. = FALSE)
  if (!any(monitored[, standardized]))
    stop("No species is monitored by the standardized dataset.", call. = FALSE)

  ## complete count table: monitored species x cells of the dataset, zeros filled
  full <- do.call(rbind, lapply(dsets, function(d) {
    expand.grid(species = species[monitored[, d]],
                cell = cells$cell[cells$dataset == d], dataset = d,
                stringsAsFactors = FALSE)
  }))
  full$count <- 0
  m <- match(paste(full$species, full$cell, full$dataset, sep = "\r"),
             paste(obs$species, obs$cell, obs$dataset, sep = "\r"))
  full$count[!is.na(m)] <- obs$count[m[!is.na(m)]]
  tot_sp <- tapply(full$count, list(full$species, full$dataset), sum)
  if (any(tot_sp == 0, na.rm = TRUE))
    warning("Some species have only zero counts in a dataset: check that they are ",
            "monitored by it (structural zeros).", call. = FALSE)

  I <- length(species)
  prop <- if (model != "no_habitat") {
    Vc <- as.matrix(cells[, habitats]); Vc / rowSums(Vc)
  } else NULL

  ## standardized cells: Poisson, long format
  cs <- cells[is_std, ]
  fs <- full[full$dataset == standardized, ]
  ci <- match(fs$cell, cs$cell)
  jd <- list(I = I, J = length(sites), Ls = nrow(fs), ys = fs$count,
             sps = match(fs$species, species), sts = match(cs$site, sites)[ci],
             cs = ci, Erel = cs$effort / cs$effort[1])
  if (!is.null(prop)) {
    jd$H <- length(habitats)
    jd$props <- prop[is_std, , drop = FALSE]
  }

  if (length(opp)) {
    ## detection P[i, std]: free for monitored species except the reference one
    mon_std <- which(monitored[, standardized])
    free <- mon_std[-1]
    pidx <- rep(length(free) + 1L, I)
    pidx[free] <- seq_along(free)
    jd$nP <- length(free)
    jd$pidx <- pidx

    co <- cells[!is_std, ]
    fo <- full[full$dataset == opp, ]
    Y <- matrix(0, nrow(co), I)
    Y[cbind(match(fo$cell, co$cell), match(fo$species, species))] <- fo$count
    propo <- if (!is.null(prop)) prop[!is_std, , drop = FALSE] else NULL
    sto <- match(co$site, sites)
    if (effort == "marginalized") {
      keep <- rowSums(Y) > 0               # empty cells carry no information
      jd$Co <- sum(keep)
      jd$Yo <- Y[keep, , drop = FALSE]
      jd$To <- rowSums(Y)[keep]
      jd$sto <- sto[keep]
      if (!is.null(propo)) jd$propo <- propo[keep, , drop = FALSE]
    } else {
      idx <- which(Y >= 0, arr.ind = TRUE)  # all cells x species
      jd$Co <- nrow(co)
      jd$Lo <- nrow(idx)
      jd$yo <- Y[idx]
      jd$spo <- idx[, 2]
      jd$co <- idx[, 1]
      jd$sto <- sto
      if (!is.null(propo)) jd$propo <- propo
    }
  }
  list(jags_data = jd, species = species, sites = sites, monitored = monitored,
       opportunistic = if (length(opp)) opp else NULL)
}


.jags_code <- function(model, effort) {
  hab <- model != "no_habitat"
  opp <- model != "standardized_only"
  code <- "model {
  for (i in 1:I) { for (j in 1:J) { N[i, j] ~ dgamma(1.0E-11, 1.0E-11) } }"
  if (hab) code <- paste0(code, "
  for (i in 1:I) {
    S[i, 1] <- 1
    for (h in 2:H) { S[i, h] ~ dgamma(1.0E-4, 1.0E-4) }
  }")
  if (opp) code <- paste0(code, "
  # detection in the standardized dataset (1 in the opportunistic one)
  P1[nP + 1] <- 1
  for (m in 1:nP) { P1[m] ~ dgamma(0.01, 0.01) }
  for (i in 1:I) { P[i] <- P1[pidx[i]] }")
  if (opp && hab) code <- paste0(code, "
  # observer preferences in the opportunistic dataset (1 in the standardized one)
  qo[1] <- 1
  for (h in 2:H) { q2[h - 1] ~ dgamma(1.0E-4, 1.0E-4)
                   qo[h] <- q2[h - 1] }")

  ## standardized data: Poisson with known effort (q = 1)
  pstd <- if (opp) " * P[sps[l]]" else ""
  code <- paste0(code, "
  for (l in 1:Ls) {
    ys[l] ~ dpois(N[sps[l], sts[l]] * Erel[cs[l]]", pstd,
    if (hab) " *\n                  inprod(S[sps[l], 1:H], props[cs[l], 1:H])", ")
  }")

  if (opp && effort == "marginalized") {
    ## the unknown effort of each opportunistic cell is integrated out: the
    ## counts of a cell are multinomial given their total (see ?fun_abund_hab_estim)
    mu <- if (hab) "N[i, sto[c]] * inprod(S[i, 1:H], qo[1:H] * propo[c, 1:H])"
          else "N[i, sto[c]]"
    code <- paste0(code, "
  for (c in 1:Co) {
    for (i in 1:I) { mu[c, i] <- ", mu, " }
    Yo[c, 1:I] ~ dmulti(mu[c, 1:I] / sum(mu[c, 1:I]), To[c])
  }")
  } else if (opp) {
    ## one effort per opportunistic cell; intensity normalised by the observer
    ## preference of the cell (eq. 1 of Coron et al. 2018), which mixes better
    ## than the authors' form sum_h q S V / V
    lam <- if (hab) "N[spo[l], sto[co[l]]] * E2[co[l]] *
                  inprod(S[spo[l], 1:H], qo[1:H] * propo[co[l], 1:H]) /
                  inprod(qo[1:H], propo[co[l], 1:H])"
           else "N[spo[l], sto[co[l]]] * E2[co[l]]"
    code <- paste0(code, "
  for (c in 1:Co) { E2[c] ~ dgamma(1.0E-11, 1.0E-11) }
  for (l in 1:Lo) {
    yo[l] ~ dpois(", lam, ")
  }")
  }
  paste0(code, "\n}\n")
}


.jags_has <- function(node, jd, model, effort) {
  switch(node,
         N = TRUE,
         S = model != "no_habitat",
         P1 = model != "standardized_only" && isTRUE(jd$nP > 0),
         q2 = model == "habitat",
         E2 = model != "standardized_only" && effort == "estimated")
}


.jags_inits <- function(jd, n.chains, seed) {
  ## crude starting values for N: mean standardized count per species x site
  tot <- tapply(jd$ys / jd$Erel[jd$cs],
                list(factor(jd$sps, 1:jd$I), factor(jd$sts, 1:jd$J)), mean)
  if (!is.null(jd$Yo)) {
    ## species seen only in the opportunistic data: use their relative frequency
    fo <- rowsum(jd$Yo, factor(jd$sto, 1:jd$J))
    fo <- t(fo / pmax(rowSums(fo), 1)) * max(tot, 1, na.rm = TRUE)
    tot[is.na(tot)] <- fo[is.na(tot)]
  }
  tot[is.na(tot)] <- 0
  if (!is.null(seed)) set.seed(seed)
  lapply(seq_len(n.chains), function(ch) {
    ini <- list(N = (tot + 0.1) * exp(stats::rnorm(length(tot), 0, 0.1)),
                .RNG.name = "base::Mersenne-Twister",
                .RNG.seed = if (is.null(seed)) sample.int(1e6, 1) else seed + ch)
    if (!is.null(jd$H)) ini$S <- cbind(NA, matrix(1, jd$I, jd$H - 1))
    if (!is.null(jd$nP) && jd$nP > 0) ini$P1 <- c(rep(1, jd$nP), NA)
    if (!is.null(jd$propo)) ini$q2 <- rep(1, jd$H - 1)
    if (!is.null(jd$Lo)) ini$E2 <- rep(1, jd$Co)
    ini
  })
}


.mcmc_diag <- function(samples) {
  ## constant nodes (fixed by the identifiability constraints) are ignored
  M <- .draws(samples)
  keep <- apply(M, 2, function(x) stats::var(x) > 0)
  samples <- samples[, keep, drop = FALSE]
  n_eff <- coda::effectiveSize(samples)
  rhat <- rep(NA_real_, length(n_eff))
  if (coda::nchain(samples) > 1) {
    rhat <- tryCatch(coda::gelman.diag(samples, autoburnin = FALSE,
                                       multivariate = FALSE)$psrf[, 1],
                     error = function(e) rep(NA_real_, length(n_eff)))
  }
  names(rhat) <- names(n_eff)
  max_rhat <- suppressWarnings(max(rhat[is.finite(rhat)]))
  list(rhat = rhat, n_eff = n_eff, max_rhat = max_rhat)
}
