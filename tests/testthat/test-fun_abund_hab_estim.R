## Bayesian model of Coron et al. (2018): small simulated data sets, short
## MCMC runs (structure of the outputs only, not convergence).

small_sim <- function() fun_simul_hab(n_species = 4, n_sites = 6, n_hab = 2,
                                      n_std = 5, n_opp = 8, seed = 1)
hab <- c("hab1", "hab2")

test_that("fun_simul_hab returns data in the expected format", {
  sim <- small_sim()
  expect_named(sim, c("obs", "cells", "sites", "truth"))
  expect_named(sim$obs, c("species", "cell", "dataset", "count"))
  expect_true(all(c("cell", "dataset", "site", "effort", hab) %in% names(sim$cells)))
  expect_equal(nrow(sim$obs), 4 * 6 * (5 + 8))
  expect_equal(nrow(sim$sites), 6)
  expect_equal(as.numeric(tapply(sim$truth$rel_abund$rel_abund,
                                 sim$truth$rel_abund$species, sum)), rep(1, 4))
  expect_error(fun_simul_hab(n_hab = 1), "n_hab")
})

test_that("input data are checked", {
  sim <- small_sim()
  prep <- function(...) vgi4bio:::.prep_hab_data(..., model = "habitat",
                                                 effort = "marginalized")
  expect_error(prep(sim$obs[, -4], sim$cells, "STD", hab), "count")
  expect_error(prep(sim$obs, sim$cells[, -5], "STD", hab), "hab1")
  expect_error(prep(sim$obs, sim$cells, "XXX", hab), "standardized")
  c2 <- sim$cells; c2$effort[1] <- NA
  expect_error(prep(sim$obs, c2, "STD", hab), "strictly positive")
  o2 <- sim$obs; o2$count[1] <- -1
  expect_error(prep(o2, sim$cells, "STD", hab), "non-negative")
  o3 <- sim$obs[!(sim$obs$species == "sp01" & sim$obs$dataset == "OPP"), ]
  expect_error(prep(o3, sim$cells, "STD", hab), "opportunistic")
  o4 <- rbind(sim$obs, sim$obs[1, ])
  expect_error(prep(o4, sim$cells, "STD", hab), "Duplicated")
  ## missing cells are filled with zero counts
  o5 <- sim$obs[sim$obs$count > 0, ]
  d5 <- prep(o5, sim$cells, "STD", hab)
  expect_equal(sum(d5$jags_data$ys), sum(sim$obs$count[sim$obs$dataset == "STD"]))
})

test_that("JAGS engine: the three models run and the outputs are consistent", {
  skip_if_not_installed("rjags")
  sim <- small_sim()
  fit <- suppressWarnings(fun_abund_hab_estim(
    sim$obs, sim$cells, "STD", hab, n.chains = 2, n.adapt = 100,
    n.burnin = 200, n.iter = 400, thin = 2, seed = 1))
  expect_s3_class(fit, "vgi4bio_bayes")
  expect_output(print(fit), "habitat")
  expect_equal(coda::nchain(fit$samples), 2)

  ra <- fun_rel_abund(fit, sites = sim$sites)
  expect_equal(nrow(ra), 4 * 6)
  expect_equal(as.numeric(tapply(ra$rel_abund, ra$species, sum)), rep(1, 4))
  expect_true(all(ra$lower <= ra$rel_abund & ra$rel_abund <= ra$upper))
  expect_error(fun_rel_abund(fit), "sites")
  expect_error(fun_rel_abund(fit, scale = "none", sites = sim$sites), "sum")

  hp <- fun_hab_pref(fit)
  expect_equal(nrow(hp), 4 * 2)
  expect_equal(hp$pref[hp$habitat == "hab1"], rep(1, 4))
  hs <- fun_hab_pref(fit, scale = "sum")
  expect_equal(as.numeric(tapply(hs$pref, hs$species, sum)), rep(1, 4), tolerance = 1e-8)
  expect_equal(nrow(fun_hab_pref(fit, who = "observers")), 2)

  pr <- fun_pred_std(fit, sim$cells[sim$cells$dataset == "STD", ])
  expect_equal(nrow(pr), 4 * 6)
  expect_true(all(pr$pred > 0))

  ## the estimates are correlated with the truth, even with a short run
  mm <- merge(ra, sim$truth$rel_abund, by = c("species", "site"))
  expect_gt(cor(mm$rel_abund.x, mm$rel_abund.y), 0.8)

  fit0 <- suppressWarnings(fun_abund_hab_estim(
    sim$obs, sim$cells, "STD", hab, model = "no_habitat", n.chains = 2,
    n.adapt = 100, n.burnin = 100, n.iter = 200, seed = 1))
  expect_equal(nrow(fun_rel_abund(fit0)), 4 * 6)          # equal site areas
  expect_error(fun_hab_pref(fit0), "no habitat")
  fits <- suppressWarnings(fun_abund_hab_estim(
    sim$obs, sim$cells, "STD", hab, model = "standardized_only", n.chains = 2,
    n.adapt = 100, n.burnin = 100, n.iter = 200, seed = 1))
  expect_null(fits$opportunistic)
  expect_error(fun_hab_pref(fits, who = "observers"), "habitat")
  fite <- suppressWarnings(fun_abund_hab_estim(
    sim$obs, sim$cells, "STD", hab, effort = "estimated", n.chains = 2,
    n.adapt = 100, n.burnin = 100, n.iter = 200, seed = 1))
  expect_true(any(grepl("^E2", coda::varnames(fite$samples))))
  expect_error(fun_abund_hab_estim(sim$obs, sim$cells, "STD", hab,
                                   effort = "estimated", engine = "nimble_hmc"),
               "marginalized")
})

test_that("nimble HMC engine converges and recovers the simulated values", {
  ## long test (C++ compilation + MCMC: several minutes)
  skip_if_not(identical(Sys.getenv("VGI4BIO_LONG_TESTS"), "true"),
              "long test: set VGI4BIO_LONG_TESTS=true to run it")
  skip_if_not_installed("nimble", "1.4.0")
  skip_if_not_installed("nimbleHMC")
  sim <- fun_simul_hab(n_species = 6, n_sites = 8, n_std = 6, n_opp = 10, seed = 3)
  fit <- fun_abund_hab_estim(sim$obs, sim$cells, "STD", hab,
                             engine = "nimble_hmc", n.chains = 2,
                             n.burnin = 300, n.iter = 300, thin = 1, seed = 1)
  expect_lt(max(fit$rhat, na.rm = TRUE), 1.1)
  ra <- fun_rel_abund(fit, sites = sim$sites)
  mm <- merge(ra, sim$truth$rel_abund, by = c("species", "site"))
  expect_gt(cor(mm$rel_abund.x, mm$rel_abund.y), 0.9)
  hp <- fun_hab_pref(fit)
  hp <- hp[hp$habitat == "hab2", ]
  expect_gt(cor(hp$pref, sim$truth$S[, 2] / sim$truth$S[, 1]), 0.9)
})
