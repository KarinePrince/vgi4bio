## Reproduces Table 2 of Giraud et al. (2016, Biometrics 72:649-658):
## median (IQR) of the per-species Pearson correlation between the relative
## abundances estimated by the ACT + LPO model and the STOC reference, on the
## 24 sites surveyed by STOC.
## The full model takes a few minutes to fit: run with
##   Sys.setenv(VGI4BIO_LONG_TESTS = "true"); devtools::test()

test_that("ACT + LPO model reproduces Giraud et al. (2016), Table 2", {
  skip_if_not(identical(Sys.getenv("VGI4BIO_LONG_TESTS"), "true"),
              "long test: set VGI4BIO_LONG_TESTS=true to run it")

  inACT <- Aquitaine$Species_in_ACT == 1
  df <- rbind(
    data.frame(species = Aquitaine$Species[inACT], site = Aquitaine$Site[inACT],
               N_obs = Aquitaine$X_ACT[inACT], effort = Aquitaine$Effort_ACT[inACT],
               dataset = "ACT"),
    data.frame(species = Aquitaine$Species, site = Aquitaine$Site,
               N_obs = Aquitaine$X_LPO, effort = 1, dataset = "LPO"))
  ra <- fun_rel_abund(fun_abund_estim(df, standardized = "ACT"))

  ref <- Aquitaine[Aquitaine$Site_in_STOC == 1, ]
  ref$N_ref <- ref$X_STOC / ref$Effort_STOC
  m <- merge(ra, ref[, c("Species", "Site", "N_ref")],
             by.x = c("species", "site"), by.y = c("Species", "Site"))

  act <- m[m$in_standardized, ]
  nonact <- m[!m$in_standardized, ]
  expect_equal(length(unique(act$species)), 13L)
  expect_equal(length(unique(act$site)), 24L)

  # published: 0.55 (0.38-0.68) for the 13 ACT species
  expect_equal(unname(round(fun_pearsons(act$rel_abund, act$N_ref, act), 2)),
               c(0.38, 0.55, 0.68))
  # published: 0.35 (0.19-0.47) for the 21 species not monitored by ACT
  expect_equal(unname(round(fun_pearsons(nonact$rel_abund, nonact$N_ref, nonact), 2)),
               c(0.19, 0.35, 0.47))
})
