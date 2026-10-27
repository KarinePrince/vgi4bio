## subset: 4 species monitored by ACT + 1 species recorded only by LPO
make_df2 <- function() {
  sp <- c("ALAARV", "COLPAL", "TURMER", "STRDEC", "AEGCAU")
  st <- levels(Aquitaine$Site)[1:15]
  aq <- droplevels(Aquitaine[Aquitaine$Species %in% sp & Aquitaine$Site %in% st, ])
  a1 <- aq[aq$Species_in_ACT == 1, ]
  rbind(
    data.frame(species = a1$Species, site = a1$Site, N_obs = a1$X_ACT,
               effort = a1$Effort_ACT, dataset = "ACT"),
    data.frame(species = aq$Species, site = aq$Site, N_obs = aq$X_LPO,
               effort = 1, dataset = "LPO"))
}

test_that("fun_rel_abund returns one row per species x site, summing to 1", {
  df <- make_df2()
  fit <- fun_abund_estim(df, standardized = "ACT")
  ra <- fun_rel_abund(fit)
  expect_named(ra, c("species", "site", "rel_abund", "se_log", "lower",
                     "upper", "in_standardized"))
  expect_equal(nrow(ra), 5 * 15)
  expect_equal(as.numeric(tapply(ra$rel_abund, ra$species, sum)), rep(1, 5))
  expect_true(all(ra$lower <= ra$rel_abund & ra$rel_abund <= ra$upper))
  expect_true(all(ra$se_log > 0))
  expect_equal(unique(ra$in_standardized[ra$species == "AEGCAU"]), FALSE)
  expect_equal(unique(ra$in_standardized[ra$species == "TURMER"]), TRUE)
})

test_that("scale = 'none' matches predict() on the standardized dataset", {
  df <- make_df2()
  fit <- fun_abund_estim(df, standardized = "ACT")
  expect_message(rn <- fun_rel_abund(fit, scale = "none"), "AEGCAU")
  expect_true(all(is.na(rn$rel_abund[rn$species == "AEGCAU"])))
  st <- levels(fit$mod$data$site)
  nd <- data.frame(species = "TURMER", site = st, dataset = "ACT", effort = 1)
  pr <- stats::predict(fit$mod, newdata = nd, type = "link", se.fit = TRUE)
  expect_equal(unname(exp(pr$fit)), rn$rel_abund[rn$species == "TURMER"],
               tolerance = 1e-8)
  expect_equal(unname(pr$se.fit), rn$se_log[rn$species == "TURMER"],
               tolerance = 1e-8)
})

test_that("fun_rel_abund options and input checks", {
  df <- make_df2()
  fit <- fun_abund_estim(df, standardized = "ACT")
  ra <- fun_rel_abund(fit, species = "COLPAL")
  expect_equal(unique(as.character(ra$species)), "COLPAL")
  ra90 <- fun_rel_abund(fit, level = 0.90, species = "COLPAL")
  expect_true(all(ra90$upper < ra$upper))
  expect_error(fun_rel_abund(fit$mod), "fun_abund_estim")
  expect_error(fun_rel_abund(fit, species = "XXX"), "Unknown species")
  expect_error(fun_rel_abund(fit, level = 95), "level")
})
