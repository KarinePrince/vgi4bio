## small subset of the Aquitaine data, to keep tests fast
make_df <- function() {
  sp <- c("ALAARV", "COLPAL", "TURMER", "STRDEC")
  st <- levels(Aquitaine$Site)[1:15]
  aq <- droplevels(Aquitaine[Aquitaine$Species %in% sp & Aquitaine$Site %in% st, ])
  rbind(
    data.frame(species = aq$Species, site = aq$Site, N_obs = aq$X_ACT,
               effort = aq$Effort_ACT, dataset = "ACT"),
    data.frame(species = aq$Species, site = aq$Site, N_obs = aq$X_LPO,
               effort = 1, dataset = "LPO"))
}

test_that("Aquitaine data set has the documented structure", {
  expect_s3_class(Aquitaine, "data.frame")
  expect_equal(dim(Aquitaine), c(2142L, 12L))
  expect_equal(nlevels(Aquitaine$Species), 34L)
  expect_equal(nlevels(Aquitaine$Site), 63L)
})

test_that("fun_abund_estim fits the model and returns the expected objects", {
  df <- make_df()
  res <- fun_abund_estim(df, standardized = "ACT")
  expect_named(res, c("mod", "N_pred", "standardized"))
  expect_equal(res$standardized, "ACT")
  expect_s3_class(res$mod, "glm")
  expect_true(res$mod$converged)
  expect_null(res$N_pred)
  expect_equal(res$mod$family$family, "quasipoisson")
  expect_equal(fun_abund_estim(df, "ACT", family = "poisson")$mod$family$family, "poisson")
})

test_that("predictions are returned on the requested scale", {
  df <- make_df()
  nd <- transform(df[df$dataset == "ACT", ], effort = 1)
  link <- fun_abund_estim(df, "ACT", newdata = nd)$N_pred
  resp <- fun_abund_estim(df, "ACT", newdata = nd, type = "response")$N_pred
  expect_length(resp, nrow(nd))
  expect_true(all(resp >= 0))
  expect_equal(unname(exp(link)), unname(resp), tolerance = 1e-8)
})

test_that("fun_abund_estim validates its inputs", {
  df <- make_df()
  expect_error(fun_abund_estim(df[, setdiff(names(df), "effort")], "ACT"), "effort")
  df_na <- df; df_na$N_obs[1] <- NA
  expect_error(fun_abund_estim(df_na, "ACT"), "NA")
  df_0 <- df; df_0$effort[1] <- 0
  expect_error(fun_abund_estim(df_0, "ACT"), "strictly positive")
  expect_error(fun_abund_estim(df, "ACT", family = "binomial"))
  expect_error(fun_abund_estim(df, "ACT", newdata = df[, c("species", "site")]), "newdata")
  expect_error(fun_abund_estim(df, "STOC"), "not a value")
  expect_error(fun_abund_estim(df, c("ACT", "LPO")), "Only one")
  expect_warning(fun_abund_estim(df), "assuming 'ACT'")
})

test_that("a warning is issued for species with only zeros in a dataset", {
  df <- make_df()
  df$N_obs[df$species == "ALAARV" & df$dataset == "ACT"] <- 0
  expect_warning(fun_abund_estim(df, "ACT"), "ALAARV in ACT")
})

test_that("the standardized dataset does not depend on dataset names", {
  df <- make_df()
  df2 <- df
  df2$dataset <- ifelse(df2$dataset == "ACT", "zACT", "LPO")
  ra1 <- fun_rel_abund(fun_abund_estim(df, "ACT"))
  ra2 <- fun_rel_abund(fun_abund_estim(df2, "zACT"))
  expect_equal(ra1$rel_abund, ra2$rel_abund, tolerance = 1e-6)
})

test_that("fun_pearsons and fun_braycurt return quartiles", {
  nd <- data.frame(species = rep(c("a", "b", "c"), each = 4))
  pred <- c(1, 2, 3, 4,  2, 2, 5, 1,  1, 1, 1, 2)
  # identical predictions and references: perfect agreement
  expect_equal(unname(fun_pearsons(pred, pred, nd)), c(1, 1, 1))
  expect_equal(unname(fun_braycurt(pred, pred, nd)), c(0, 0, 0))
  # proportional predictions: Bray-Curtis on relative abundances is 0
  expect_equal(unname(fun_braycurt(pred, 3 * pred, nd)), c(0, 0, 0))
  expect_named(fun_pearsons(pred, rev(pred), nd), c("25%", "50%", "75%"))
  expect_error(fun_pearsons(pred, pred), "newdata")
  expect_error(fun_braycurt(pred, pred[-1], nd), "same length")
})
