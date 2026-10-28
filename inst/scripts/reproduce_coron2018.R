## ===========================================================================
## Reproduction of Table 1 of Coron et al. (2018) with vgi4bio
##
## Coron C., Calenge C., Giraud C. & Julliard R. (2018). Bayesian estimation
## of species relative abundances and habitat preferences using opportunistic
## data. Environmental and Ecological Statistics, 25(1), 71-93.
##
## Input files (from the authors' material, not distributed with vgi4bio):
##   * datasetsFormer.RData       (objects act, lpo, stoc, reductACT)
##   * VolHabACTLPOSTOC7.csv      (area of the 7 habitat types in each ACT
##                                 point, LPO commune and STOC point)
##
## Usage, from a terminal (one model per R session lets models run in
## parallel on a multi-core machine):
##   Rscript reproduce_coron2018.R <data_dir> <model> <engine> [quick]
##     <model>  : habitat | no_habitat | standardized_only | table
##                ("table" gathers the saved fits and prints Table 1)
##     <engine> : nimble_hmc (recommended) | jags
##     quick    : very short MCMC, only to check that everything runs
##
## Computing time: a converged fit takes many hours (one NUTS iteration takes
## ~20 s on these data, one JAGS iteration ~4 s but JAGS needs tens of
## thousands of iterations). Fits are saved in <data_dir>/vgi4bio_fits/.
## ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript reproduce_coron2018.R <data_dir> <model> [engine] [quick]")
data_dir <- args[1]
what <- args[2]
engine <- if (length(args) >= 3) args[3] else "nimble_hmc"
quick <- length(args) >= 4 && args[4] == "quick"

library(vgi4bio)
out_dir <- file.path(data_dir, "vgi4bio_fits")
dir.create(out_dir, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. Data preparation (same choices as the authors' RealData.Rnw)
## ---------------------------------------------------------------------------
prep_coron <- function(data_dir) {
  e <- new.env()
  load(file.path(data_dir, "datasetsFormer.RData"), envir = e)
  act <- e$act
  lpo <- if (methods::is(e$lpo, "Spatial")) e$lpo@data else e$lpo
  stoc <- if (methods::is(e$stoc, "Spatial")) e$stoc@data else e$stoc
  V <- utils::read.csv(file.path(data_dir, "VolHabACTLPOSTOC7.csv"))
  hab <- c("urban", "intensive_agri", "open_natural", "farmland",
           "mixed_forest", "deciduous_forest", "coniferous_forest")
  names(V) <- c(hab, "cell")
  V$cell <- as.character(V$cell)

  ## ACT (standardized): counts per listening point; effort = number of years
  a <- stats::aggregate(Nombre.contact.auditif ~ esp + pt, data = act, sum)
  obs_act <- data.frame(species = a$esp, cell = a$pt, dataset = "ACT",
                        count = a$Nombre.contact.auditif)
  eff <- tapply(act$Annee, as.character(act$pt), function(x) length(unique(x)))
  qa <- tapply(as.character(act$Code.carte.IGN), as.character(act$pt), `[`, 1)
  cells_act <- data.frame(cell = names(eff), dataset = "ACT",
                          site = qa[names(eff)], effort = as.numeric(eff))

  ## LPO (opportunistic): counts per commune; quadrat of the commune = quadrat
  ## of its first record, as in RealData.Rnw
  l <- stats::aggregate(eff ~ code + insee, data = lpo, sum)
  obs_lpo <- data.frame(species = l$code, cell = as.character(l$insee),
                        dataset = "LPO", count = l$eff)
  ql <- tapply(as.character(lpo$quad), as.character(lpo$insee), `[`, 1)
  cells_lpo <- data.frame(cell = names(ql), dataset = "LPO", site = ql,
                          effort = NA)

  cells <- rbind(cells_act, cells_lpo)
  cells <- cbind(cells, V[match(cells$cell, V$cell), hab])
  if (anyNA(cells[, hab])) stop("Some ACT/LPO cells have no habitat data.")

  ## STOC (independent validation data): counts per point and quadrat
  s <- stats::aggregate(eff ~ esp + carrept, data = stoc, sum)
  stoc_obs <- data.frame(species = s$esp, cell = s$carrept, count = s$eff)
  ## (as.character: unused factor levels would give NA efforts)
  es <- tapply(stoc$annee, as.character(stoc$carrept), function(x) length(unique(x)))
  qs <- tapply(as.character(stoc$quad), as.character(stoc$carrept), `[`, 1)
  stoc_cells <- data.frame(cell = names(es), site = qs[names(es)],
                           effort = as.numeric(es))
  stoc_cells <- cbind(stoc_cells, V[match(stoc_cells$cell, V$cell), hab])
  if (anyNA(stoc_cells[, hab])) {
    message(sum(!stats::complete.cases(stoc_cells[, hab])),
            " STOC points without habitat data are dropped.")
    stoc_cells <- stoc_cells[stats::complete.cases(stoc_cells[, hab]), ]
  }
  list(obs = rbind(obs_act, obs_lpo), cells = cells, stoc_obs = stoc_obs,
       stoc_cells = stoc_cells, habitats = hab)
}

d <- prep_coron(data_dir)

## ---------------------------------------------------------------------------
## 2. Model fit (one model per call)
## ---------------------------------------------------------------------------
mcmc <- if (engine == "nimble_hmc") {
  if (quick) list(n.chains = 1, n.burnin = 10, n.iter = 10, thin = 1)
  else list(n.chains = 3, n.burnin = 1000, n.iter = 1000, thin = 1)
} else {
  if (quick) list(n.chains = 1, n.adapt = 10, n.burnin = 0, n.iter = 10, thin = 1)
  else list(n.chains = 3, n.adapt = 1000, n.burnin = 20000, n.iter = 60000, thin = 30)
}

if (what %in% c("habitat", "no_habitat", "standardized_only")) {
  t0 <- Sys.time()
  fit <- do.call(fun_abund_hab_estim,
                 c(list(obs = d$obs, cells = d$cells, standardized = "ACT",
                        habitats = d$habitats, model = what, engine = engine,
                        seed = 1, quiet = FALSE), mcmc))
  fit$time <- difftime(Sys.time(), t0, units = "hours")
  print(fit)
  cat("Computing time:", format(fit$time), "\n")
  saveRDS(fit, file.path(out_dir, paste0("fit_", what, "_", engine,
                                         if (quick) "_quick", ".rds")))
}

## ---------------------------------------------------------------------------
## 3. Table 1: per-species Pearson correlation between predicted and observed
##    STOC counts per quadrat (median and quartiles)
## ---------------------------------------------------------------------------
if (what == "table") {
  ## observed STOC counts per unit of effort, per species and quadrat
  st <- merge(d$stoc_obs, d$stoc_cells[, c("cell", "site", "effort")], by = "cell")
  obs_q <- stats::aggregate(count ~ species + site, data = st, sum)
  eff_q <- tapply(d$stoc_cells$effort, d$stoc_cells$site, sum)
  obs_q$obs <- obs_q$count / eff_q[obs_q$site]

  published <- list(
    habitat = list(ACT = c(0.27, 0.49, 0.57), nonACT = c(0.13, 0.35, 0.52)),
    standardized_only = list(ACT = c(0.08, 0.28, 0.47), nonACT = NULL),
    no_habitat = list(ACT = c(0.31, 0.44, 0.68), nonACT = c(0.17, 0.39, 0.47)))

  res <- NULL
  for (m in names(published)) {
    f <- file.path(out_dir, paste0("fit_", m, "_", engine,
                                   if (quick) "_quick", ".rds"))
    if (!file.exists(f)) { message("No saved fit for ", m, ": skipped."); next }
    fit <- readRDS(f)
    pr <- suppressMessages(fun_pred_std(fit, d$stoc_cells))
    mm <- merge(pr, obs_q, by = c("species", "site"))
    for (g in c("ACT", "nonACT")) {
      if (is.null(published[[m]][[g]])) next
      sub <- mm[mm$in_standardized == (g == "ACT"), ]
      q <- fun_pearsons(sub$pred, sub$obs, sub)
      res <- rbind(res, data.frame(
        model = m, species = g, n_species = length(unique(sub$species)),
        max_rhat = if (all(is.na(fit$rhat))) NA else round(max(fit$rhat, na.rm = TRUE), 2),
        vgi4bio = sprintf("%.2f (%.2f-%.2f)", q[2], q[1], q[3]),
        published = do.call(sprintf, c("%.2f (%.2f-%.2f)",
                                       as.list(published[[m]][[g]][c(2, 1, 3)])))))
    }
  }
  print(res, row.names = FALSE)
  utils::write.csv(res, file.path(out_dir, paste0("table1_", engine,
                                                  if (quick) "_quick", ".csv")),
                   row.names = FALSE)
  cat("\nPublished values: Coron et al. (2018), Table 1 / supplementary file",
      "RealData (1 chain of 500 JAGS iterations).\n")
}
