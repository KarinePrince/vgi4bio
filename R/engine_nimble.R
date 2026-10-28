## ---------------------------------------------------------------------------
## nimble HMC engine for fun_abund_hab_estim(engine = "nimble_hmc")
##
## Same model, constraints and priors as the JAGS engine (opportunistic
## efforts marginalized). Only the MCMC algorithm differs:
## * the positive parameters (N, S, P, q) are sampled on the log scale. The
##   prior of log(X) when X ~ Gamma(a, b) is given by the custom distribution
##   dvgi4bioLogGamma below, so the priors are exactly those of the authors;
## * all parameters are updated jointly by the No-U-Turn sampler (NUTS) of
##   nimbleHMC.
## On simulated data, NUTS reached convergence (max Rhat 1.05) with 100
## warm-up + 200 iterations per chain, where JAGS needed about 80000
## iterations (about 4 times longer in wall-clock time). A previous attempt
## with nimble's block random walk samplers (per species) did not improve on
## JAGS, because the slowest direction is global (q vs all the S[i, ]); it
## was removed.
## ---------------------------------------------------------------------------

.nimble_code <- function(model, nP) {
  hab <- model != "no_habitat"
  opp <- model != "standardized_only"
  txt <- "{
  for (i in 1:I) { for (j in 1:J) {
    lN[i, j] ~ dvgi4bioLogGamma(1.0E-11, 1.0E-11)
    N[i, j] <- exp(lN[i, j])
  } }"
  if (hab) txt <- paste0(txt, "
  for (i in 1:I) {
    S[i, 1] <- 1
    for (h in 2:H) {
      lS[i, h - 1] ~ dvgi4bioLogGamma(1.0E-4, 1.0E-4)
      S[i, h] <- exp(lS[i, h - 1])
    }
  }")
  if (opp) {
    if (nP > 0) txt <- paste0(txt, "
  for (m in 1:nP) {
    lP[m] ~ dvgi4bioLogGamma(0.01, 0.01)
    P1[m] <- exp(lP[m])
  }")
    txt <- paste0(txt, "
  P1[nP + 1] <- 1
  for (i in 1:I) { P[i] <- P1[pidx[i]] }")
  }
  if (opp && hab) txt <- paste0(txt, "
  qo[1] <- 1
  for (h in 2:H) {
    lq[h - 1] ~ dvgi4bioLogGamma(1.0E-4, 1.0E-4)
    q2[h - 1] <- exp(lq[h - 1])
    qo[h] <- q2[h - 1]
  }")
  pstd <- if (opp) " * P[sps[l]]" else ""
  txt <- paste0(txt, "
  for (l in 1:Ls) {
    ys[l] ~ dpois(N[sps[l], sts[l]] * Erel[cs[l]]", pstd,
    if (hab) " * inprod(S[sps[l], 1:H], props[cs[l], 1:H])", ")
  }")
  if (opp) {
    mu <- if (hab) "N[i, sto[c]] * inprod(S[i, 1:H], qo[1:H] * propo[c, 1:H])"
          else "N[i, sto[c]]"
    txt <- paste0(txt, "
  for (c in 1:Co) {
    for (i in 1:I) { mu[c, i] <- ", mu, " }
    pr[c, 1:I] <- mu[c, 1:I] / sum(mu[c, 1:I])
    Yo[c, 1:I] ~ dmulti(pr[c, 1:I], To[c])
  }")
  }
  parse(text = paste0(txt, "\n}"))[[1]]
}


.fit_nimble <- function(jd, model, n.chains, n.burnin, n.iter, thin, seed,
                        quiet, inits) {
  if (!requireNamespace("nimble", quietly = TRUE) ||
      !requireNamespace("nimbleHMC", quietly = TRUE))
    stop("Packages 'nimble' (>= 1.4.0) and 'nimbleHMC' are required for ",
         "engine = \"nimble_hmc\".", call. = FALSE)
  ## nimble's model builder evaluates some internal code in the search path:
  ## it only works when the package is attached (known nimble limitation)
  if (!"package:nimble" %in% search())
    suppressPackageStartupMessages(attachNamespace("nimble"))
  ## same for nimbleHMC, whose samplers are looked up by name
  if (!"package:nimbleHMC" %in% search())
    suppressPackageStartupMessages(attachNamespace("nimbleHMC"))

  ## log of a Gamma(shape, rate) variable. nimble's compiler only finds
  ## user-defined distributions in the global environment: they are put
  ## there for the duration of the fit, under names unlikely to clash, and
  ## removed afterwards.
  ## (nimble DSL code, parsed from text so that R CMD check does not analyse
  ## it as R code)
  dfun <- eval(parse(text = "nimble::nimbleFunction(
    run = function(x = double(0), shape = double(0), rate = double(0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      lp <- shape * log(rate) - lgamma(shape) + shape * x - rate * exp(x)
      if (log) return(lp) else return(exp(lp))
    }, buildDerivs = 'run')"))
  rfun <- eval(parse(text = "nimble::nimbleFunction(
    run = function(n = integer(0), shape = double(0), rate = double(0)) {
      returnType(double(0))
      return(log(rgamma(1, shape, rate)))
    })"))
  genv <- globalenv()
  nm <- c("dvgi4bioLogGamma", "rvgi4bioLogGamma")
  old <- mget(nm, envir = genv, ifnotfound = list(NULL, NULL))
  assign(nm[1], dfun, envir = genv)
  assign(nm[2], rfun, envir = genv)
  on.exit({
    for (k in 1:2) {
      if (is.null(old[[k]])) rm(list = nm[k], envir = genv)
      else assign(nm[k], old[[k]], envir = genv)
    }
  }, add = TRUE)
  ## (re-)registering an already registered distribution only gives a message
  suppressMessages(nimble::registerDistributions(
      list(dvgi4bioLogGamma = list(BUGSdist = "dvgi4bioLogGamma(shape, rate)",
                                   types = c("value = double(0)", "shape = double(0)",
                                             "rate = double(0)"))),
      verbose = FALSE))

  nP <- if (is.null(jd$nP)) 0 else jd$nP
  code <- .nimble_code(model, nP)
  data_names <- intersect(c("ys", "Yo"), names(jd))
  data <- jd[data_names]
  const <- jd[setdiff(names(jd), data_names)]

  ## initial values on the log scale
  ini <- lapply(inits, function(x) {
    out <- list(lN = log(x$N))
    if (!is.null(x$S)) out$lS <- matrix(0, jd$I, jd$H - 1)
    if (!is.null(x$P1)) out$lP <- rep(0, nP)
    if (!is.null(x$q2)) out$lq <- rep(0, jd$H - 1)
    out
  })

  run <- function() {
    m <- nimble::nimbleModel(code, constants = const, data = data,
                             inits = ini[[1]], check = FALSE, calculate = FALSE,
                             buildDerivs = TRUE)
    monitors <- c("N", if (model != "no_habitat") "S",
                  if (model == "habitat") "q2", if (nP > 0) "P1")
    ## No-U-Turn Hamiltonian Monte Carlo on all the (log-scale) parameters
    ## jointly: it uses the gradient of the log-posterior, which lets it
    ## follow the correlations between parameters (N[i, ] and S[i, ] within a
    ## species, and the global one between the observer preferences q and all
    ## the S[i, ]), along which the univariate updates of JAGS mix slowly.
    conf <- nimble::configureMCMC(m, monitors = monitors, thin = thin,
                                  nodes = NULL, print = FALSE)
    nimbleHMC::addHMC(conf, target = m$getNodeNames(stochOnly = TRUE,
                                                    includeData = FALSE),
                      type = "NUTS", print = FALSE)
    mcmc <- nimble::buildMCMC(conf)
    cm <- nimble::compileNimble(m)
    cmcmc <- nimble::compileNimble(mcmc, project = m)
    nimble::runMCMC(cmcmc, niter = n.burnin + n.iter, nburnin = n.burnin,
                    nchains = n.chains, inits = ini,
                    setSeed = if (is.null(seed)) FALSE else seed + seq_len(n.chains),
                    progressBar = !quiet, samplesAsCodaMCMC = TRUE)
  }
  samples <- if (quiet) suppressMessages(run()) else run()
  if (n.chains == 1) samples <- coda::mcmc.list(samples)
  ## node names as in the JAGS engine ("N[1,2]" rather than "N[1, 2]")
  for (k in seq_along(samples))
    coda::varnames(samples[[k]]) <- gsub(" ", "", coda::varnames(samples[[k]]))
  samples
}
