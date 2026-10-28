# vgi4bio 0.1.0

* First archived release (Zenodo). Citation files `CITATION.cff` and
  `inst/CITATION` added: see `citation("vgi4bio")`.

## Changes since 0.0.0.9002

* New `fun_abund_hab_estim()`: Bayesian model of Coron et al. (2018),
  estimating relative abundances and habitat preferences of species and
  observers from standardized and opportunistic data described at the cell
  level. Three models (`"habitat"`, `"no_habitat"`, `"standardized_only"`),
  same likelihood, constraints and priors as the authors' JAGS code.
  - `effort = "marginalized"` (default): the unknown efforts of the
    opportunistic cells are integrated out (multinomial counts given the cell
    total), same posterior for the other parameters, far fewer parameters.
  - `engine = "jags"` (default) or `"nimble_hmc"` (No-U-Turn sampler of
    nimbleHMC), which converges with far fewer iterations.
  - Convergence diagnostics (Rhat, effective sample size) with a warning
    when Rhat > 1.1.
* New `fun_hab_pref()` (habitat preferences of species and observers),
  `fun_pred_std()` (predicted counts for an independent standardized survey)
  and `fun_simul_hab()` (simulation following the design of Coron et al.).
* `fun_rel_abund()` also handles the Bayesian models (posterior mean and
  credible intervals).
* New script `inst/scripts/reproduce_coron2018.R` reproducing Table 1 of
  Coron et al. (2018) from the authors' data files.
* New vignette `vignette("habitat")`.
* `Aquitaine` documentation completed from the web appendix of Giraud et al.

# vgi4bio 0.0.0.9002

* New `fun_rel_abund()`: extracts the relative abundance of each species on
  each site from a fitted model, with delta-method confidence intervals,
  either rescaled to sum to one per species (`scale = "sum"`) or per unit of
  effort of the standardized dataset (`scale = "none"`).
* `fun_abund_estim()` gains a `standardized` argument naming the standardized
  dataset (known effort). Previously, the dataset fixed as reference was the
  first one in alphabetical order, which gave wrong estimates when the
  opportunistic dataset came first. Without `standardized`, a warning is now
  issued. The result also contains `standardized`.
* New long test reproducing Table 2 of Giraud et al. (2016)
  (run with `VGI4BIO_LONG_TESTS=true`).
* Vignette: extraction of relative abundances and comparison with STOC.

# vgi4bio 0.0.0.9001

* Package passes `R CMD check`: valid DESCRIPTION, complete documentation,
  documented `Aquitaine` data set, working tests.
* `fun_pearsons()` and `fun_braycurt()` are exported and take `newdata`
  explicitly.
* Warning for species with only zero counts in a dataset (structural zeros).
