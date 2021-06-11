## If script is run in the REPL interactively, use the local `spatq` package and
## manually set a range of replicates to fit, because they were probably not
## passed as command line arguments when R was started. Otherwise (typically on
## Hyak via SLURM), use the installed version of `spatq` and read in the
## replicate numbers from the command line arguments.
if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
  repl_arg <- c(1, 5)
} else {
  library(spatq)
  repl_arg <- as.numeric(commandArgs(trailingOnly = TRUE))
}
library(tidyverse)

## Where are we working?
root_dir <- "."
## Which simulation study are we fitting?
study <- "prefintensity"
## What range of replicates are going to be fit?
repls <- repl_arg[1]:repl_arg[2]
## How many years to fit?
max_T <- 15
## Tune the optimization routine
optcontrol <- list(eval.max = 1000L, iter.max = 750L)
## Names of the operating models
opmods <- 1:6
## Names of the estimation models
estmods <- "model"
estd <- specify_estimated(beta = TRUE,
                          lambda = FALSE,
                          obs_lik = 1L)

spec_df <- cross_df(list(study = study,
                         repl = repls,
                         opmod = opmods,
                         estmod = estmods))

design_estimator <- function(setup, max_T = 15) {
  mat <- matrix(setup$data$catch_obs, ncol = max_T)
  index <- colMeans(mat)
  names(index) <- rep_along(index, "Index")
  sd <- apply(mat, 2, sd) / sqrt(nrow(mat))

  list(value = index,
       sd = sd,
       unbiased = list(value = index,
                       sd = sd))
}

for (i in 1:nrow(spec_df)) {
  spec <- spatq_simstudyspec(as.list(spec_df[i, ]))
  setup <- spatq_simsetup(study = spec$study,
                          repl = spec$repl,
                          opmod = spec$opmod,
                          sub_df = data.frame(vessel_idx = 2, n = 0),
                          max_T = max_T,
                          spec_estd = estd)

  ## Calculate and save design estimator
  dest <- design_estimator(setup, max_T)
  spec0 <- spec
  spec0$estmod <- "design"
  save_index(spec0, dest, max_T)

  ## Now simplest model-based estimator
  obj <- spatq_obj(setup)
  fit <- spatq_fit(obj)
  fit <- spatq_fit(obj, fit = fit, method = "BFGS")
  fit <- spatq_fit(obj, fit = fit)
  sdr <- sdreport_spatq(obj)
  sdr$unbiased <- list(value = sdr$value,
                       sd = sdr$sd)
  rep <- report_spatq(obj)
  lpb <- obj$env$last.par.best

  save_fit(spec, fit, lpb, rep, sdr)
  save_index(spec, sdr, max_T)
}