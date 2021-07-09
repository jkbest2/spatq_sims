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
estmods <- c("design", "model", "survey")
estd_mod <- function(em) {
  switch(em,
         design = specify_estimated(beta = FALSE,
                                    obs_lik = 1L),
         model = specify_estimated(beta = TRUE,
                                   lambda = FALSE,
                                   obs_lik = 1L),
         survey = specify_estimated(beta = TRUE,
                                    omega = TRUE,
                                    lambda = FALSE,
                                    obs_lik = 1L))
}

design_estimator <- function(catch_df, max_T = 15) {
  mat <- matrix(catch_df$catch_biomass, ncol = max_T)
  index <- colMeans(mat)
  names(index) <- rep_along(index, "Index")
  sd <- apply(mat, 2, sd) / sqrt(nrow(mat))

  list(value = index,
       sd = sd,
       unbiased = list(value = index,
                       sd = sd))
}

has_index <- function(spec) {
  file.exists(res_file_paths(spec$study,
                             spec$repl,
                             spec$opmod,
                             spec$estmod,
                             spec$root_dir)$index_feather)
}

spec_df <- cross_df(list(study = study,
                         repl = repls,
                         opmod = opmods,
                         estmod = estmods))

for (i in 1:nrow(spec_df)) {
  spec <- spatq_simstudyspec(as.list(spec_df[i, ]))

  ## Make sure that results directory exists
  rd <- file.path(study_dir(spec$study, root_dir),
                  "results",
                  repl_dir(spec$repl))
  if (!dir.exists(rd)) dir.create(rd)

  if (spec$estmod == "design") {
    catch_df <- read_catch(spec$study,
                           repl = spec$repl,
                           opmod = spec$opmod,
                           root_dir = root_dir,
                           feather = TRUE) %>%
      filter(vessel == "Survey",
             year <= max_T)

    ## Calculate and save design estimator
    dest <- design_estimator(catch_df, max_T)
    if (!dir.exists())
    save_index(spec, dest, max_T)
  } else {
    estd <- estd_mod(spec$estmod)
    setup <- spatq_simsetup(study = spec$study,
                            repl = spec$repl,
                            opmod = spec$opmod,
                            sub_df = data.frame(vessel_idx = 2, n = 0),
                            max_T = max_T,
                            spec_estd = estd)
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
}
