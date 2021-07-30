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
estmods <- c("survey",     # Survey-only
             "spatial_ab", # All data, spatial abundance
             "spatial_q")  # All data, spatial abundance + catchability

## List all possible combinations of OM/EM in given replicate range
specify_fits <- function(study, repls, opmods, estmods, root_dir = ".") {
  create_res_dir(study, repls)

  res_paths <- all_res_file_paths(study, repls, opmods, estmods, root_dir)

  df <- cross_df(list(estmod = estmods, opmod = opmods, repl = repls)) %>%
    mutate(study = study,
           Rdata = res_paths$rdata,
           sub_df = map(estmod, em_subsample),
           estd = map(estmod, em_estd),
           root_dir = root_dir)
  df
}

## Which still need to be fit?
fits_todo <- function(fit_spec, result_root = "prefintensity/results") {
  fit_spec %>%
    filter(!file.exists(Rdata))
}

## How many observations to use from each ?
## specify_subset <- function(estmod) {
##   sub_df <- switch(estmod,
##                    ## Don't use any fish-dep data for survey index
##                    survey = data.frame(vessel_idx = 2, n = 0),
##                    spatial_ab = data.frame(vessel_idx = 2, n = 4000),
##                    spatial_q = data.frame(vessel_idx = 2, n = 4000))
##   sub_df
## }

## Specify which parameters to estimate for each estimation model; don't
## estimate catchability parameters if using a single survey vessel.
## estmod_pars <- function(estmod) {
##   switch(estmod,
##          survey = specify_estimated(beta = TRUE,
##                                     gamma = FALSE,
##                                     omega = list(omega_n = TRUE,
##                                                  omega_w = FALSE),
##                                     epsilon = FALSE,
##                                     lambda = FALSE, # Survey-only
##                                     eta = FALSE,
##                                     phi = FALSE,
##                                     psi = FALSE,
##                                     kappa_map =
##                                       c(1, NA, NA, NA, NA, NA, NA, NA),
##                                     obs_lik = 1L),
##          spatial_ab = specify_estimated(beta = TRUE,
##                                         gamma = FALSE,
##                                         omega = list(omega_n = TRUE,
##                                                      omega_w = FALSE),
##                                         epsilon = FALSE,
##                                         phi = FALSE,
##                                         psi = FALSE,
##                                         kappa_map =
##                                           c(1, NA, NA, NA, NA, NA, NA, NA),
##                                         obs_lik = 1L),
##          spatial_q = specify_estimated(beta = TRUE,
##                                        gamma = FALSE,
##                                        omega = list(omega_n = TRUE,
##                                                     omega_w = FALSE),
##                                        epsilon = FALSE,
##                                        lambda = TRUE,
##                                        eta = FALSE,
##                                        phi = list(phi_n = TRUE,
##                                                   phi_w = FALSE),
##                                        psi = FALSE,
##                                        kappa_map =
##                                          c(1, NA, NA, NA, 1, NA, NA, NA),
##                                        obs_lik = 1L))
## }

fit_list <- fits_todo(specify_fits(study = study,
                                   repls = repls,
                                   opmods = opmods,
                                   estmods = estmods,
                                   root_dir = "."))

## Iterate over rows to fit each model
for (idx in seq_len(nrow(fit_list))) {
  spec <- spatq_simstudyspec(as.list(fit_list[idx, ]))

  setup <- spatq_simsetup(repl = spec$repl,
                          study,
                          spec$opmod,
                          spec$sub_df[[1]],
                          max_T = max_T,
                          root_dir = root_dir,
                          index_step = 1,
                          spec_estd = spec$estd[[1]])
  obj <- spatq_obj(setup,
                   runSymbolicAnalysis = TRUE,
                   normalize = TRUE,
                   silent = TRUE)

  fit <- tryCatch({
    ## Fit with large number of iterations and do it twice so more likely to
    ## reach optimum. Previous fits have ended early and/or with large
    ## gradient components, and many of these did not have PD Hessians
    fit <- spatq_fit(obj = obj, control = optcontrol)
    fit <- spatq_fit(obj = obj, fit = fit, control = optcontrol)
    fit},
    error = function(e) list(fail = TRUE))
  lpb <- tryCatch(
    gather_nvec(obj$env$last.par.best),
    error = function(e) list(fail = TRUE))
  rep <- tryCatch(
    report_spatq(obj),
    error = function(e) list(fail = TRUE))
  sdr <- tryCatch(
    sdreport_spatq(obj),
    error = function(e) list(fail = TRUE))

  save_fit(spec, fit, lpb, rep, sdr)
  save_index(spec, sdr, feather = FALSE)
}
