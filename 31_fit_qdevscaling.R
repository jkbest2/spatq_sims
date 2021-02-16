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

## What range of replicates are we looking for?
repls <- as.numeric(commandArgs(trailingOnly = TRUE))

study <- "qdevscaling"

## List all possible combinations of OM/EM in given replicate range
specify_fits <- function(repl = repls[1]:repls[2],
                         estmod = c("survey",
                                    "spatial_ab",
                                    "spatial_q"),
                         opmod = factor(1:5),
                         root = "qdevscaling/results") {
  if (!file.exists(root)) dir.create(root)

  df <- cross_df(list(estmod = estmod, opmod = opmod, repl = repl)) %>%
    mutate(Rdata = get_rdata_filename(repl, estmod, opmod, root),
           index = get_indexcsv_filename(repl, estmod, opmod, root),
           sub_df = map(estmod, specify_subset),
           estd = map(estmod, estmod_pars))
  df
}

get_repl_str <- function(repl) {
  str_pad(repl, 2, pad = 0)
}

get_rdata_filename <- function(repl, estmod, opmod, root = "qdevscaling/results") {
  repl_str <- get_repl_str(repl)
  file.path(root, paste0(repl_str, "_", opmod, "_", estmod, ".Rdata"))
}

get_indexcsv_filename <- function(repl, estmod, opmod, root = "qdevscaling/results") {
  repl_str <- get_repl_str(repl)
  file.path(root, paste0(repl_str, "_", opmod, "_", estmod, "_index.csv"))
}

## Which still need to be fit?
fits_todo <- function(fit_spec, result_root = "qdevscaling/results") {
  fit_spec %>%
    filter(!file.exists(Rdata))
}

## How many observations to use from each ?
specify_subset <- function(estmod) {
  sub_df <- switch(estmod,
                   ## Don't use any fish-dep data for survey index
                   survey = data.frame(vessel_idx = 2, n = 0),
                   spatial_ab = NULL,
                   spatial_q = NULL)
  sub_df
}

## Specify which parameters to estimate for each estimation model; don't
## estimate catchability parameters if using a single survey vessel.
estmod_pars <- function(estmod) {
  switch(estmod,
         survey = specify_estimated(beta = TRUE,
                                    gamma = FALSE,
                                    omega = list(omega_n = TRUE,
                                                 omega_w = FALSE),
                                    epsilon = FALSE,
                                    lambda = FALSE, # Survey-only
                                    eta = FALSE,
                                    phi = FALSE,
                                    psi = FALSE,
                                    kappa_map =
                                      c(1, NA, NA, NA, NA, NA, NA, NA),
                                    obs_lik = 1L),
         spatial_ab = specify_estimated(beta = TRUE,
                                        gamma = FALSE,
                                        omega = list(omega_n = TRUE,
                                                     omega_w = FALSE),
                                        epsilon = FALSE,
                                        phi = FALSE,
                                        psi = FALSE,
                                        kappa_map =
                                          c(1, NA, NA, NA, NA, NA, NA, NA),
                                        obs_lik = 1L),
         spatial_q = specify_estimated(beta = TRUE,
                                       gamma = FALSE,
                                       omega = list(omega_n = TRUE, omega_w = FALSE),
                                       epsilon = FALSE,
                                       lambda = TRUE,
                                       eta = FALSE,
                                       phi = list(phi_n = TRUE, phi_w = FALSE),
                                       psi = FALSE,
                                       kappa_map =
                                         c(1, NA, NA, NA, 2, NA, NA, NA),
                                       obs_lik = 1L))
}

main <- function(max_T = 15,
                 optcontrol = list(eval.max = 1000L, iter.max = 750L)) {
  fit_list <- fits_todo(specify_fits())

  for (idx in seq_len(nrow(fit_list))) {
    spec <- fit_list[idx, ]

    setup <- spatq_simsetup(repl = spec$repl,
                            study,
                            spec$opmod,
                            spec$sub_df[[1]],
                            max_T = max_T,
                            root_dir = "qdevscaling",
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

    saveRDS(list(spec = spec, fit = fit, lpb = lpb, rep = rep, sdr = sdr),
            spec$Rdata)

    ## Read true population state and calculate index
    true_index <- read_popstate(spec$repl, spec$opmod) %>%
      rename(year = time,
             raw_true = pop) %>%
      filter(year <= max_T) %>%
      mutate(index_true = rescale_index(raw_true)$index)

    if (!("fail" %in% names(sdr))) {
      ## Organize details for estimated index
      which_index <- which(names(sdr$value) == "Index")
      est_index <- tibble(repl = spec$repl,
                          opmod = spec$opmod,
                          estmod = spec$estmod,
                          year = 1:max_T,
                          raw_est = sdr$value[which_index],
                          index_est = rescale_index(raw_est)$index,
                          raw_unb = sdr$unbiased$value[which_index],
                          index_unb = rescale_index(raw_unb)$index,
                          raw_sd = sdr$sd[which_index],
                          index_sd = rescale_index(raw_est, sd_raw)$sd,
                          raw_unb_sd = sdr$unbiased$sd,
                          unb_sd = rescale_index(raw_unb, raw_unb_sd)$sd)
    } else {
      est_index <- tibble(repl = spec$repl,
                          opmod = spec$opmod,
                          estmod = spec$estmod,
                          year = 1:max_T,
                          raw_est = rep(NA, max_T),
                          index_est = rep(NA, max_T),
                          raw_unb = rep(NA, max_T),
                          index_unb = rep(NA, max_T),
                          raw_sd = rep(NA, max_T),
                          index_sd = rep(NA, max_T),
                          raw_unb_sd = rep(NA, max_T),
                          unb_sd = rep(NA, max_T))

    }
    ## Join and write to CSV file
    index_df <- left_join(est_index, true_index, by = "year")
    write_csv(index_df, spec$index)
  }
}

main()
