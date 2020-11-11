library(spatq)
## devtools::load_all("~/dev/spatq", helpers = FALSE)
library(tidyverse)

## What range of replicates are we looking for?
repls <- as.numeric(commandArgs(trailingOnly = TRUE))

## List all possible combinations of OM/EM in given replicate range
specify_fits <- function(repl = repls[1]:repls[2],
                         estmod = c("surv",
                                    "spatab",
                                    "spatq"),
                         opmod = c(## "combo",
                                   "pref",
                                   "spat"),
                         root = "results") {
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

get_rdata_filename <- function(repl, estmod, opmod, root = "results") {
  repl_str <- get_repl_str(repl)
  file.path(root, paste0(repl_str, "_", opmod, "_", estmod, ".Rdata"))
}

get_indexcsv_filename <- function(repl, estmod, opmod, root = "results") {
  repl_str <- get_repl_str(repl)
  file.path(root, paste0(repl_str, "_", opmod, "_", estmod, "_index.csv"))
}

## Which still need to be fit?
fits_todo <- function(fit_spec, result_root = "results") {
  fit_spec %>%
    filter(!file.exists(Rdata))
}

## How many observations to use from each ?
specify_subset <- function(estmod) {
  sub_df <- switch(estmod,
                   ## Don't use any fish-dep data for survey index
                   surv = data.frame(vessel_idx = 2, n = 0),
                   spatab = NULL,
                   spatq = NULL)
  sub_df
}

make_setup <- function(repl, estmod, opmod, root = "results") {

}
## Specify which parameters to estimate for each estimation model; don't
## estimate catchability parameters if using a single survey vessel.
estmod_pars <- function(estmod) {
  switch(estmod,
         surv = specify_estimated(beta = TRUE,
                                  gamma = FALSE,
                                  omega = list(omega_n = TRUE, omega_w = FALSE),
                                  epsilon = FALSE,
                                  lambda = FALSE, # Survey-only
                                  eta = FALSE,
                                  phi = FALSE,
                                  psi = FALSE,
                                  kappa_map =
                                    c(1, NA, NA, NA, NA, NA, NA, NA)),
         spatab = specify_estimated(beta = TRUE,
                                    gamma = FALSE,
                                    omega = list(omega_n = TRUE, omega_w = FALSE),
                                    epsilon = FALSE,
                                    phi = FALSE,
                                    psi = FALSE,
                                    kappa_map =
                                      c(1, NA, NA, NA, NA, NA, NA, NA)),
         spatq = specify_estimated(beta = TRUE,
                                   gamma = FALSE,
                                   omega = list(omega_n = TRUE, omega_w = FALSE),
                                   epsilon = FALSE,
                                   lambda = TRUE,
                                   eta = FALSE,
                                   phi = list(phi_n = TRUE, phi_w = FALSE),
                                   psi = FALSE,
                                   kappa_map =
                                     c(1, NA, NA, NA, 2, NA, NA, NA)))
}

get_prev_em <- function(estmod) {
  switch(estmod,
         surv = NULL,
         spatab = NULL,
         spatq = NULL)
         ## spatab = "surv",
         ## spatq = "spatab")
}

get_prev_fit <- function(repl, estmod, opmod, root = "results") {
  prev_em <- get_prev_em(estmod)
  if (is.null(prev_em)) return(NULL)

  fn <- get_rdata_filename(repl, prev_em, opmod, root)
  rd <- readRDS(fn)

  ## Return both fixed and random parameter values
  return(c(rd$sdr$par.fixed, rd$sdr$par.random))
}

main <- function(max_T = 15) {
  fit_list <- fits_todo(specify_fits())

  for (idx in seq_len(nrow(fit_list))) {
    spec <- fit_list[idx, ]

    prev_fit <- get_prev_fit(spec$repl, spec$estmod, spec$opmod, root)
    setup <- spatq_simsetup(repl = spec$repl,
                            spec$opmod,
                            spec$sub_df[[1]],
                            max_T = max_T,
                            index_step = 1,
                            spec_estd = spec$estd[[1]])
    setup <- update_setup(setup, prev_fit, spec$estd[[1]])
    obj <- spatq_obj(setup,
                     runSymbolicAnalysis = TRUE,
                     normalize = TRUE,
                     silent = TRUE)

    optctl <- spatq_optcontrol(maxopts = 4)
    fit <- tryCatch(
      fit_spatq(obj = obj,
                fit = NULL,
                optcontrol = optctl),
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
                          sd_raw = sdr$sd[which_index],
                          sd = rescale_index(raw_est, sd_raw)$sd)
    } else {
      est_index <- tibble(repl = spec$repl,
                          opmod = spec$opmod,
                          estmod = spec$estmod,
                          year = 1:max_T,
                          raw_est = rep(NA, max_T),
                          index_est = rep(NA, max_T),
                          sd_raw = rep(NA, max_T),
                          sd = rep(NA, max_T))
    }
    ## Join and write to CSV file
    index_df <- left_join(est_index, true_index, by = "year")
    write_csv(index_df, spec$index)
  }
}

main()
