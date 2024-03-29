## library(spatq)
devtools::load_all("~/src/spatq", helpers = FALSE)
library(tidyverse)

specify_fits <- function(repl = 1:50,
                         estmod = c("spatq", "all", "indep"),
                         opmod = "combo",
                         root = "~/spatq_simruns/results") {
  if (!file.exists(root)) dir.create(root)
  df <- cross_df(list(estmod = estmod, opmod = opmod, repl = repl)) %>%
    mutate(repl_str = str_pad(repl, 2, pad = 0),
           Rdata = file.path(root, paste0(repl_str, "_", opmod, "_",
                                          estmod, ".Rdata")),
           index = file.path(root, paste0(repl_str, "_", opmod,
                                          "_", estmod, "_index.csv")),
           sub_df = map(estmod, specify_subset),
           estd = map(estmod, estmod_pars))
  df
}

fits_todo <- function(fit_spec, result_root = "results") {
  fit_spec %>%
    filter(!file.exists(Rdata))
}

specify_subset <- function(estmod) {
  sub_df <- switch(estmod,
                   spatq = NULL,
                   all = NULL,
                   indep = data.frame(vessel_idx = 2, n = 0))
  sub_df
}

## Specify which parameters to estimate for each estimation model; don't
## estimate catchability parameters if using a single survey vessel.
estmod_pars <- function(estmod) {
  switch(estmod,
         spatq = specify_estimated(beta = TRUE,
                                   gamma = FALSE,
                                   omega = TRUE,
                                   epsilon = TRUE,
                                   lambda = TRUE,
                                   eta = FALSE,
                                   phi = TRUE,
                                   psi = FALSE),
         all = specify_estimated(beta = TRUE,
                                 gamma = FALSE,
                                 omega = TRUE,
                                 epsilon = TRUE,
                                 lambda = TRUE,
                                 eta = FALSE,
                                 phi = FALSE,
                                 psi = FALSE),
         indep = specify_estimated(beta = TRUE,
                                   gamma = FALSE,
                                   omega = TRUE,
                                   epsilon = TRUE,
                                   lambda = FALSE,
                                   eta = FALSE,
                                   phi = FALSE,
                                   psi = FALSE))
}

main <- function(max_T = 15) {
  fit_list <- fits_todo(specify_fits(opmod = "pref", root = "~/spatq_simruns"),
                        result_root = "~/spatq_simruns")

  for (idx in seq_len(nrow(fit_list))) {
    spec <- fit_list[idx, ]
    obj <- make_sim_adfun(repl = spec$repl,
                          sc = spec$opmod,
                          sub_df = spec$sub_df[[1]],
                          root_dir = "~/gscratch/spatq_sims",
                          max_T = max_T,
                          spec_estd = spec$estd[[1]],
                          runSymbolicAnalysis = TRUE,
                          silent = TRUE,
                          normalize = TRUE)

    fit <- fit_spatq(obj)
    fit <- fit_spatq(obj, fit)
    lpb <- gather_nvec(obj$env$last.par.best)
    rep <- report_spatq(obj)
    sdr <- sdreport_spatq(obj)

    saveRDS(list(spec = spec, fit = fit, lpb = lpb, rep = rep, sdr = sdr),
            spec$Rdata)

    ## Read true population state and calculate index
    true_index <- read_popstate(spec$repl, spec$opmod) %>%
      rename(year = time,
             raw_true = pop) %>%
      filter(year <= max_T) %>%
      mutate(index_true = rescale_index(raw_true))

    ## Organize details for estimated index
    which_index <- which(names(sdr$value) == "Index")
    est_index <- tibble(repl = spec$repl,
                        opmod = spec$opmod,
                        estmod = spec$estmod,
                        year = 1:max_T,
                        raw_est = sdr$value[which_index],
                        index_est = rescale_index(raw_est),
                        sd_raw = sdr$sd[which_index],
                        sd = sd_raw / attr(index_est, "scale"))

    ## Join and write to CSV file
    index_df <- left_join(est_index, true_index, by = "year")
    write_csv(index_df, spec$index)
  }
}

main()
