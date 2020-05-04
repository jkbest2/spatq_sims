library(spatq)
library(tidyverse)

specify_fits <- function(repl = 1:100,
                         data = c("all", "indep"),
                         opmod = "combo",
                         root = "results") {
  if (!file.exists(root)) dir.create(root)
  df <- cross_df(list(data = data, opmod = opmod, repl = repl)) %>%
    mutate(repl_str = str_pad(repl, 2, pad = 0),
           Rdata = file.path(root, paste0(repl_str, "_", opmod, "_",
                                          data, ".Rdata")),
           index = file.path(root, paste0(repl_str, "_", opmod,
                                          "_", data, "_index.csv")),
           sub_df = map(data, specify_subset),
           map_pars = map(data, specify_map_pars))
  df
}

fits_todo <- function(fit_spec, result_root = "results") {
  fit_spec %>%
    filter(!file.exists(Rdata))
}

specify_subset <- function(estmod) {
  sub_df <- switch(estmod,
                   all = NULL,
                   indep = data.frame(vessel_idx = 2, n = 0))
  sub_df
}

specify_map_pars <- function(estmod) {
  switch(estmod,
         all = c(
           "gamma_n", "gamma_w"
           ## ,"omega_n", "omega_w"
          , "epsilon1_n", "epsilon1_w"
          , "eta_n", "eta_w"
           ## ,"phi_n", "phi_w"
          , "psi1_n", "psi1_w"
           ## ,"log_tau"
          , "log_kappa"
         ),
         indep = c(
           "gamma_n", "gamma_w"
           ## , "omega_n", "omega_w"
          , "epsilon1_n", "epsilon1_w"
          , "lambda_n", "lambda_w"
          , "eta_n", "eta_w"
          , "phi_n", "phi_w"
          , "psi1_n", "psi1_w"
           ## , "log_tau"
         , "log_kappa"
         ))
}

main <- function(max_T = 15) {
  fit_list <- fits_todo(specify_fits())

  for (idx in seq_len(nrow(fit_list))) {
    spec <- fit_list[idx, ]
    obj <- make_sim_adfun(repl = spec$repl,
                          sc = spec$opmod,
                          sub_df = spec$sub_df[[1]],
                          root_dir = ".",
                          max_T = max_T,
                          map_pars = spec$map_pars[[1]],
                          runSymbolicAnalysis = TRUE,
                          silent = TRUE,
                          normalize = FALSE)

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
                        estmod = spec$data,
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
