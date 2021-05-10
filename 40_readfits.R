## if script is run in the REPL interactively, use the local `spatq` package and
## manually set a range of replicates to fit, because they were probably not
## passed as command line arguments when R was started. Otherwise (typically on
## Hyak via SLURM), use the installed version of `spatq` and read in the
## replicate numbers from the command line arguments.
if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
  repl_arg <- c(1, 25)
} else {
  library(spatq)
  repl_arg <- as.numeric(commandArgs(trailingOnly = TRUE))
}
library(tidyverse)

## Which simulation study are we fitting?
studies <- c("qdevscaling", "sharedq", "prefintensity")
## What range of replicates are going to be fit?
repls <- repl_arg[1]:repl_arg[2]
## How many years to fit?
max_T <- 15
## Names of the operating models
opmods <- 1:6
## Names of the estimation models
estmods <- c("survey",     # Survey-only
             "spatial_ab", # All data, spatial abundance
             "spatial_q")  # All data, spatial abundance + catchability

results <- lapply(studies,
                  all_res_file_paths,
                  repls = repls,
                  opmods = opmods,
                  estmods = estmods,
                  root_dir = ".")
names(results) <- studies

get_convcode <- function(fitlist) {
  pluck(fitlist, "fit", "convergence")
}

get_mgc <- function(fitlist) {
  pluck(fitlist, "fit", attr_getter("mgc"))
}

get_om_parval <- function(study, opmod = 1:6) {
  studyvals <- switch(study,
                      qdevscaling = 10 ^ seq(-3, -0.5, 0.5),
                      sharedq = seq(0, 1, 0.2),
                      prefintensity = c(0, 1, 2, 4, 8, 16))
  studyvals[opmod]
}

get_om_parlabel <- function(study) {
  switch(study,
         qdevscaling = "log catchability deviation SD",
         sharedq = "Prop shared catchabilty dev",
         prefintensity = "Preference power")
}

## read_all_indices <- function(csv_list) {
##   map_df(csv_list, read_index_csv)
## }

evaluate_bias <- function(index_df) {
  index_df %>%
    group_by(opmod, estmod) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(log(raw_unb) ~ repl + log(raw_true),
                                data = .x)),
           coef = map(mod, coef),
           delta = map_dbl(coef, pluck, "log(raw_true)")) %>%
    select(opmod, estmod, delta) %>%
    mutate(parval = map2_dbl(study, opmod, get_om_parval))
}

plot_bias2 <- function(index_df) {
  bias_df <- evaluate_bias(index_df)
  ggplot(bias_df, aes(x = parval, y = delta, color = estmod)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    labs(x = get_om_parlabel(index_df$study[1]),
         y = "δ bias metric")
}

evaluate_rmse <- function(index_df) {
  index_df %>%
    mutate(sq_err = (index_unb - index_true)^2) %>%
    group_by(opmod, estmod) %>%
    summarize(rmse = sqrt(mean(sq_err)), .groups = "drop") %>%
    select(opmod, estmod, rmse) %>%
    mutate(parval = map2_dbl(study, opmod, get_om_parval))
}

plot_rmse2 <- function(index_df) {
  rmse_df <- evaluate_rmse(index_df)
  ggplot(rmse_df, aes(x = parval, y = rmse, color = estmod)) +
    geom_line() +
    geom_point() +
    labs(x = get_om_parlabel(index_df$study[1]),
         y = "RMSE")
}

evaluate_calibration <- function(index_df) {
  index_df %>%
    mutate(pnorm = pnorm(index_true, index_unb, unb_sd)) %>%
    ggplot(aes(x = pnorm, y = stat(density), fill = estmod)) +
    geom_histogram(breaks = seq(0.0, 1.0, 0.1)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_grid(opmod ~ estmod) +
    labs(x = "Quantile") +
    guides(fill = FALSE) +
    theme_minimal() +
    theme(axis.line.x = element_line(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

plot_index_devs <- function(index_df) {
  index_df %>%
    mutate(dev = index_est - index_true) %>%
    ggplot(aes(x = year, y = dev, color = estmod, group = repl)) +
    geom_line(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid(estmod ~ opmod) +
    guides(color = FALSE)
}

plot_bias <- function(index_df) {
  index_df %>%
    ggplot(aes(x = index_true, y = index_est, color = estmod, group = repl)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_grid(estmod ~ opmod) +
    coord_fixed() +
    guides(color = FALSE)
}

eval_dir <- "evaluation"
pdonly <- FALSE

for (study in studies) {
  index_df <- read_all_indices(results[[study]]$indexcsv) %>%
    filter(complete.cases(.)) %>%
    mutate(study = study,
           parval = map2_dbl(study, opmod, get_om_parval))

  convg_df <- map_dfr(results[[study]]$rdata,
                      ~ readRDS(.)$spec) %>%
    select(estmod, opmod, repl, Rdata) %>%
    mutate(repl = factor(repl, levels = 1:50),
            convcode = map_dbl(results[[study]]$rdata,
                              ~ readRDS(.)$fit$convergence),
            outer_mgc = map_dbl(results[[study]]$rdata,
                                ~ max(readRDS(.)$fit$grad)),
            pdhess = map_lgl(results[[study]]$rdata,
                            function(fn) {
              pdhess <- readRDS(fn)$sdr$pdHess
              if(is.null(pdhess)) pdhess <- FALSE
              return(pdhess)
            })) %>%
    select(estmod, opmod, repl, pdhess) %>%
    mutate(opmod = factor(opmod, levels = opmods),
           estmod = factor(estmod, levels = estmods))

  index_df <- left_join(index_df, convg_df,
                        by = c("estmod", "opmod", "repl"))

  pdhess_df <- index_df %>%
    group_by(estmod, opmod, repl) %>%
    summarize(pdhess = all(pdhess)) %>%
    summarize(pdhess = sum(pdhess)) %>%
    pivot_wider(names_from = opmod, values_from = pdhess)

  if (pdonly) {
      index_df <- filter(index_df, pdhess)
  }

  if (!file.exists(eval_dir))
    dir.create(eval_dir)

  if (!file.exists(file.path(eval_dir, study)))
    dir.create(file.path(eval_dir, study))

  bias_df <- evaluate_bias(index_df)
  bias_wide <- bias_df %>%
    pivot_wider(names_from = opmod,
                values_from = delta)
  rmse_df <- evaluate_rmse(index_df)
  rmse_wide <- rmse_df %>%
    pivot_wider(names_from = opmod,
                values_from = rmse)

  bias_plot <- plot_bias(index_df)
  bias2_plot <- plot_bias2(index_df)
  calibration_plot <- evaluate_calibration(index_df)
  index_devs <- plot_index_devs(index_df)
  rmse_plot <- plot_rmse2(index_df)

  write_csv(pdhess_df, file.path(eval_dir, study, "pdhess_wide.csv"))
  write_csv(bias_df, file.path(eval_dir, study, "bias.csv"))
  write_csv(bias_wide, file.path(eval_dir, study, "bias_wide.csv"))
  write_csv(rmse_df, file.path(eval_dir, study, "rmse.csv"))
  write_csv(rmse_wide, file.path(eval_dir, study, "rmse_wide.csv"))
  ggsave(file.path(eval_dir, study, "bias_plot.svg"), bias_plot)
  ggsave(file.path(eval_dir, study, "bias2_plot.svg"), bias2_plot)
  ggsave(file.path(eval_dir, study, "calibration.svg"), calibration_plot)
  ggsave(file.path(eval_dir, study, "index_devs.svg"), index_devs)
  ggsave(file.path(eval_dir, study, "rmse_plot.svg"), rmse_plot)
  ggsave(file.path(eval_dir, study, "bias_plot.png"), width = 7, height = 7, bias_plot)
  ggsave(file.path(eval_dir, study, "bias2_plot.png"), bias2_plot)
  ggsave(file.path(eval_dir, study, "calibration.png"), width = 7, height = 7, calibration_plot)
  ggsave(file.path(eval_dir, study, "index_devs.png"), width = 7, height = 7, index_devs)
  ggsave(file.path(eval_dir, study, "rmse_plot.png"), rmse_plot)
}
