library(tidyverse)

list_Rdata <- function(root_dir = "results") {
  file.path(root_dir,
            list.files(root_dir, pattern = "\\.Rdata$"))
}

list_csvs <- function(root_dir = "results") {
  file.path(root_dir,
            list.files(root_dir, pattern = "\\.csv$"))
}

get_convcode <- function(fitlist) {
  pluck(fitlist, "fit", "convergence")
}

get_mgc <- function(fitlist) {
  pluck(fitlist, "fit", attr_getter("mgc"))
}

read_index_csv <- function(filename) {
  opmods <- c("combo", "pref", "spat")
  estmods <- c("surv", "spatab", "spatq")
  read_csv(filename,
           col_types = cols(repl = col_factor(levels = 1:50),
                            opmod = col_factor(levels = opmods),
                            estmod = col_factor(levels = estmods),
                            year = col_integer(),
                            raw_est = col_double(),
                            index_est = col_double(),
                            sd_raw = col_double(),
                            sd = col_double(),
                            raw_true = col_double(),
                            index_true = col_double()))
}

read_all_indices <- function(root_dir = "results") {
  filenames <- list_csvs(root_dir)
  map_df(filenames, read_index_csv)
}

evaluate_bias <- function(index_df) {
  index_df %>%
    group_by(opmod, estmod) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(log(raw_est) ~ repl + log(raw_true), data = .x)),
           coef = map(mod, coef),
           delta = map_dbl(coef, pluck, "log(raw_true)")) %>%
    select(opmod, estmod, delta)
}

evaluate_rmse <- function(index_df) {
  index_df %>%
    mutate(sq_err = (index_est - index_true)^2) %>%
    group_by(opmod, estmod) %>%
    summarize(rmse = sqrt(mean(sq_err))) %>%
    select(opmod, estmod, rmse)
}

evaluate_calibration <- function(index_df) {
  index_df %>%
    mutate(pnorm = pnorm(index_true, index_est, sd)) %>%
    ggplot(aes(x = pnorm, y = stat(density), fill = estmod)) +
    geom_histogram(breaks = seq(0.0, 1.0, 0.1)) +
    geom_hline(yintercept = 1) +
    facet_grid(opmod ~ estmod) +
    guides(fill = FALSE)
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

eval_main <- function(results_dir = "results",
                      eval_dir = "results/evaluation",
                      pdonly = TRUE) {
  index_df <- read_all_indices(results_dir) %>% filter(complete.cases(.))

  convg_df <- map_dfr(list_Rdata(), ~ readRDS(.)$spec) %>%
    select(estmod, opmod, repl, Rdata) %>%
    mutate(repl = factor(repl, levels = 1:50),
           convcode = map_dbl(Rdata, ~ readRDS(.)$fit$convergence),
           outer_mgc = map_dbl(Rdata, ~ max(readRDS(.)$fit$grad)),
           pdhess = map_lgl(Rdata, function(fn) {
             pdhess <- readRDS(fn)$sdr$pdHess
             if(is.null(pdhess)) pdhess <- FALSE
             return(pdhess)
           })) %>%
    select(estmod, opmod, repl, pdhess)

  if (pdonly) {
    index_df <- left_join(index_df, convg_df,
                          by = c("estmod", "opmod", "repl")) %>%
      filter(pdhess)
  }

  if (!file.exists(eval_dir))
    dir.create(eval_dir)

  bias_df <- evaluate_bias(index_df)
  rmse_df <- evaluate_rmse(index_df)

  bias_plot <- plot_bias(index_df)
  calibration_plot <- evaluate_calibration(index_df)
  index_devs <- plot_index_devs(index_df)

  write_csv(bias_df, file.path(eval_dir, "bias.csv"))
  write_csv(rmse_df, file.path(eval_dir, "rmse.csv"))
  ggsave(file.path(eval_dir, "bias_plot.pdf"), bias_plot)
  ggsave(file.path(eval_dir, "calibration.pdf"), calibration_plot)
  ggsave(file.path(eval_dir, "index_devs.pdf"), index_devs)
  ggsave(file.path(eval_dir, "bias_plot.png"), width = 7, height = 7, bias_plot)
  ggsave(file.path(eval_dir, "calibration.png"), width = 7, height = 7, calibration_plot)
  ggsave(file.path(eval_dir, "index_devs.png"), width = 7, height = 7, index_devs)
}

eval_main()
