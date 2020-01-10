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

mgc_getter <- attr_getter("mgc")
get_mgc <- function(fitlist) {
  pluck(fitlist, "fit", mgc_getter)
}

read_index_csv <- function(filename) {
  scenarios <- c("combo", "pref", "spatq")
  data_specs <- c("all", "indep")
  read_csv(filename,
           col_types = cols(repl = col_factor(levels = 1:100),
                            scenario = col_factor(levels = scenarios),
                            data_spec = col_factor(levels = data_specs),
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
    group_by(scenario, data_spec) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(log(raw_est) ~ repl + log(raw_true), data = .x)),
           coef = map(mod, coef),
           delta = map_dbl(coef, pluck, "log(raw_true)")) %>%
    select(scenario, data_spec, delta)
}

evaluate_rmse <- function(index_df) {
  index_df %>%
    mutate(sq_err = (index_est - index_true)^2) %>%
    group_by(scenario, data_spec) %>%
    summarize(rmse = sqrt(mean(sq_err))) %>%
    select(scenario, data_spec, rmse)
}

evaluate_calibration <- function(index_df) {
  index_df %>%
    mutate(pnorm = pnorm(index_true, index_est, sd)) %>%
    ggplot(aes(x = pnorm, y = stat(density), fill = data_spec)) +
    geom_histogram(breaks = seq(0.0, 1.0, 0.1)) +
    geom_hline(yintercept = 1) +
    facet_grid(scenario ~ data_spec)
}

plot_index_devs <- function(index_df) {
  index_df %>%
    mutate(dev = index_est - index_true) %>%
    ggplot(aes(x = year, y = dev, color = data_spec, group = repl)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid(scenario ~ data_spec)
}

eval_main <- function(results_dir = "results", eval_dir = "results/evaluation") {
  index_df <- read_all_indices(results_dir)

  if (!file.exists(eval_dir))
    dir.create(eval_dir)

  bias_df <- evaluate_bias(index_df)
  rmse_df <- evaluate_rmse(index_df)

  calibration_plot <- evaluate_calibration(index_df)
  index_devs <- plot_index_devs(index_df)

  write_csv(bias_df, file.path(eval_dir, "bias.csv"))
  write_csv(rmse_df, file.path(eval_dir, "rmse.csv"))
  ggsave(file.path(eval_dir, "calibration.pdf"), calibration_plot)
  ggsave(file.path(eval_dir, "index_devs.pdf"), index_devs)
}

eval_main()
