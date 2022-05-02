## if script is run in the REPL interactively, use the local `spatq` package and
## manually set a range of replicates to fit, because they were probably not
## passed as command line arguments when R was started. Otherwise (typically on
## Hyak via SLURM), use the installed version of `spatq` and read in the
## replicate numbers from the command line arguments.
if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
} else {
  library(spatq)
}
library(tidyverse)

root_dir <- "."
## Which simulation study are we fitting?
studies <- c("qdevscaling",
             ## "sharedq",
             "prefintensity",
             "densdepq",
             "habq",
             "bycatch")
## What range of replicates are going to be fit?
repls <- factor(1:100)
## How many years to fit?
max_T <- 15
## Names of the operating models
opmods <- factor(1:6)
## Names of the estimation models; can't use a factor because it's not
## recognized by om_* functions yet.
estmods <- c("model",              # Non-spatial survey only
             "survey_spt",         # Spatiotemporal survey
             "sptemp_ab",          # All data, spatial abundance
             "sptemp_q")           # All data, spatial abundance + catchability

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
                      prefintensity = c(0, 1, 2, 4, 8, 16),
                      densdepq = seq(0, 1.25, 0.25),
                      habq = 2.0 ^ (-2:3),
                      bycatch = seq(0, 1, 0.2))
  studyvals[opmod]
}

get_om_parlabel <- function(study) {
  switch(study,
         qdevscaling = "log catchability deviation SD",
         sharedq = "Prop shared catchabilty dev",
         prefintensity = "Preference power",
         densdepq = "Density dependent multiplier",
         habq = "Rocky habitat preference",
         bycatch = "Comm catchability reduction")
}

evaluate_bias <- function(index_df) {
  index_df %>%
    group_by(opmod, estmod) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(log(raw_unb) ~ repl + log(raw_true),
                                data = .x)),
           coef = map(mod, coef),
           delta = map_dbl(coef, pluck, "log(raw_true)"),
           sigma = map_dbl(mod, sigma),
           ci = map(mod, confint, parm = "log(raw_true)"),
           ci_lower = map_dbl(ci, pluck, 1),
           ci_upper = map_dbl(ci, pluck, 2)) %>%
    select(opmod, estmod, delta, sigma, ci_lower, ci_upper) %>%
    mutate(parval = map2_dbl(study, opmod, get_om_parval),
           estmod = factor(estmod, levels = estmods))
}

plot_bias2 <- function(index_df) {
  study <- index_df$study[1]
  opmods <- unique(index_df$opmod)
  estmods <- unique(index_df$estmod)

  bias_df <- evaluate_bias(index_df) %>%
    mutate(parval = factor(parval))

  x_pos <- position_dodge(width = 0.5)

  ggplot(bias_df, aes(x = parval, y = delta,
                      ymin = ci_lower, ymax = ci_upper,
                      color = estmod, group = estmod)) +
    geom_line(position = x_pos) +
    geom_point(position = x_pos) +
    geom_errorbar(position = x_pos, width = 1 / length(estmods)) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    scale_x_discrete(labels = signif(get_om_parval(study, opmods), 2)) +
    labs(x = get_om_parlabel(study),
         y = "δ bias metric")
}

evaluate_rmse <- function(index_df) {
  index_df %>%
    mutate(sq_err = (index_unb - index_true)^2) %>%
    group_by(opmod, estmod) %>%
    summarize(rmse = sqrt(mean(sq_err)), .groups = "drop") %>%
    select(opmod, estmod, rmse) %>%
    mutate(parval = map2_dbl(study, opmod, get_om_parval),
           estmod = factor(estmod, levels = estmods))
}

plot_rmse2 <- function(index_df) {
  study <- index_df$study[1]
  opmods <- unique(index_df$opmod)
  estmods <- unique(index_df$estmod)
  rmse_df <- evaluate_rmse(index_df) %>%
    mutate(parval = factor(parval))
  ggplot(rmse_df, aes(x = parval, y = rmse, color = estmod, group = estmod)) +
    geom_line() +
    geom_point() +
    scale_x_discrete(labels = signif(get_om_parval(study, opmods), 2)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.1), c(0, 0))) +
    labs(x = get_om_parlabel(study),
         y = "RMSE")
}

evaluate_calibration <- function(index_df) {
  index_df %>%
    mutate(pnorm = pnorm(index_true, index_unb, unb_sd),
           estmod = factor(estmod, levels = estmods)) %>%
    ggplot(aes(x = pnorm, y = stat(density), fill = estmod)) +
    geom_histogram(breaks = seq(0.0, 1.0, 0.1)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_grid(opmod ~ estmod) +
    labs(x = "Quantile") +
    guides(fill = "none") +
    theme_minimal() +
    theme(axis.line.x = element_line(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

plot_index_devs <- function(index_df) {
  index_df %>%
    mutate(dev = index_est - index_true,
           estmod = factor(estmod, levels = estmods)) %>%
    ggplot(aes(x = year, y = dev, color = estmod, group = repl)) +
    geom_line(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid(estmod ~ opmod) +
    guides(color = "none")
}

plot_bias <- function(index_df) {
  index_df %>%
    mutate(estmod = factor(estmod, levels = estmods)) %>%
    ggplot(aes(x = index_true, y = index_est, color = estmod, group = repl)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_grid(estmod ~ opmod) +
    coord_fixed() +
    guides(color = "none")
}

evaluate_mapcor <- function(res_df) {
  res_df %>%
    select(study, repl, opmod, estmod, spec, res) %>%
    group_by(study, repl, opmod) %>%
    ## Nest here so that the true populations are read once each instead of once
    ## for each estimation model
    nest() %>%
    rowwise() %>%
    mutate(true = list(read_popstate(study, repl, opmod, root_dir = root_dir, filetype = "h5"))) %>%
    unnest(cols = data) %>%
    mutate(est = map(res, read_estpop),
           cor = map2_dbl(true, est, map_correlation, method = "spearman")) %>%
    select(study, repl, opmod, estmod, cor)
}

plot_mapcor <- function(cor_df) {
    ggplot(cor_df, aes(x = estmod, y = cor, color = estmod)) +
    geom_boxplot() +
      facet_grid(~ opmod) +
      guides(color = "none")
}

eval_dir <- "evaluation"
pdonly <- TRUE

for (study in studies) {
  spec_df <- cross_df(list(study = study,
                           repl = repls,
                           opmod = opmods,
                           estmod = estmods)) %>%
    rowwise() %>%
    mutate(spec = list(spatq_simstudyspec(study, repl = repl, opmod = opmod, estmod = estmod)),
           paths = list(res_file_paths(study, repl, opmod, estmod, root_dir)),
           has_rdata = file.exists(paths$rdata))

  index_df <- map_df(spec_df$spec, read_index, estmods = estmods)

  res_df <- spec_df %>%
    filter(has_rdata) %>%
    mutate(res = list(read_rdata(spec)),
           pdhess = posdefhess(res))

  pdhess_df <- res_df %>%
    group_by(opmod, estmod) %>%
    summarize(pdhess = sum(pdhess),
              .groups = "drop") %>%
    pivot_wider(names_from = opmod, values_from = pdhess)

  index_df <- res_df %>%
    select(study, repl, opmod, estmod, pdhess) %>%
    right_join(index_df, by = c("study", "repl", "opmod", "estmod"))

  if (pdonly)
      index_df <- filter(index_df, pdhess)

  if (!file.exists(eval_dir))
    dir.create(eval_dir)

  if (!file.exists(file.path(eval_dir, study)))
    dir.create(file.path(eval_dir, study))

  bias_df <- evaluate_bias(index_df)
  bias_wide <- bias_df %>%
    ungroup() %>%
    select(-opmod, -ci_lower, -ci_upper, -sigma) %>%
    pivot_wider(names_from = parval,
                values_from = delta)

  rmse_df <- evaluate_rmse(index_df)
  rmse_wide <- rmse_df %>%
    select(-opmod) %>%
    pivot_wider(names_from = parval,
                values_from = rmse)

  bias_plot <- plot_bias(index_df)
  bias2_plot <- plot_bias2(index_df)
  calibration_plot <- evaluate_calibration(index_df)
  index_devs <- plot_index_devs(index_df)
  rmse_plot <- plot_rmse2(index_df)

  cor_df <- evaluate_mapcor(res_df)
  cor_plot <- plot_mapcor(cor_df)

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
  ggsave(file.path(eval_dir, study, "mapcor_plot.svg"), cor_plot)
  ggsave(file.path(eval_dir, study, "bias_plot.png"), width = 7, height = 7, bias_plot)
  ggsave(file.path(eval_dir, study, "bias2_plot.png"), width = 10, height = 7, bias2_plot)
  ggsave(file.path(eval_dir, study, "calibration.png"), width = 7, height = 7, calibration_plot)
  ggsave(file.path(eval_dir, study, "index_devs.png"), width = 7, height = 7, index_devs)
  ggsave(file.path(eval_dir, study, "rmse_plot.png"), width = 7, height = 7, rmse_plot)
  ggsave(file.path(eval_dir, study, "mapcor_plot.png"), width = 10, height = 7, cor_plot)
}

#rmarkdown::render("41_postproc_report.Rmd")
