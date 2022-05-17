## Fit convenience functions ----------------------------------------------------
get_convcode <- function(fitlist) {
  pluck(fitlist, "fit", "convergence")
}

get_mgc <- function(fitlist) {
  pluck(fitlist, "fit", attr_getter("mgc"))
}

## Plotting convenience functions -----------------------------------------------
get_om_parval <- function(study, opmod = 1:6) {
  qds <- 10 ^ (seq(-3, -0.5, 0.5))
  qdcv <- sqrt(exp(qds ^ 2) - 1)
  studyvals <- switch(study,
                      qdevscaling = qdcv,
                      sharedq = seq(0, 1, 0.2),
                      prefintensity = c(0, 1, 2, 4, 8, 16),
                      densdepq = seq(0, 1.25, 0.25),
                      habq = 2 ^ (-2:3),
                      bycatch = seq(0, 1, 0.2))
  studyvals[opmod]
}

get_om_parlabel <- function(study) {
  switch(study,
         qdevscaling = "Spatial catchability CV",
         sharedq = "Proportion shared catchability devs",
         prefintensity = "Preference power",
         densdepq = "Density dependent multiplier",
         habq = "Relative habitat preference",
         bycatch = "Catchability reduction in bycatch areas")
}

get_study_title <- function(study) {
  switch(study,
         qdevscaling = "a. Spatial catchability variability",
         habq = "b. Habitat-dependent catchability",
         sharedq = "x. Shared catchability deviations",
         prefintensity = "c. Preference intensity",
         densdepq = "d. Density-dependent catchability",
         bycatch = "e. Bycatch avoidance")
}

get_base_om <- function(study, opmod) {
  switch(study,
         qdevscaling = 1,
         habq = NA,
         prefintensity = 1,
         densdepq = 1)
}

## Evaluate each metric ---------------------------------------------------------
evaluate_bias <- function(index_df) {
  index_df %>%
    group_by(study, opmod, estmod) %>%
    nest() %>%
    mutate(mod = map(data, ~ lm(log(raw_unb) ~ repl + log(raw_true),
                                data = .x)),
           coef = map(mod, coef),
           delta = map_dbl(coef, pluck, "log(raw_true)"),
           sigma = map_dbl(mod, sigma),
           ci = map(mod, confint, parm = "log(raw_true)"),
           ci_lower = map_dbl(ci, pluck, 1),
           ci_upper = map_dbl(ci, pluck, 2)) %>%
    select(study, opmod, estmod, delta, sigma, ci_lower, ci_upper) %>%
    mutate(parval = map2_dbl(study, opmod, get_om_parval),
           estmod = factor(estmod, levels = estmods))
}

evaluate_rmse <- function(index_df) {
  index_df %>%
    mutate(sq_err = (index_unb - index_true)^2) %>%
    group_by(study, opmod, estmod) %>%
    summarize(rmse = sqrt(mean(sq_err)), .groups = "drop") %>%
    select(study, opmod, estmod, rmse) %>%
    mutate(parval = map2_dbl(study, opmod, get_om_parval),
           estmod = factor(estmod, levels = estmods))
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
    mutate(
      est = map(res, read_estpop),
      cor = map2_dbl(true, est, map_correlation, method = "spearman", byyear = TRUE),
      estmod = factor(estmod, levels = estmods),
      parval = map2_dbl(study, opmod, get_om_parval)
    ) %>%
    select(study, repl, opmod, parval, estmod, cor)
}

## Plot each metric -------------------------------------------------------------
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
         y = "δ bias metric",
         color = "Estimation\nmodel")
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
         y = "RMSE",
         color = "Estimation\nmodel")
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

plot_mapcor <- function(cor_df) {
  study <- cor_df$study[1]
  ggplot(cor_df, aes(x = factor(signif(parval, 2)), y = cor, color = estmod)) +
    geom_boxplot() +
    labs(x = get_om_parlabel(study),
         y = "Spearman correlation",
         color = "Estimation\nmodel")
}

## Create plots and tables ------------------------------------------------------
postproc <- function(study, pdonly = TRUE) {
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
    select(-opmod, -sigma, -ci_lower, -ci_upper) %>%
    pivot_wider(names_from = parval,
                values_from = delta)

  rmse_df <- evaluate_rmse(index_df)
  rmse_wide <- rmse_df %>%
    select(-opmod) %>%
    pivot_wider(names_from = parval,
                values_from = rmse)

  mapcor_df <- evaluate_mapcor(res_df)

  bias_plot <- plot_bias(index_df)
  bias2_plot <- plot_bias2(index_df)
  calibration_plot <- evaluate_calibration(index_df)
  index_devs <- plot_index_devs(index_df)
  rmse_plot <- plot_rmse2(index_df)
  mapcor_plot <- plot_mapcor(mapcor_df)

  list(study = study,
       pdhess_df = pdhess_df,
       bias_df = bias_df,
       bias_wide = bias_wide,
       rmse_df = rmse_df,
       rmse_wide = rmse_wide,
       bias_plot = bias_plot,
       bias2_plot = bias2_plot,
       calibration_plot = calibration_plot,
       index_devs = index_devs,
       rmse_plot = rmse_plot,
       mapcor_plot = mapcor_plot)
}

save_tables <- function(studypost, eval_dir) {
  study <- studypost$study
  write_csv(studypost$pdhess_df, file.path(eval_dir, study, "pdhess_wide.csv"))
  write_csv(studypost$bias_df, file.path(eval_dir, study, "bias.csv"))
  write_csv(studypost$bias_wide, file.path(eval_dir, study, "bias_wide.csv"))
  write_csv(studypost$rmse_df, file.path(eval_dir, study, "rmse.csv"))
  write_csv(studypost$rmse_wide, file.path(eval_dir, study, "rmse_wide.csv"))
  invisible(TRUE)
}

save_plots <- function(studypost, eval_dir, width = 6, height = 4) {
  study <- studypost$study
  ggsave(file.path(eval_dir, study, "bias_plot.svg"),
         studypost$bias_plot,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "bias2_plot.svg"),
         studypost$bias2_plot,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "calibration.svg"),
         studypost$calibration_plot,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "index_devs.svg"),
         studypost$index_devs,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "rmse_plot.svg"),
         studypost$rmse_plot,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "bias_plot.png"),
         studypost$bias_plot,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "bias2_plot.png"),
         studypost$bias2_plot,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "calibration.png"),
         studypost$calibration_plot,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "index_devs.png"),
         studypost$index_devs,
         width = width, height = height)
  ggsave(file.path(eval_dir, study, "rmse_plot.png"),
         studypost$rmse_plot,
         width = width, height = height)
  invisible(TRUE)
}

## Get ylims for each plot type -------------------------------------------------
get_bias2_ylims <- function(studypost) {
  bias_data <- map(studypost, pluck, "bias2_plot", "data")
  bias_min <- map_dbl(bias_data, ~ min(.$ci_lower))
  bias_max <- map_dbl(bias_data, ~ max(.$ci_upper))
  c(min(bias_min), max(bias_max))
}

get_rmse_ylims <- function(studypost) {
  rmse_data <- map(studypost, pluck, "rmse_plot", "data")
  rmse_max <- map_dbl(rmse_data, ~ max(.$rmse))
  c(0, max(rmse_max))
}

get_mapcor_ylims <- function(studypost) {
  mapcor_data <- map(studypost, pluck, "mapcor_plot", "data")
  mapcor_min <- map_dbl(mapcor_data, ~ min(.$cor))
  c(min(mapcor_min), 1)
}

## Fix ylims and add title ------------------------------------------------------
refine_bias2_plot <- function(studypost, study, ylims) {
  studypost[[study]]$bias2_plot +
    ylim(ylims) +
    ylab(NULL) +
    ggtitle(get_study_title(study))
}

refine_rmse_plot <- function(studypost, study, ylims) {
  studypost[[study]]$rmse_plot +
    ylim(ylims) +
    ylab(NULL) +
    ggtitle(get_study_title(study))
}

refine_mapcor_plot <- function(studypost, study, ylims) {
  studypost[[study]]$mapcor_plot +
    ylim(ylims) +
    ylab(NULL) +
    ggtitle(get_study_title(study))
}

## Patchwork plot to hold y-axis label ------------------------------------------
make_ylab_plot <- function(label) {
  ggplot(data.frame(l = label, x = 1, y = 1)) +
      geom_text(aes(x, y, label = l), angle = 90) +
      theme_void() +
      coord_cartesian(clip = "off")
}
