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
library(patchwork)

source("40_postproc_funs.R")

root_dir <- "."
eval_dir <- "figs"
## Which simulation study are we fitting?
studies <- c("qdevscaling",
             "prefintensity",
             "densdepq",
             "habq",
             "bycatch")
names(studies) <- studies # So `lapply` outputs named list
## What range of replicates are going to be fit?
repls <- factor(1:100)
## How many years to fit?
max_T <- 15
## Names of the operating models
opmods <- factor(1:6)
## Names of the estimation models; can't use a factor because it's not
## recognized by om_* functions yet.
estmods <- c(## "model",              # Non-spatial survey only
             ## "survey",             # Survey-only
             "survey_spt",         # Spatiotemporal survey
             "sptemp_ab",         # All data, spatial abundance
             "sptemp_q")          # All data, spatial abundance + catchability

## Evaluate metrics and save ----------------------------------------------------
studypost <- lapply(studies, postproc)
saveRDS(studypost, file.path(eval_dir, "studypost.Rdata"))
sapply(studypost, save_tables, eval_dir = eval_dir)
sapply(studypost, save_plots, eval_dir = eval_dir)

## Collect metrics into a manuscript-ready figure -------------------------------
## Including each of the `collect_` functions in this script because they're
## only used for MS figures and they depend on which studies are being used.
collect_bias2_plots <- function(studypost) {
  ## Get y-axis limits so everything is on the same scale
  bias2_ylims <- get_bias2_ylims(studypost)

  ## make blank plot with y-label only
  ylab_plot <- make_ylab_plot("δ bias metric")

  ## Make the plot
  ylab_plot +
    (refine_bias2_plot(studypost, "qdevscaling", bias2_ylims) +
     refine_bias2_plot(studypost, "habq", bias2_ylims)) /
    (refine_bias2_plot(studypost, "prefintensity", bias2_ylims) +
     refine_bias2_plot(studypost, "densdepq", bias2_ylims)) /
    (refine_bias2_plot(studypost, "bycatch", bias2_ylims) +
     plot_spacer()) +
    plot_layout(widths = c(1, 49),
                guides = "collect") &
    theme(legend.position = "bottom")
}

bias_all <- collect_bias2_plots(studypost)
ggsave(file.path(eval_dir, "bias_all.svg"),
       bias_all,
       width = 8, height = 8)

collect_rmse_plots <- function(studypost) {
  rmse_ylims <- get_rmse_ylims(studypost)

  ylab_plot <- make_ylab_plot("RMSE")

  ylab_plot +
    (refine_rmse_plot(studypost, "qdevscaling", rmse_ylims) +
     refine_rmse_plot(studypost, "habq", rmse_ylims)) /
    (refine_rmse_plot(studypost, "prefintensity", rmse_ylims) +
     refine_rmse_plot(studypost, "densdepq", rmse_ylims)) /
    (refine_rmse_plot(studypost, "bycatch", rmse_ylims) +
     plot_spacer()) +
    plot_layout(widths = c(1, 49),
                guides = "collect") &
    theme(legend.position = "bottom")
}

rmse_all <- collect_rmse_plots(studypost)
ggsave(file.path(eval_dir, "rmse_all.svg"),
       rmse_all,
       width = 8, height = 8)

collect_calibration_plots <- function(studypost) {
  studypost$qdevscaling$calibration_plot +
    ggtitle(get_study_title(studypost$qdevscaling$study)) +
    studypost$habq$calibration_plot +
    ggtitle(get_study_title(studypost$habq$study)) +
    studypost$prefintensity$calibration_plot +
    ggtitle(get_study_title(studypost$prefintensity$study)) +
    studypost$densdepq$calibration_plot +
    ggtitle(get_study_title(studypost$densdepq$study)) +
    studypost$bycatch$calibration_plot +
    ggtitle(get_study_title(studypost$bycatch$study)) +
    plot_layout(ncol = 2, nrow = 3, byrow = TRUE,
                guides = "collect")
}

calibration_all <- collect_calibration_plots(studypost)
ggsave(file.path(eval_dir, "calibration_all.svg"),
       calibration_all,
       width = 8, height = 8)

collect_mapcor_plots <- function(studypost) {
  mapcor_ylims <- get_mapcor_ylims(studypost)

  ylab_plot <- make_ylab_plot("Spearman rank correlation")

  ylab_plot +
    (refine_mapcor_plot(studypost, "qdevscaling", mapcor_ylims) +
     refine_mapcor_plot(studypost, "habq", mapcor_ylims)) /
    (refine_mapcor_plot(studypost, "prefintensity", mapcor_ylims) +
     refine_mapcor_plot(studypost, "densdepq", mapcor_ylims)) /
    (refine_mapcor_plot(studypost, "bycatch", mapcor_ylims) +
     plot_spacer()) +
    plot_layout(widths = c(1, 49),
                guides = "collect") &
    theme(legend.position = "bottom")
}

mapcor_all <- collect_mapcor_plots(studypost)
ggsave(file.path(eval_dir, "mapcor_all.svg"),
       mapcor_all,
       width = 8, height = 8)
