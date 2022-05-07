library(spatq)
library(tidyverse)
library(patchwork)
library(hdf5r)

source("40_postproc_funs.R")

root_dir <- "."
eval_dir <- "om_results"
if (!dir.exists(eval_dir)) dir.create(eval_dir)
## Which simulation study are we fitting?
studies <- c("qdevscaling",
             "habq",
             "prefintensity",
             "densdepq",
             "bycatch")
names(studies) <- studies # So `lapply` outputs named list
## What range of replicates are going to be fit?
repls <- factor(1:100)
## How many years to fit?
max_T <- 15
## Names of the operating models
opmods <- factor(1:6)

spec_df <- cross_df(list(study = studies,
                         repl = repls,
                         opmod = opmods)) %>%
  rowwise() %>%
  mutate(paths = list(sim_file_paths(study, repl, opmod, root_dir)),
         has_pop_feather = file.exists(paths$pop_feather),
         has_pop_h5 = file.exists(paths$pop_h5))

pop_df <- spec_df %>%
  filter(has_pop_feather) %>%
  rowwise() %>%
  mutate(pop = list(read_popstate(study, repl, opmod, root_dir = root_dir, filetype = "feather")),
         studyf = factor(study, levels = studies)) %>%
  select(-has_pop_feather, -has_pop_h5) %>%
  unnest(pop) %>%
  filter(year <= 15)

traj_plot <- function(data) {
  data %>%
  mutate(parval = map2_dbl(study, opmod, get_om_parval),
         parval_lab = signif(parval, 2)) %>%
  ggplot(aes(x = year, y = pop, group = repl)) +
    geom_line(alpha = 0.2) +
    facet_grid(study ~ parval_lab) +
    ylim(0, 100) +
    scale_x_continuous(name = "Year", c(5, 10, 15)) +
    theme_bw()
}

pop_df_nest <- pop_df %>%
  group_by(studyf) %>%
  nest()
traj_subplots <- map(pop_df_nest$data, traj_plot)
pop_traj <- reduce(traj_subplots, `+`) +
  plot_layout(ncol = 1)
ggsave(file.path(eval_dir, "pop_trajectories.svg"),
       pop_traj,
       width = 8, height = 6)

study_depl_df <- pop_df %>%
  filter(year == 15) %>%
  group_by(study) %>%
  summarize(max = max(pop),
            min = min(pop),
            mean = mean(pop),
            median = median(pop),
            sd = sd(pop))
write_csv(study_depl_df,
          file.path(eval_dir, "study_depl.csv"))

opmod_depl_df <- pop_df %>%
  filter(year == 15) %>%
  group_by(study, opmod) %>%
  summarize(max = max(pop),
            min = min(pop),
            mean = mean(pop),
            median = median(pop),
            sd = sd(pop))
write_csv(opmod_depl_df,
          file.path(eval_dir, "opmod_depl.csv"))

depl_fig <- pop_df %>%
  filter(year == 15) %>%
  group_by(study, opmod) %>%
  ggplot(aes(x = opmod, y = pop)) +
  geom_boxplot() +
  xlab("Operating model") +
  ylab("Year 15 population") +
  facet_wrap(~ studyf, ncol = 1) +
  theme_bw()

ggsave(file.path(eval_dir, "depl_fig.svg"),
       width = 8, height = 6)

depl_boxplot <- function(data, ylims) {
  data %>%
    filter(year == 15) %>%
    ggplot(aes(x = factor(parval_lab), y = pop)) +
    geom_boxplot() +
    xlab(get_om_parlabel(data$study[1])) +
    ylim(ylims) +
    ylab("Final abundance") +
    theme_bw()
}

depl_df_nest <- pop_df %>%
  mutate(parval = map2_dbl(study, opmod, get_om_parval),
         parval_lab = signif(parval, 2)) %>%
  group_by(studyf) %>%
  nest()
ylims <- c(min(pop_df[pop_df$year == 15, "pop"]),
           max(pop_df[pop_df$year == 15, "pop"]))

depl_subplots <- map(depl_df_nest$data, depl_boxplot, ylims = ylims)
depl_plot <- reduce(depl_subplots, `+`) +
  plot_layout(ncol = 1)
ggsave(file.path(eval_dir, "depl_fig.svg"),
       width = 8, height = 8)


abundq_df1 <- spec_df %>%
  filter(study %in% c("qdevscaling", "prefintensity", "densdepq")) %>%
  mutate(pop = list(read_popstate(study, repl, opmod, root_dir, "h5")))
