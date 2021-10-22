library(spatq)
library(tidyverse)

root_dir <- "."
eval_dir <- "om_results"
## Which simulation study are we fitting?
studies <- c("qdevscaling",
             "habq",
             "prefintensity",
             "densdepq")
names(studies) <- studies # So `lapply` outputs named list
## What range of replicates are going to be fit?
repls <- factor(1:100)
## How many years to fit?
max_T <- 15
## Names of the operating models
opmods <- factor(1:6)

get_om_parval <- function(study, opmod = 1:6) {
  qds <- 10 ^ (seq(-3, -0.5, 0.5))
  qdcv <- sqrt(exp(qds ^ 2) - 1)
  studyvals <- switch(study,
                      qdevscaling = qdcv,
                      habq = 1.75 ^ (-2:3),
                      prefintensity = c(0, 1, 2, 4, 8, 16),
                      densdepq = seq(0, 1.25, 0.25))
  studyvals[opmod]
}

get_om_parlabel <- function(study) {
  switch(study,
         qdevscaling = "Spatial catchability CV",
         habq = "Relative habitat preference",
         prefintensity = "Preference power",
         densdepq = "Density dependent multiplier")
}

get_study_title <- function(study) {
  switch(study,
         qdevscaling = "a. Spatial catchability variability",
         habq = "b. Habitat-dependent catchability",
         prefintensity = "c. Preference intensity",
         densdepq = "d. Density-dependent catchability")
}


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
  mutate(pop = list(read_popstate(study, repl, opmod, root_dir = root_dir, filetype = "feather"))) %>%
  select(-has_pop_feather, -has_pop_h5) %>%
  unnest(pop) %>%
  filter(year <= 15)

pop_traj <- ggplot(pop_df, aes(x = year, y = pop, group = repl)) +
  geom_line() +
  facet_grid(study ~ opmod) +
  ylim(0, 100)
ggsave(file.path(eval_dir, "pop_trajectories.svg"),
       pop_traj,
       width = 8, height = 6)
