if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
} else {
  library(spatq)
}
library(tidyverse)

estmods <- c("design",
             "model",
             "survey")#,
             ## "survey_spt")

ft <- function(estmod) {
  switch(estmod,
         design = "feather",
         model = "feather",
         survey = "feather",
         survey_spt = "csv")
}

em_df <- cross_df(list(
  study = "prefintensity",
  repl = 1:100,
  opmod = 1:6,
  estmod = estmods)) %>%
  rowwise() %>%
  mutate(spec = list(spatq_simstudyspec(study, repl, opmod, estmod)),
         filetype = ft(estmod))

index_df <- map2_df(em_df$spec, em_df$filetype,
                    read_index, estmods = estmods)

bias_df <- index_df %>%
  group_by(opmod, estmod) %>%
  nest() %>%
  mutate(bias_res = map(data, bias_metric),
         delta = map_dbl(bias_res, pluck, "delta")) %>%
  select(opmod, estmod, delta)
bias_wide <- bias_df %>%
  pivot_wider(names_from = opmod, values_from = delta)

ggplot(bias_df, aes(x = as.numeric(opmod), y = delta, color = estmod)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed")

ggplot(index_df, aes(x = index_true, y = index_unb)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point() +
  facet_grid(opmod ~ estmod) +
  geom_smooth(method = lm)

ggplot(index_df, aes(x = log(raw_true), y = log(raw_est))) +
  geom_abline(slope = 1, intercept = seq(-20, 0, 2), alpha = 0.4, linetype = "dashed") +
  geom_point() +
  geom_smooth(method = lm) +
  facet_grid(opmod ~ estmod, scales = "free")
