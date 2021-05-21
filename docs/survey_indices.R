devtools::load_all("~/dev/spatq", helpers = FALSE)
library(tidyverse)

spec_ls <- cross(list(study = c("qdevscaling",
                                "sharedq",
                                "prefintensity"),
                      repl = 1:25,
                      opmod = 1:6,
                      estmod = "survey",
                      root_dir = "archive")) %>%
  map(spatq_simstudyspec)

index_df <- map_df(spec_ls, read_index)

bias_df <- index_df %>%
  select(study, repl, opmod, estmod, year, raw_est, raw_unb, raw_true) %>%
  group_by(study, opmod, estmod) %>%
  summarize(raw_lm = list(lm(log(raw_est) ~ log(raw_true) + repl)),
            unb_lm = list(lm(log(raw_unb) ~ log(raw_true) + repl)),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(delta_raw = coef(raw_lm)["log(raw_true)"],
         delta_unb = coef(unb_lm)["log(raw_true)"])

bias_df %>%
  select(study, opmod, estmod, delta_raw, delta_unb) %>%
  ggplot(aes(x = delta_raw, y = delta_unb)) +
  geom_point() +
  facet_wrap(~ study) +
  geom_abline(slope = 1, intercept = 0)

bias2_df <- index_df %>%
  select(study, repl, opmod, estmod, year, raw_est, raw_unb, raw_true) %>%
  group_by(study, opmod, estmod) %>%
  summarize(raw_lm = list(lm(log(raw_true) ~ log(raw_est) + repl)),
            unb_lm = list(lm(log(raw_true) ~ log(raw_unb) + repl)),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(delta_raw = coef(raw_lm)["log(raw_est)"],
         delta_unb = coef(unb_lm)["log(raw_unb)"])

