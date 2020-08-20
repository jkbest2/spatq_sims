library(hdf5r)
library(tidyverse)

poparray_file <- function(repl, opmod, root_dir = ".") {
  repl_str <- str_pad(repl, width = 2, pad = 0)
  repl_dir <- paste0("repl_", repl_str)
  h5_file <- paste0("pop_", repl_str, ".h5")
  file.path(root_dir, repl_dir, h5_file)
}

read_poparray <- function(repl, opmod, root_dir = ".") {
  fn <- poparray_file(repl, opmod, root_dir)
  file <- H5File$new(fn, mode = "r")
  grp <- file[[opmod]]
  pop <- grp[["popstate"]]
  pop$read()
}

survey_stations <- function(st_s1 = seq(2, 98, 4), st_s2 = seq(5, 95, 10)) {
  stations <- matrix(0, 100, 100)
  stations[st_s1, st_s2] <- 1
  stations
}

## Survey stations: 2:4:98, 5:10:95
naive_index <- function(popmat,
                        prop_indep = 250 / 2750,
                        iprob_fun = survey_stations,
                        dprob_fun = identity) {
  iprob <- iprob_fun()
  iprob <- iprob / sum(iprob)

  dprob <- dprob_fun(popmat)
  dprob <- dprob / sum(dprob)

  I_indep <- sum(iprob * popmat)
  I_dep <- sum(dprob * popmat)

  prop_indep * I_indep + (1 - prop_indep) * I_dep
}

indep_index <- function(popmat, iprob_fun = survey_stations) {
  naive_index(popmat, prop_indep = 1,
              iprob_fun = iprob_fun,
              dprob_fun = identity)
}

naive_index_df <- function(poparray, prop_indep, max_T = 25) {
  index <- apply(poparray[, , seq_len(max_T)], 3,
                 naive_index, prop_indep = prop_indep)
  df <- tibble(year = seq_along(index),
               index_raw = index) %>%
    filter(year <= max_T) %>%
    mutate(index_exp = rescale_index(index_raw)) %>%
    select(year, index_raw, index_exp)
  df
}

max_T <- 15

index_df <- read_all_indices() %>%
  filter(year <= max_T)

naive_df <- cross_df(list(opmod = "combo",
                          prop_indep = c(1, 250 / 2750),
                          repl = factor(1:5, levels = 1:5))) %>%
  mutate(estmod = factor(ifelse(prop_indep < 1, "all", "indep")),
         poparray = map2(repl, opmod, read_poparray),
         index_df = map2(poparray, prop_indep,
                         ~ naive_index_df(.x, .y, max_T = max_T))) %>%
  select(opmod, estmod, repl, index_df) %>%
  unnest(index_df)

fullindex_df <- right_join(index_df, naive_df,
                           by = c("repl", "estmod", "year")) %>%
  mutate(exp_dev = index_exp - index_true,
         estexp_dev = index_est - index_exp) %>%
  select(repl, opmod = opmod.x, estmod, year, index_est, sd, index_true,
         index_exp, index_raw, raw_est, exp_dev, estexp_dev)

## Evaluate overall bias
bias_df <- fullindex_df %>%
  group_by(opmod, estmod) %>%
  nest() %>%
  mutate(mod_est = map(data,
                       ~ lm(log(raw_est) ~ log(index_true),
                            data = .x)),
         coef_est = map(mod_est, coef),
         delta_est = map_dbl(coef_est, pluck, "log(index_true)"),
         mod_exp = map(data,
                       ~ lm(log(index_raw) ~ 0 + log(index_true),
                            data = .x)),
         coef_exp = map(mod_exp, coef),
         delta_exp = map_dbl(coef_exp, pluck, "log(index_true)")) %>%
select(opmod, estmod, delta_est, delta_exp)

## Evaluate bias by year; funky behavior at intermediate years
biasyear_df <- fullindex_df %>%
  group_by(opmod, estmod) %>%
  nest() %>%
  mutate(mod_est = map(data,
                       ~ lm(log(index_est) ~ factor(year):log(index_true),
                            data = .x)),
         coef_est = map(mod_est, coef),
         mod_exp = map(data,
                       ~ lm(log(index_exp) ~ factor(year):log(index_true),
                            data = .x)),
         coef_exp = map(mod_exp, coef)) #%>%
  select(opmod, estmod, coef_est, coef_exp) %>%
  unnest(c(coef_est, coef_exp)) #%>%
  mutate(year = seq_len(max_T)) %>%
  pivot_longer(c(coef_est, coef_exp), "delta_type")

biasyear_df %>%
  ggplot(aes(x = year, y = value, color = estmod, linetype = delta_type)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
  facet_grid(estmod ~ ., scales = "free")


ggplot(fullindex_df,
       aes(x = year, y = exp_dev, group = repl, color = estmod)) +
  geom_line() +
  geom_hline(yintercept = 0,  linetype = 2, alpha = 0.5) +
  facet_grid(estmod ~ ., scales = "free")
