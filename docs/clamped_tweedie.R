library(tweedie)
library(tidyverse)

pop_trend <- exp(seq(0, -1.5, length = 25))
p <- 1.84
phi <- 1.2

repls <- 1:100
years <- 1:25
trips <- 1:2500

df <- cross_df(list(trip = trips,
                    year = years,
                    repl = repls))


cpue_df <- tibble(repl = 1:100,
                  data = replicate(100,
{
  tibble(year = rep(1:25, each = n)) %>%
    rowwise() %>%
    mutate(pop = pop_trend[year],
           mu = 0.2 * pop,
           catch = rtweedie(1, p, mu, phi),
           clamp = min(pop, catch),
           sat = clamp < catch) %>%
    group_by(year) %>%
    summarize(cpue = mean(catch),
              clamp_cpue = mean(clamp),
              pop = mean(pop),
              psat = mean(sat))
}, simplify = FALSE)) %>%
  unnest(cols = c(data))

mod_clamp <- lm(log(clamp_cpue) ~ log(pop) + factor(repl), data = cpue_df)
mod_cpue <- lm(log(cpue) ~ log(pop) + factor(repl), data = cpue_df)

coef(mod_clamp)[2]
coef(mod_cpue)[2]

ggplot(cpue_df, aes(x = factor(year), y = psat)) +
  geom_boxplot() +
  geom_jitter(width = 0.25)
