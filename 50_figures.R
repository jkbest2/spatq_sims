library(tidyverse)

Iest <- rescale_index(sdr_sim$value)
Isd <- sdr_sim$sd / attr(Iest, "scale")

Itrue <- rescale_index(read_popstate(2, "pref")$pop)

I_df <- tibble(year = rep(1:25, 2),
               type = rep(c("est", "true"), each = 25),
               index = c(Iest, Itrue))
I_sddf <- tibble(year = 1:25,
                 type = rep("est", 25),
                 index = Iest,
                 sd = Isd,
                 low_ci = Iest + qnorm(0.025) * sd,
                 high_ci = Iest + qnorm(0.975) * sd)

ggplot(I_df, aes(x = year, y = index, color = type)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci), data = I_sddf)
