if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
} else {
  library(spatq)
}
library(tidyverse)

spec_df <- cross_df(list(study = "prefintensity",
                         repl = factor(1:25),
                         opmod = factor(1:6),
                         estmod = c("survey",
                                    "survey_spt"),
                         root_dir = ".")) %>%
  rowwise() %>%
  mutate(spec = list(spatq_simstudyspec(study, repl, opmod, estmod, root_dir)))

index_df <- map_df(spec_df$spec,
                   read_index,
                   estmods = c("survey", "survey_spt"))

index_df %>%
  group_by(estmod, opmod, repl) %>%
  summarize(pdhess = !any(is.na(raw_est))) %>%
  summarize(sum(pdhess))

bias_df <- index_df %>%
  group_by(estmod, opmod) %>%
  nest() %>%
  mutate(bias = map(data, bias_metric)) %>%
  mutate(delta = map_dbl(bias, "delta"))

bias_df %>%
  ungroup() %>%
ggplot(aes(x = as.numeric(opmod), y = delta, color = estmod)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed")

bias_wide <- bias_df %>%
  select(opmod, estmod, delta) %>%
  pivot_wider(names_from = opmod, values_from = delta)


##  res_ls <- map(spec_df$spec, read_rdata)

repl <- 6
pop <- read_popstate("prefintensity", repl, 6, filetype = "h5")
pop <- log(pop[, , 1:15])

beta_hat <- apply(pop, 3, mean)
plot(beta_hat)

p2 <- sweep(pop, 3, beta_hat)
dev.new()
par(mfrow = c(4, 4))
for (yr in 1:15) {
  image(p2[, , yr], asp = 1)
}

omega_hat <- apply(p2, c(1, 2), mean)
dev.new()
image(omega_hat, asp = 1)

epsilon_hat <- p2
dev.new()
par(mfrow = c(4, 4))
for (yr in 1:15) {
  epsilon_hat[, , yr] <- p2[, , yr] - omega_hat
  image(epsilon_hat[, , yr], asp = 1)
}
hist(epsilon_hat)

res_spt <- read_rdata(spatq_simstudyspec("prefintensity", repl, 6, "survey_spt"))
eps <- gather_nvec(res_spt$lpb)$epsilon_n[[1]]
mesh <- generate_mesh()
proj_df <- create_index_df(1, 15)
proj <- generate_projection(mesh, proj_df, group = proj_df$year)
spt <- project_field(eps, proj, c(100, 100, 15))
plot_field(res_spt$lpb, par_regex = "epsilon")

## setup <- spatq_simsetup("prefintensity", 6, 6,
##                         sub_df = data.frame(vessel_idx = 2, n = 0),
##                         root_dir = ".",
##                         max_T = 15,
##                         index_step = 1,
##                         spec_estd = specify_estimated(beta = TRUE,
##                                                       omega = list(omega_n = TRUE,
##                                                                    omega_w = FALSE),
##                                                       epsilon = list(epsilon_n = TRUE,
##                                                                      epsilon_w = FALSE),
##                                                       lambda = FALSE,
##                                                       kappa_map = c(1, NA, 1, NA, NA, NA, NA, NA),
##                                                       obs_lik = 1L))
