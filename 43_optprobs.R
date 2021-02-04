devtools::load_all("~/dev/spatq", helpers = FALSE)
## library(spatq)
library(TMB)
library(tidyverse)

list_Rdata <- function(root_dir = "results") {
  file.path(root_dir,
            list.files(root_dir, pattern = "\\.Rdata$"))
}

list_csvs <- function(root_dir = "results") {
  file.path(root_dir,
            list.files(root_dir, pattern = "\\.csv$"))
}

get_convcode <- function(fitlist) {
  pluck(fitlist, "fit", "convergence")
}

get_mgc <- function(fitlist) {
  pluck(fitlist, "fit", attr_getter("mgc"))
}

fit_list <- map(list_Rdata("results"), readRDS)
sdr_list <- map(fit_list, pluck, "sdr")

spec_df <- map_dfr(fit_list, pluck, "spec") %>%
  mutate(opmod = factor(opmod, levels = c("pref", "spat", "combo")),
    estmod = factor(estmod, levels = c("spatial_q", "spatial_ab", "survey"))) %>%
  select(estmod, opmod, repl) %>%
  mutate(pdhess = map_lgl(fit_list, pluck, "sdr", "pdHess"))

spec_df %>%
  ggplot(aes(x = opmod, y = estmod, fill = pdhess, color = pdhess)) +
  geom_point(size = 10, shape = "square") +
  facet_wrap(~ repl) +
  coord_equal()

spec_df %>%
  group_by(opmod, estmod) %>%
  summarize(pdhess = mean(pdhess)) %>%
  ggplot(aes(x = opmod, y = estmod, fill = pdhess)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c()

spec_df %>%
  filter(estmod == "spatial_q") %>%
  ggplot(aes(x = repl, y = opmod, fill = pdhess, color = pdhess)) +
  geom_tile() +
  coord_equal() +
  scale_x_continuous(breaks = 1:25, labels = 1:25)

nonpd_idx <- which(!spec_df$pdhess)

spatq_pd <- fit_list[[35]]

fixpar_corrplot(spatq_pd$sdr)
plot_field(spatq_pd$sdr$par.random, colorbar = TRUE)

pdidx_df <- spec_df %>%
  mutate(idx = seq_len(nrow(spec_df))) %>%
  group_by(estmod) %>%
  summarize(pd = list(idx[pdhess]),
            nonpd = list(idx[!pdhess]))

##############
npd <- 1
simspec <- spec_df[npd, ] # 31 is the first non-PD spatial_ab fit
spab_spec <- specify_estimated(beta = TRUE, omega = TRUE,
                               obs_lik = 1L)
spab_setup <- spatq_simsetup(repl = simspec$repl,
                             sc = simspec$opmod,
                             max_T = 15,
                             spec_estd = spab_spec)
## spab_setup <- update_setup(spab_setup, fit_list[[npd]]$sdr, spab_spec)
spab_obj <- spatq_obj(spab_setup, normalize = FALSE)

spab_fit <- spatq_fit(spab_obj,
                      control = list(eval.max = 1000L, iter.max = 750L))
## spab_fit2 <- spatq_fit(spab_obj, spab_fit1)
## spab_fit3 <- spatq_fit(spab_obj, spab_fit2)

## spab_fitn <- spab_fit2
## while (spab_fitn$convergence == 1) {
##   spab_fitn <- spatq_fit(spab_obj, spab_fitn)
## }

spab_sdr <- sdreport_spatq(spab_obj)

#############
npd <- 2
simspec <- spec_df[npd, ]
spq_spec <- specify_estimated(beta = TRUE, omega = TRUE, phi = TRUE,
                              ## kappa_map = c(1, NA,
                              ##               NA, NA,
                              ##               1, NA,
                              ##               NA, NA),
                              obs_lik = 1L)
spq_setup <- spatq_simsetup(repl = simspec$repl,
                            sc = simspec$opmod,
                            max_T = 15,
                            spec_estd = spq_spec)
spq_setup <- update_setup(spq_setup, fit_list[[npd]]$sdr, spq_spec)
spq_obj <- spatq_obj(spq_setup)

spq_fit <- spatq_fit(spq_obj,
                     fit_list[[npd]]$fit,
                     control = list(eval.max = 1000L, iter.max = 750L))
spq_fit1 <- spatq_fit(spq_obj, spq_fit)
## spq_fit2 <- spatq_fit(spq_obj, spq_fit1, method = "BFGS")
spq_sdr <- sdreport_spatq(spq_obj)


#===========================================
spq_oldfit <- fit_list[[npd]]$fit
spq_oldfit$par[18] <- log(pars_kappa(200))
spq_longfit <- spatq_fit(spq_obj, spq_oldfit)
spq_longfit1 <- spatq_fit(spq_obj, spq_longfit)


############
############
rhos <- map(fit_list, pluck, "rep", "rho_sp")
rho_omega <- map_dbl(rhos, pluck, 1)
rho_phi

rho_df <- spec_df %>%
  mutate(rho_omega = map_dbl(fit_list, pluck, "rep", "rho_sp", 1),
         rho_phi = ifelse(estmod == "spatial_q",
                           map_dbl(fit_list, pluck, "rep", "rho_sp", 5),
                           NA)) %>%
  mutate(sig_omega = map_dbl(fit_list, pluck, "rep", "sigma_sp", 1),
         sig_phi = ifelse(estmod == "spatial_q",
                          map_dbl(fit_list, pluck, "rep", "sigma_sp", 5),
                          NA)) %>%
  mutate(beta1 = map_dbl(fit_list, pluck, "fit", "par", 1))

rho_df %>%
  ggplot(aes(x = opmod, y = estmod, fill = pdhess, color = pdhess)) +
  geom_tile() +
  facet_wrap(~ repl)

rho_df %>%
  ggplot(aes(x = rho_omega, color = estmod, fill = estmod)) +
  geom_density(alpha = 0.2)
  ## facet_wrap(~ estmod)

rho_df %>%
  filter(estmod == "spatial_q") %>%
  ggplot(aes(x = rho_phi)) +
  geom_histogram() +
  facet_wrap(~ pdhess)

rho_df %>%
  ggplot(aes(x = rho_omega, y = sig_omega, color = pdhess)) +
  geom_point()

rho_df %>%
  ggplot(aes(x = sig_omega, y = beta1, color = pdhess)) +
  geom_point() +
  facet_wrap(~ estmod)

rho_df %>%
  ggplot(aes(x = rho_phi, y = sig_phi, color = pdhess)) +
  geom_point() +
  facet_wrap(~ opmod)

rho_df %>%
  ggplot(aes(x = repl, y = opmod, fill = pdhess, color = pdhess)) +
  geom_tile() +
  facet_grid(estmod ~ .)



pf <- function(fl_idx) {
  plot_field(fit_list[[fl_idx]]$sdr$par.random, colorbar = TRUE)
}
