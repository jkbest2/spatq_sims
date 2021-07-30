devtools::load_all("~/dev/spatq", helpers = FALSE)


spec0 <- specify_estimated(beta = list(beta_n = TRUE,
                                       beta_w = FALSE),
                           omega = FALSE,
                           epsilon = FALSE,
                           lambda = FALSE,
                           kappa_map = c(NA, NA,   # omega
                                         NA, NA,   # epsilon
                                         NA, NA,   # phi
                                         NA, NA),  # psi
                           obs_lik = 1L)
sub_df <- data.frame(vessel_idx = 2, n = 0)

simspec_df <- cross_df(list(study = "prefintensity",
                            repl = 1:25,
                            opmod = 1:6)) %>%
  rowwise() %>%
  mutate(spec = list(spatq_simsetup(study = study,
                                    repl = repl,
                                    opmod = opmod,
                                    sub_df = sub_df,
                                    max_T = 15,
                                    index_step = 1,
                                    spec_estd = spec0)))

setup0 <- spatq_simsetup(study = "prefintensity",
                         repl = 6,
                         opmod = 6,
                         sub_df = data.frame(vessel_idx = 2,
                                             n = 0),
                         max_T = 10,
                         index_step = 1,
                         spec_estd = specify_estimated(beta = list(beta_n = TRUE,
                                                                   beta_w = FALSE),
                                                       omega = FALSE,
                                                       epsilon = FALSE,
                                                       lambda = FALSE,
                                                       kappa_map = c(NA, NA,   # omega
                                                                     NA, NA,   # epsilon
                                                                     NA, NA,   # phi
                                                                     NA, NA),  # psi
                                                       obs_lik = 1L))

obj0 <- spatq_obj(setup0)
fit0 <- spatq_fit(obj0)
fit0 <- spatq_fit(obj0, fit0)
sdr0 <- sdreport_spatq(obj0)

setup <- spatq_simsetup(study = "prefintensity",
                        repl = 6,
                        opmod = 6,
                        sub_df = data.frame(vessel_idx = 2,
                                            n = 0),
                        max_T = 10,
                        index_step = 1,
                        spec_estd = specify_estimated(beta = list(beta_n = TRUE,
                                                                  beta_w = FALSE),
                                                      omega = list(omega_n = TRUE,
                                                                   omega_w = FALSE),
                                                      epsilon = list(epsilon_n = TRUE,
                                                                     epsilon_w = FALSE),
                                                      lambda = FALSE,
                                                      kappa_map = c(1, NA,    # omega
                                                                    1, NA,    # epsilon
                                                                    NA, NA,   # phi
                                                                    NA, NA),  # psi
                                                      obs_lik = 1L))


obj <- spatq_obj(setup)

sim <- obj$simulate()

fit <- spatq_fit(obj)
fit <- spatq_fit(obj, fit)
fit <- spatq_fit(obj, fit)
fit <- spatq_fit(obj, fit, method = "BFGS")
sdr <- sdreport_spatq(obj)
rep <- report_spatq(obj)

index <- rescale_index(sdr$unbiased$value, sdr$unbiased$sd)

pop <- read_popstate(study = "prefintensity",
                     repl = 6,
                     opmod = 6,
                     filetype = "feather") %>%
  filter(year <= 10) %>%
  mutate(pop_index = rescale_index(pop)$index)
  ##        est_index = index$index,
  ##        est_sd = index$sd)
