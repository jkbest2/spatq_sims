library(spatq)
library(tidyverse)

## survey_sdr <- function(...) {
##   estd <- specify_estimated(...)
##   obj <- make_sim_adfun(1, "combo",
##                         sub_df = tibble(vessel_idx = 2, n = 0),
##                         max_T = 15,
##                         spec_estd = estd)
##   ## fit <- fit_spatq(obj)
##   ## fit <- fit_spatq(obj, fit)
##   fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(maxit = 500))
##   rep <- report_spatq(obj)
##   sdr <- sdreport_spatq(obj)
##   list(obj = obj,
##        rep = rep,
##        sdr = sdr)
## }

## sdr0 <- survey_sdr(omega = list(log_kappa = FALSE, log_tau = FALSE),
##                    lambda = FALSE)
## sdr0k <- survey_sdr(omega = list(log_kappa = TRUE, log_tau = FALSE),
##                     lambda = FALSE)
## sdr0t <- survey_sdr(omega = list(log_kappa = FALSE, log_tau = TRUE),
##                     lambda = FALSE)

sub_df <- tibble(vessel_idx = 2, n = 0)

estd0 <- specify_estimated(omega = list(log_kappa = FALSE, log_tau = FALSE),
                           lambda = FALSE)
estdk <- specify_estimated(omega = list(log_kappa = TRUE, log_tau = FALSE),
                           lambda = FALSE)
estdt <- specify_estimated(omega = list(log_kappa = FALSE, log_tau = TRUE),
                           lambda = FALSE)

res <- cross_df(list(repl = 1:100, spec = list(estd0, estdk, estdt))) %>%
  mutate(obj = map2(repl, spec,
                    ~ make_sim_adfun(.x, "combo", sub_df, spec_estd = .y)),
         fit = map(obj,
                   ~ nlminb(.x$par, .x$fn, .x$gr,
                            control = list(maxit = 1000))),
         sdr = map(obj,
                   ~ sdreport_spatq)) %>%
  select(-obj)

saveRDS(res, "survey_res.Rdata")

obj0 <- make_sim_adfun(53, "combo", sub_df, max_T = 15, spec_estd = estd0)
fit0 <- fit_spatq(obj0)
sdr0 <- sdreport_spatq(obj0)

obj0k <- make_sim_adfun(53, "combo", sub_df, max_T = 15, spec_estd = estdk)
fit0k <- fit_spatq(obj0k)
rep0k <- report_spatq(obj0k)
sdr0k <- sdreport_spatq(obj0k)


estd <- specify_estimated(omega = list(omega_n = list(log_kappa = TRUE, log_tau = FALSE),
                                       omega_w = list(log_kappa = FALSE, log_tau = FALSE)),
                          lambda = FALSE)
obj <- make_sim_adfun(55, "combo", sub_df, max_T = 15, spec_estd = estd)
fit <- fit_spatq(obj)
rep <- report_spatq(obj)
sdr <- sdreport_spatq(obj)

library(VAST)
catch_df <- read_catch(67, "combo")
catch_df <- subsample_catch(catch_df, sub_df)

grid_extrap <- expand.grid(Lat = seq(5, 95, 10),
                           Lon = seq(5, 95, 10),
                           Area_km2 = 100)

vextrap <- make_extrapolation_info(Region = "user",
                                   zone = NA,
                                   input_grid = grid_extrap,
                                   observations_LL = data.frame(Lat = catch_df$s1,
                                                                Lon = catch_df$s2))

vspat <- make_spatial_info(n_x = 450,
                           Lon_i = catch_df$s1,
                           Lat_i = catch_df$s2,
                           Extrapolation_List = vextrap,
                           fine_scale = TRUE)

vdat <- make_data(b_i = catch_df$catch_biomass,
                  a_i = 1,
                  c_iz = 1,
                  t_iz = catch_df$time,
                  FieldConfig = c(Omega1 = 1, Epsilon1 = 0,
                                  Omega2 = 1, Epsilon2 = 0),
                  ObsModel_ez = c(PosDist = 1, Link = 1),
                  spatial_list = vspat, Aniso = FALSE)



## Use observation locations from simulated data set. Only using survey
## information so just the survey grid
obj0 <- make_sim_adfun(53, "combo", sub_df, max_T = 15, spec_estd = estd0)
sim0 <- obj0$simulate()

## Read in simulated data set, replace catch with simulated observations
sim_catch <- read_catch(53, "combo") %>%
  subsample_catch(sub_df) %>%
  filter(time <= 15) %>%
  mutate(catch_biomass = sim0$catch_obs)

mesh <- generate_mesh()
fem <- spatq::generate_fem(mesh)
index_df <- create_index_df(T = 15)
estd <- specify_estimated(omega = TRUE, lambda = FALSE)

sim_data <- prepare_data(catch_df = sim_catch,
                         index_df = index_df,
                         mesh = mesh, fem = fem)
sim_pars <- prepare_pars(sim_data, mesh)
sim_map <- prepare_map(sim_pars, estd)
sim_random <- prepare_random(sim_map)
sim_obj <- prepare_adfun(sim_data, sim_pars, sim_map, sim_random,
                         silent = FALSE)

sim_fit <- fit_spatq(sim_obj)
sim_rep <- report_spatq(sim_obj)
sim_sdr <- sdreport_spatq(sim_obj)
