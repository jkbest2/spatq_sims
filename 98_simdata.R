library(tidyverse)
library(INLA)
library(TMB)
devtools::load_all("~/src/spatq", helpers = FALSE)

compile("src/spatq_simplified.cpp")
dyn.load(dynlib("src/spatq_simplified"))

source("97_debug_fns.R")

set.seed(12345)



## Parameter values
n_year <- 10
n_obs <- 2000

beta1 <- rep(0.5, n_year)
beta2 <- rep(-0.5, n_year)
catch_sigma <- 0.75

rho <- 20
sigma <- 1
tau <- pars_tau(sigma, rho)
kappa <- pars_kappa(rho)



## Spatial realizations
mesh <- generate_mesh()
fem2 <- generate_fem(mesh)
fem <- inla.mesh.fem(mesh)
spde <- inla.spde2.matern(mesh,
                          prior.variance.nominal = 1,
                          prior.range.nominal = 80)

q <- tau^2 * (kappa^4 * fem$c0 + 2 * kappa^2 * fem$g1 + fem$g2)

omega_indep <- inla.qsample(2, q)
r <-  matrix(c(1, 0.9, 0.9, 1), nrow = 2)
nw_sd <- c(1, 0.2)
s <- diag(nw_sd) %*% r %*% diag(nw_sd)
Lt_nw <- chol(s)
omega <- omega_indep %*% Lt_nw
omega1 <- as.vector(omega[, 1])
omega2 <- as.vector(omega[, 2])



## Index grid
index_grid <- seq(0.5, 99.5, 1.0)
index_df <- create_index_df(step = 1, T = n_year)
index_proj <- inla.spde.make.A(mesh, as.matrix(index_locs[, c(2, 3)]))



## Simulate catch
catch_locs <- tibble(year = rep(seq(1:n_year), each = n_obs),
                     s1 = 100 * runif(n_year * n_obs),
                     s2 = 100 * runif(n_year * n_obs))

catch_proj <- inla.spde.make.A(mesh, as.matrix(catch_locs[c(2, 3)]))

sim_catch <- function(p1, p2, catch_sigma) {
  p_enc <- plogis(p1)
  catch_meanlog <- p2 - catch_sigma^2 / 2
  rbinom(1, 1, p_enc) * rlnorm(1, catch_meanlog, catch_sigma)
}

catch_df <- tibble(
  year = catch_locs$year,
  s1 = catch_locs$s1,
  s2 = catch_locs$s2,
  beta1 = beta1[year],
  beta2 = beta2[year],
  spat1 = as.vector(catch_proj %*% omega1),
  spat2 = as.vector(catch_proj %*% omega2),
  p1 = beta1 + spat1,
  p2 = beta2 + spat2,
  catch_biomass = map2_dbl(p1, p2, sim_catch, catch_sigma = catch_sigma)
)

catch_df



plot_flag <- FALSE
if (plot_flag) {
  dev.new()
  filled.contour(
    index_grid, index_grid,
    matrix(index_proj %*% omega1, nrow = 100))

  dev.new()
  filled.contour(
    index_grid, index_grid,
    matrix(index_proj %*% omega2, nrow = 100))

  true_index_df <- tibble(
    year = index_locs[, 3],
    s1 = index_locs[, 1],
    s2 = index_locs[, 2],
    beta1 = beta1[year],
    beta2 = beta2[year],
    spat1 = as.vector(index_proj %*% omega1),
    spat2 = as.vector(index_proj %*% omega2),
    proc1 = exp(beta1 + spat1),
    proc2 = exp(beta2 + spat2),
    index = proc1 * proc2)

  dev.new()
  filled.contour(
    index_grid, index_grid,
    matrix(true_index_df$index[10001:20000], nrow = 100))
}



catch_df2 <- catch_df %>%
  transmute(time = year,
            vessel_idx = 1,
            effort = 1,
            catch_biomass = catch_biomass,
            s1 = s1, s2 = s2)

estd <- specify_estimated(beta = TRUE,
                          omega = TRUE,
                          lambda = FALSE)

data <- prepare_data(catch_df2, index_df, mesh, fem2)
parameters <- prepare_pars(data, mesh)
map <- prepare_map(parameters, estd)
random <- prepare_random(map)

data <- simplify_data(data)
parameters <- simplify_pars(parameters)
map <- simplify_map(map)
random <- simplify_random(random)

data$proc_switch <- simple_proc_switch(random)
data$norm_flag <- FALSE

obj <- MakeADFun(
  data = data,
  parameters = parameters,
  map = map,
  random = random,
  silent = FALSE,
  DLL = "spatq_simplified"
)

## TODO Try without normalization: make sure that it's implemented correctly and
## make sure that it's actually faster!
if (!data$norm_flag) {
  obj <- TMB::normalize(obj, flag = "incl_data", value = FALSE)
}
if (length(random > 0)) {
  TMB::runSymbolicAnalysis(obj)
}

fit <- fit_spatq(obj)
fit <- fit_spatq(obj, fit)
rep <- report_spatq(obj)
sdr <- sdreport_spatq(obj)
