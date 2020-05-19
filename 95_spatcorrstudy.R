library(tidyverse)
library(INLA)
library(TMB)
devtools::load_all("~/src/spatq", helpers = FALSE)

compile("src/dumb.cpp")
dyn.load(dynlib("src/dumb"))

source("97_debug_fns.R")

corr_omega <- function(omega, rho) {
  r <- matrix(c(1, rho, rho, 1), nrow = 2)
  nw_sd <- c(1, 0.2)
  s <- diag(nw_sd) %*% r %*% diag(nw_sd)
  lt <- chol(s)
  omega %*% lt
}

## Simulate catch using logit-lognormal
sim_catch <- function(p1, p2, catch_sigma) {
  p_enc <- plogis(p1)
  catch_meanlog <- p2 - catch_sigma^2 / 2
  rbinom(1, 1, p_enc) * rlnorm(1, catch_meanlog, catch_sigma)
}

## Spatial prerequisites
mesh <- generate_mesh()
fem <- generate_fem(mesh)

index_grid <- seq(0.5, 99.5, 1)
index_locs <- as.matrix(expand.grid(index_grid, index_grid))
index_proj <- inla.spde.make.A(mesh, index_locs)

n_obs <- 2000
catch_locs <- tibble(s1 = 100 * runif(n_obs),
                     s2 = 100 * runif(n_obs))
catch_proj <- inla.spde.make.A(mesh, as.matrix(catch_locs))

rho <- 20
sigma <- 1
tau <- pars_tau(sigma^2, rho)
kappa <- pars_kappa(rho)

q <- tau^2 * (kappa^4 * fem$M0 + 2 * kappa^2 * fem$M1 + fem$M2)

n_repl = 1
beta1 <- 0.5
beta2 <- -0.5
nw_sd <- c(1, 0.2)
catch_sigma <- 0.75

pars_logit <- list(beta1 = beta1,
                   beta2 = beta2,
                   omega1 = rep(0, mesh$n),
                   omega2 = rep(0, mesh$n),
                   log_kappa = rep(log(kappa), 2),
                   log_tau = log(tau / nw_sd),
                   log_sigma = log(catch_sigma))
pars_pois <- list(beta1 = log(-log(1 - plogis(beta1))),
                  beta2 = beta2 -
                    log(-log(1 - plogis(beta1))) +
                    log(plogis(beta1)),
                  omega1 = rep(0, mesh$n),
                  omega2 = rep(0, mesh$n),
                  log_kappa = rep(log(kappa), 2),
                  log_tau = log(tau / nw_sd),
                  log_sigma = log(catch_sigma))

omega_df <- tibble(repl = seq_len(n_repl),
                   omega_indep = replicate(n_repl,
                                           inla.qsample(2, q),
                                           simplify = FALSE))

sim_df <- cross_df(list(em = c("logit", "pois"),
                        rho = c(0.0, 0.8),
                        repl = seq_len(n_repl))) %>%
  left_join(omega_df, by = "repl") %>%
  mutate(omega = map2(omega_indep, rho, corr_omega),
         p1 = map(omega, ~ .[, 1] + beta1),
         p2 = map(omega, ~ .[, 2] + beta2),
         obs = map2(p1, p2, sim_catch, catch_sigma = catch_sigma),
         data = map2(em, obs, ~ list(obs = .y,
                                     A = catch_proj,
                                     IA = index_proj,
                                     spde = fem,
                                     link_fn = .x == "pois")),
         pars = map(em, ~ if (. == "logit") pars_logit else pars_pois),
         pdhess = NA) # %>%
  ## select(repl, em, rho, data, pars, pdhess)


## Use a for loop here because only want to fit one at a time anyway and can't
## hold all of `obj` in memory anyway
## for (i in seq_len(nrow(sim_df))) {
i <- 1
    obj <- MakeADFun(data = sim_df$data[[i]],
                     parameters = sim_df$pars[[i]],
                     random = c("omega1", "omega2"),
                     DLL = "dumb")
    fit <- optim(obj$par, obj$fn, obj$gr, method = "CG")
    fit <- optim(fit$par, obj$fn, obj$gr, method = "BFGS")
    fit <- optim(fit$par, obj$fn, obj$gr, method = "BFGS")
    sdr <- sdreport(obj)
    sim_df$pdhess[i] <- sdr$pdHess
## }
