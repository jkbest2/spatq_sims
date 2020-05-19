library(tidyverse)
library(INLA)
library(TMB)
devtools::load_all("~/src/spatq", helpers = FALSE)

compile("src/dumb.cpp")
dyn.load(dynlib("src/dumb"))

source("97_debug_fns.R")

## Spatial prerequisites
mesh <- generate_mesh()
fem <- generate_fem(mesh)

index_grid <- seq(0.5, 99.5, 1)
index_locs <- as.matrix(expand.grid(index_grid, index_grid))
index_proj <- inla.spde.make.A(mesh, index_locs)

rho <- 20
sigma <- 1
tau <- pars_tau(sigma^2, rho)
kappa <- pars_kappa(rho)

q <- tau^2 * (kappa^4 * fem$M0 + 2 * kappa^2 * fem$M1 + fem$M2)

## Parameters
n_obs <- 2000

catch_locs <- tibble(s1 = 100 * runif(n_obs),
                     s2 = 100 * runif(n_obs))
catch_proj <- inla.spde.make.A(mesh, as.matrix(catch_locs))

sim_catch <- function(p1, p2, catch_sigma) {
  p_enc <- plogis(p1)
  catch_meanlog <- p2 - catch_sigma^2 / 2
  rbinom(1, 1, p_enc) * rlnorm(1, catch_meanlog, catch_sigma)
}

plot_flag <- FALSE
inds1 <- inds2 <- seq(0.5, 99.5, 1)

if (plot_flag) {
  index1 <- matrix(as.vector(beta[1] + index_proj %*% omega1), nrow = 100)
  index2 <- matrix(as.vector(beta[2] + index_proj %*% omega2), nrow = 100)
  dev.new()
  filled.contour(inds1, inds2, plogis(index1), asp = 1)
  filled.contour(inds1, inds2, index1, asp = 1)
  dev.new()
  filled.contour(inds1, inds2, exp(index2), asp = 1)
  filled.contour(inds1, inds2, index2, asp = 1)
}

set.seed(12345)

new_omega_indep <- function(q)
  omega_indep <- inla.qsample(2, q)

new_replicate <- function(omega_indep, rho = 0, pois = FALSE) {
  r <-  matrix(c(1, rho, rho, 1), nrow = 2)
  nw_sd <- c(1, 0.2)
  s <- diag(nw_sd) %*% r %*% diag(nw_sd)
  Lt_nw <- chol(s)
  omega <- omega_indep %*% Lt_nw
  omega1 <- as.vector(omega[, 1])
  omega2 <- as.vector(omega[, 2])

  beta <- c(0.5, -0.5)
  catch_sigma <- 0.75
  spat1 <- as.vector(catch_proj %*% omega1)
  spat2 <- as.vector(catch_proj %*% omega2)

  p1 <- beta[1] + spat1
  p2 <- beta[2] + spat2
  obs <- map2_dbl(p1, p2, sim_catch, catch_sigma)

  res <- list()

  if (!pois) {
  res$data <- list(obs = obs,
                   A = catch_proj,
                   IA = index_proj,
                   spde = fem,
                   link_fn = 0)

  res$pars <- list(beta1 = beta[1],
                   beta2 = beta[2],
                   omega1 = rep_along(omega1, 0),
                   omega2 = rep_along(omega2, 0),
                   log_kappa = rep(log(kappa), 2),
                   log_tau = log(tau * nw_sd),
                   log_sigma = log(catch_sigma))
  } else {
  res$data <- list(obs = obs,
                   A = catch_proj,
                   IA = index_proj,
                   spde = fem,
                   link_fn = 1)
  res$pars <- list(beta1 = log(-log(1 - plogis(beta[1]))),
                   beta2 = beta[2] -
                     log(-log(1 - plogis(beta[1]))) +
                     log(plogis(beta[1])),
                   omega1 = rep_along(omega1, 0),
                   omega2 = rep_along(omega2, 0),
                   log_kappa = rep(log(kappa), 2),
                   log_tau = log(tau * nw_sd),
                   log_sigma = log(catch_sigma))
  }

  res$map <- list()
  res$rand <- c("omega1", "omega2")

  return(res)
}

make_obj <- function(spec, ...) {
  MakeADFun(data = spec$data,
            parameters = spec$pars,
            map = spec$map,
            random = spec$rand,
            ...,
            DLL = "dumb")
}

fit_obj <- function(obj) {
  fit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")
  fit <- optim(fit$par, obj$fn, obj$gr, method = "BFGS")
  sdr <- sdreport(obj, getJointPrecision = TRUE)

  res <- tibble_row(pdhess = sdr$pdHess,
                    joint_cond = kappa(sdr$jointPrecision),
                    fixed_cond = kappa(sdr$cov.fixed),
                    par_fixed = list(sdr$par.fixed),
                    outer_mgc = max(sdr$gradient.fixed),
                    par_random = list(sdr$par.random))
}

res_fn <- function(spec){
  tryCatch({
    obj <- make_obj(spec)
    fit_obj(obj)},
    error = function(e)
      tibble_row(pdhess = NA,
                 joint_cond = NA,
                 fixed_cond = NA,
                 par_fixed = NA,
                 outer_mgc = NA,
                 par_random = NA))
}

n_repl <- 1

set.seed(666)
omega_df <- tibble(repl = seq_len(n_repl),
                   omega_indep = map(repl, ~ new_omega_indep(q)))
sim_df <- cross_df(list(pois = c(FALSE, TRUE),
                        rho = c(0),
                        repl = seq_len(n_repl))) %>%
  left_join(omega_df, by = "repl") %>%
  mutate(spec = pmap(list(omega_indep, rho, pois),
                     ~ new_replicate(..1, ..2, ..3)))

res_df <- map_dfr(sim_df$spec, res_fn) %>%
  mutate(repl = sim_df$repl,
         rho = sim_df$rho,
         pois = sim_df$pois)

summary_df <- res_df %>%
  group_by(pois, rho) %>%
  summarize(nofit  = sum(is.na(pdhess)),
            pdhess = sum(!is.na(pdhess) & pdhess)) %>%
  ungroup() %>%
  mutate(notpd = n_repl - nofit - pdhess,
         link = ifelse(pois, "Poisson", "logit"),
         proc_corr = factor(rho)) %>%
  select(link, proc_corr, nofit, pdhess, notpd)

res_df %>%
  filter(!is.na(pdhess)) %>%
  mutate(link = ifelse(pois, "Poisson", "logit"),
         proc_corr = rho) %>%
  ggplot(aes(x = pdhess, y = fixed_cond, color = proc_corr)) +
  geom_boxplot() +
  scale_y_log10(n.breaks = 10) +
  facet_wrap(~ link) +
  labs(title = "Fixed-effect parameter Hessian condition number",
       x = "Positive definite Hessian", y = "Condition number")

fixpar_df <- res_df %>%
  filter(!is.na(pdhess)) %>%
  mutate(beta1 = map_dbl(par_fixed, pluck, "beta1"),
         beta2 = map_dbl(par_fixed, pluck, "beta2"),
         log_kappa1 = map_dbl(par_fixed, ~ log(exp(pluck(., 3)))),
         log_kappa2 = map_dbl(par_fixed, ~ log(exp(pluck(., 4)))),
         log_tau1 = map_dbl(par_fixed, ~ log(exp(pluck(., 5)))),
         log_tau2 = map_dbl(par_fixed, ~ log(exp(pluck(., 6)))),
         sigma = map_dbl(par_fixed, ~ exp(pluck(., "log_sigma"))))%>%
  select(repl, rho, pois, pdhess, beta1, beta2, log_kappa1, log_kappa2,
         log_tau1, log_tau2, sigma)

logit_truefix <- tribble(~par,          ~ val,
                        "beta1",        0.5,
                        "beta2",       -0.5,
                        "log_tau1",    log(pars_tau(1^2, 20)),
                        "log_tau2",    log(pars_tau(0.2^2, 20)),
                        "log_kappa1",  log(pars_kappa(20)),
                        "log_kappa2",  log(pars_kappa(20)),
                        "sigma",        0.75)

pois_truefix <- logit_truefix
## Convert beta1 and beta2 to Poisson link
pois_truefix$val[1:2] <- c(
  log(-log(1 - plogis(logit_truefix$val[1]))),
  logit_truefix$val[2] - log(-log(1 - plogis(logit_truefix$val[1]))) + log(plogis(logit_truefix$val[1])))

fixpar_df %>%
  gather(key = "par", value = "val", -repl, -rho, -pois, -pdhess) %>%
  filter(!pois) %>%
  mutate(rho = factor(rho)) %>%
ggplot(aes(x = pdhess, y = val, color = rho)) +
  geom_hline(aes(yintercept = val), data = logit_truefix) +
  geom_jitter(width = 0.1, height = 0) +
  facet_wrap(~ par, scales = "free") +
  labs(title = "Logit-link EM", x = "Positive-definite Hessian", y = "",
       color = "Process correlation")
ggsave("figs/logitlink.png", width = 10)

fixpar_df %>%
  gather(key = "par", value = "val", -repl, -rho, -pois, -pdhess) %>%
  filter(pois) %>%
  mutate(rho = factor(rho)) %>%
ggplot(aes(x = pdhess, y = val, color = rho)) +
  geom_hline(aes(yintercept = val), data = pois_truefix) +
  geom_jitter(width = 0.1, height = 0) +
  facet_wrap(~ par, scales = "free") +
  labs(title = "Poisson-link EM", x = "Positive-definite Hessian", y = "",
       color = "Process correlation")
ggsave("figs/poislink.png", width = 10)

## Double check nll calculation in R
## loglik <- function(obs, p_enc, pos_r, catch_sigma) {
##   if (obs == 0) {
##     nll <- -log(1 - p_enc)
##   } else {
##     nll <- -log(p_enc) - dlnorm(obs,
##                                 log(pos_r) - catch_sigma^2 / 2,
##                                 catch_sigma,
##                                 TRUE)
##   }
##   return(nll)
## }

## poislink_p_enc <- function(p1) {
##   return(1 - exp(-exp(p1)))
## }

## poislink_pos_r <- function(p1, p2) {
##   n <- exp(p1)
##   w <- exp(p2)
##   p <- 1 - exp(-n)
##   return(n * w / p)
## }

## nll_df <- tibble(obs = data$obs,
##        p_enc = map_dbl(p1, poislink_p_enc),
##        pos_r = map2_dbl(p1, p2, poislink_pos_r)) %>%
##   mutate(nll = pmap_dbl(., loglik, catch_sigma))

## sum(nll_df$nll)

## obj_check <- MakeADFun(
##   data = data,
##   parameters = pars,
##   DLL = "dumb")
## obj2$fn(obj2$par)
