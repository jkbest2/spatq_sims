library(tidyverse)
devtools::load_all("~/dev/spatq", helpers = FALSE)
## library(spatq)

get_parbounds <- function(obj) {
  lower <- rep(-Inf, length(obj$par))
  upper <- rep(Inf, length(obj$par))

  if (obj$env$data$obs_lik == 1) {
    ## If Tweedie, bound the shape parameter between 1 and 2 to ensure a
    ## compound Poisson-gamma observation likelihood
    lower[length(lower)] <- 1
    upper[length(upper)] <- 2
  }

  return (list(lower = lower, upper = upper))
}

################################################################################
repl <- 2
scenario <- "spat"
max_T <- 5L

################################################################################

fix_spec <- specify_estimated(beta = TRUE, lambda = TRUE, obs_lik = 1L)
fix_setup <- spatq_simsetup(repl, scenario,
                            max_T = max_T,
                            spec_estd = fix_spec,
                            init_fixef = FALSE)
fix_obj <- spatq_obj(fix_setup)
fix_bounds <- get_parbounds(fix_obj)

fix_fit1 <- spatq_fit(fix_obj)


fix_fit2 <- nlminb(fix_fit1$par, fix_obj$fn, fix_obj$gr,
                   control = list(eval.max = 1000L))
fix_fit2 <- attach_optdiags(fix_fit2, fix_fit1, fix_obj)
fix_fit2b <- nlminb(fix_fit1$par, fix_obj$fn, fix_obj$gr,
                    control = list(eval.max = 1000L),
                    lower = fix_bounds$lower, upper = fix_bounds$upper)
fix_fit2b <- attach_optdiags(fix_fit2b, fix_fit1, fix_obj)

fix_fit3 <- optim(fix_fit1$par, fix_obj$fn, fix_obj$gr,
                  method = "CG",
                  control = list(maxit = 1000L))
fix_fit3 <- attach_optdiags(fix_fit3, fix_fit1, fix_obj)

fix_fit1b <- optim(fix_fit3$par, fix_obj$fn, fix_obj$gr,
                   method = "L-BFGS-B",
                   control = list(maxit = 1000L),
                   lower = fix_bounds$lower, upper = fix_bounds$upper)
fix_fit1b <- attach_optdiags(fix_fit1b, fix_fit1, fix_obj)

fix_sdr <- sdreport_spatq(fix_obj)

################################################################################
surv_spec <- specify_estimated(beta = TRUE, omega = TRUE, obs_lik = 1L)
surv_setup <- spatq_simsetup(repl, scenario,
                             max_T = max_T,
                             data.frame(vessel_idx = 2, n = 0),
                             spec_estd = surv_spec)
surv_obj <- spatq_obj(surv_setup)
surv_bounds <- get_parbounds(surv_obj)

################################################################################
spat_spec <- specify_estimated(beta = TRUE, omega = TRUE, lambda = TRUE,
                               obs_lik = 1L)
spat_setup <- update_setup(fix_setup, fix_obj, spat_spec)
## spat_setup <- spatq_simsetup(2, "spat", max_T = 15L, spec_estd = spat_spec)
spat_obj <- spatq_obj(spat_setup)
spat_bounds <- get_parbounds(spat_obj)

spat_fit1 <- spatq_fit(spat_obj)

spat_fit2 <- nlminb(spat_obj$par, spat_obj$fn, spat_obj$gr,
                   control = list(trace = 1L, eval.max = 1000L))
spat_fit2 <- attach_optdiags(spat_fit2, spat_fit1, spat_obj)
spat_fit2b <- nlminb(spat_fit1$par, spat_obj$fn, spat_obj$gr,
                    control = list(eval.max = 1000L),
                    lower = spat_bounds$lower, upper = spat_bounds$upper)
spat_fit2b <- attach_optdiags(spat_fit2b, spat_fit1, spat_obj)

spat_fit3 <- optim(spat_obj$par, spat_obj$fn, spat_obj$gr,
                  method = "CG",
                  control = list(maxit = 100L))
spat_fit3 <- attach_optdiags(spat_fit3, NULL, spat_obj)

spat_fit1b <- optim(spat_obj$par, spat_obj$fn, spat_obj$gr,
                    method = "L-BFGS-B",
                    control = list(maxit = 100L, trace = 1L),
                    lower = spat_bounds$lower, upper = spat_bounds$upper)
spat_fit1b <- attach_optdiags(spat_fit1b, NULL, spat_obj)

spat_rep <- report_spatq(spat_obj)
spat_sdr <- sdreport_spatq(spat_obj)

spat2_setup <- spatq_simsetup(repl = 1, sc = "combo", spec_estd = spat_spec)
spat2_obj <- spatq_obj(spat2_setup)
spat2_fit <- spatq_fit(spat2_obj)
