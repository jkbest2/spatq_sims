devtools::load_all("~/dev/spatq", helpers = FALSE)

spec <- specify_estimated()
setup <- spatq_simsetup(1, "spat", max_T = 10)
obj0 <- prepare_adfun(setup$data, setup$parameters, setup$map, setup$random,
                      runSymbolicAnalysis = FALSE, normalize = TRUE)
## fit0 <- optim(obj0$par, obj0$fn, obj0$gr, method = "BFGS")
fit0 <- fit_spatq(obj0)
sdr0 <- sdreport_spatq(obj0)

spec1 <- specify_estimated(beta = TRUE, omega = TRUE,
                           lambda = TRUE,
                           kappa_map = c(1, 1, NA, NA, NA, NA, NA, NA))
setup1 <- update_setup(setup, obj0, spec1)
obj1 <- prepare_adfun(setup1$data, setup1$parameters, setup1$map, setup1$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
## fit1 <- optim(obj1$par, obj1$fn, obj1$gr, method = "BFGS")
fit1 <- fit_spatq(obj1)
rep1 <- report_spatq(obj1)
hess1 <- hessian_spatq(obj1, fit1)
sdr1 <- sdreport_spatq(obj1, getJointPrecision = FALSE)

spec2 <- specify_estimated(beta = TRUE,
                           omega = TRUE, epsilon = TRUE,
                           lambda = TRUE,
                           kappa_map = c(1, 2, 1, 2, NA, NA, NA, NA))
setup2 <- update_setup(setup1, obj1, spec2)
obj2 <- prepare_adfun(setup2$data, setup2$parameters, setup2$map, setup2$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fit2 <- optim(obj2$par, obj2$fn, obj2$gr, method = "BFGS")
hess2 <- hessian_spatq(obj2, fit2)
sdr2 <- sdreport_spatq(obj2, getJointPrecision = FALSE)

spec3 <- specify_estimated(beta = TRUE,
                           omega = TRUE, epsilon = TRUE,
                           lambda = TRUE,
                           phi = TRUE,
                           kappa_map = c(1, 2, 1, 2, 1, 2, NA, NA))
setup3 <- update_setup(setup1, obj1, spec3)
obj3 <- prepare_adfun(setup3$data, setup3$parameters, setup3$map, setup3$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fit3 <- optim(obj3$par, obj3$fn, obj3$gr, method = "BFGS")
hess3 <- hessian_spatq(obj3, fit3)
sdr3 <- sdreport_spatq(obj3, getJointPrecision = FALSE)

spec4 <- specify_estimated(beta = TRUE,
                           omega = TRUE, epsilon = TRUE,
                           lambda = TRUE,
                           phi = TRUE, psi = TRUE,
                           kappa_map = c(1, 2, 1, 2, 1, 2, 1, 2))
setup4 <- update_setup(setup1, obj1, spec4)
obj4 <- prepare_adfun(setup4$data, setup4$parameters, setup4$map, setup4$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fit4 <- optim(obj4$par, obj4$fn, obj4$gr, method = "BFGS")
hess4 <- hessian_spatq(obj4, fit4)
sdr4 <- sdreport_spatq(obj4, getJointPrecision = FALSE)

spec5 <- specify_estimated(beta = TRUE,
                           omega = TRUE, epsilon = TRUE,
                           lambda = TRUE,
                           phi = TRUE, psi = TRUE,
                           kappa_map = c(1, 2, 1, 2, 3, 4, 3, 4))
setup5 <- update_setup(setup1, obj1, spec5)
obj5 <- prepare_adfun(setup5$data, setup5$parameters, setup5$map, setup5$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fit5 <- optim(obj5$par, obj5$fn, obj5$gr, method = "BFGS")
## hess5 <- hessian_spatq(obj5, fit5)
sdr5 <- sdreport_spatq(obj5, getJointPrecision = FALSE)

### Try updating sims with changing parameter dimensions
## setup <- spatq_simsetup(1, "spat",
##                         sub_df = data.frame(vessel_idx = 2, n = 0),
##                         max_T = 10, spec_estd = specify_estimated(lambda = FALSE))
## obj_om <- prepare_adfun(setup_om$data,
##                         setup_om$parameters,
##                         setup_om$map,
##                         setup_om$random,
##                         runSymbolicAnalysis = FALSE,
##                         normalize = TRUE)
## fit_om <- fit_spatq(obj_om, optcontrol = spatq_optcontrol(maxopts = 3))
## sdr_om <- sdreport_spatq(obj_om)

optctl <- spatq_optcontrol(maxopts = 3)
max_T <- 10
## OM1 Survey only
spec1 <- specify_estimated(beta = TRUE, gamma = FALSE,
                           omega = TRUE, epsilon = FALSE,
                           lambda = FALSE, eta = FALSE,
                           phi = FALSE, psi = FALSE,
                           kappa_map = c(1, 2, NA, NA, NA, NA, NA, NA))
setup1 <- spatq_simsetup(1, "spat",
                         sub_df = data.frame(vessel_idx = 2, n = 0),
                         max_T = max_T,
                         spec_estd = spec1)
obj1 <- spatq_obj(setup1, silent = TRUE)
fit1 <- fit_spatq(obj1, optcontrol = optctl)
sdr1 <- sdreport_spatq(obj1)

## Survey and commercial
spec2 <- specify_estimated(beta = TRUE, gamma = FALSE,
                           omega = TRUE, epsilon = FALSE,
                           lambda = TRUE, eta = FALSE,
                           phi = FALSE, psi = FALSE,
                           kappa_map = c(1, 2, NA, NA, NA, NA, NA, NA))
setup2 <- spatq_simsetup(1, "spat",
                         max_T = max_T,
                         spec_estd = spec2)
## Update to use values from survey-only model as initial values.
setup2 <- update_setup(setup2, fit1, spec2)
obj2 <- spatq_obj(setup2, silent = TRUE)
fit2 <- fit_spatq(obj2, optcontrol = optctl)
sdr2 <- sdreport_spatq(obj2)

## Survey and commercial with spatial q
spec3 <- specify_estimated(beta = TRUE, gamma = FALSE,
                           omega = TRUE, epsilon = FALSE,
                           lambda = TRUE, eta = FALSE,
                           phi = TRUE, psi = FALSE,
                           kappa_map = c(1, 2, NA, NA, 1, 2, NA, NA))
setup3 <- spatq_simsetup(1, "spat",
                         max_T = max_T,
                         spec_estd = spec3)
## Update to use values from survey-only model as initial values.
setup3 <- update_setup(setup3, fit2, spec3)
obj3 <- spatq_obj(setup3, silent = TRUE)
fit3 <- fit_spatq(obj3, optcontrol = optctl)
sdr3 <- sdreport_spatq(obj3)



## ### Test differences with and without prior/penalty applied
## prior <- grf_pcprior(rho0 = rep(5, 8),
##                      sig0 = rep(0.1, 8),
##                      alpha_rho = rep(0.05, 8),
##                      alpha_sig = rep(0.05, 8),
##                      use_prior = c(TRUE, TRUE))
## ## Default includes `use_prior = FALSE`
## prior0 <- grf_pcprior()
## prior1 <- grf_pcprior(rho0 = rep(5, 8),
##                       sig0 = rep(0.1, 8),
##                       alpha_rho = rep(0.05, 8),
##                       alpha_sig = rep(0.05, 8),
##                       use_prior = c(TRUE, FALSE))

## spec <- specify_estimated(beta = TRUE, omega = TRUE, lambda = TRUE, phi = TRUE)

## ## limit to 3 years of data for faster running/testing
## setup <- spatqsetup_sim(1, "spat", max_T = 3, spec_estd = spec, prior = prior)
## setup0 <- spatqsetup_sim(1, "spat", max_T = 3, spec_estd = spec, prior = prior0)
## setup1 <- spatqsetup_sim(1, "spat", max_T = 3, spec_estd = spec, prior = prior1)

## objp <- prepare_adfun(setup$data, setup$parameters, setup$map, setup$random,
##                       runSymbolicAnalysis = TRUE, normalize = TRUE)
## fitp <- optim(objp$par, objp$fn, objp$gr, method = "BFGS")
## repp <- report_spatq(objp)
## sdrp <- sdreport_spatq(objp)

## obj0 <- prepare_adfun(setup0$data, setup0$parameters, setup0$map, setup0$random,
##                       runSymbolicAnalysis = TRUE, normalize = TRUE)
## fit0 <- optim(obj0$par, obj0$fn, obj0$gr, method = "BFGS")
## rep0 <- report_spatq(obj0)
## sdr0 <- sdreport_spatq(obj0)

## ## Different prior
## setup1 <- spatqsetup_sim(1, "spat", max_T = 3, spec_estd = spec, prior = prior1)
## objp1 <- prepare_adfun(setup1$data, setup1$parameters, setup1$map,
##                        setup1$random, runSymbolicAnalysis = TRUE,
##                        normalize = TRUE)
## fitp1 <- optim(objp1$par, objp1$fn, objp1$gr, method = "BFGS")
## repp1 <- report_spatq(objp1)
## sdrp1 <- sdreport_spatq(objp1)
