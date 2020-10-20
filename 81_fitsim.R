devtools::load_all("~/dev/spatq", helpers = FALSE)

prior <- grf_pcprior(rho0 = rep(5, 8),
                     sig0 = rep(0.1, 8),
                     alpha_rho = rep(0.05, 8),
                     alpha_sig = rep(0.04, 8),
                     use_prior = c(TRUE, TRUE))

spec <- specify_estimated()
setup <- spatqsetup_sim(1, "spat", max_T = 5, prior = prior)
obj0 <- prepare_adfun(setup$data, setup$parameters, setup$map, setup$random,
                      runSymbolicAnalysis = FALSE, normalize = TRUE)
fit0 <- optim(obj0$par, obj0$fn, obj0$gr, method = "BFGS")

spec1 <- specify_estimated(beta = TRUE, omega = TRUE,
                           lambda = TRUE)
setup1 <- update_setup(setup, obj0, spec1)
obj1 <- prepare_adfun(setup1$data, setup1$parameters, setup1$map, setup1$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fit1 <- optim(obj1$par, obj1$fn, obj1$gr, method = "BFGS")
hess1 <- hessian_spatq(obj1, fit1)
sdr1 <- sdreport_spatq(obj1, getJointPrecision = FALSE)

spec2 <- specify_estimated(beta = TRUE,
                           omega = TRUE, epsilon = TRUE,
                           lambda = TRUE)
setup2 <- update_setup(setup1, obj1, spec2)
obj2 <- prepare_adfun(setup2$data, setup2$parameters, setup2$map, setup2$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fit2 <- optim(obj2$par, obj2$fn, obj2$gr, method = "BFGS")
hess2 <- hessian_spatq(obj2, fit2)
sdr2 <- sdreport_spatq(obj2, getJointPrecision = FALSE)

spec4 <- specify_estimated(beta = TRUE,
                           omega = TRUE, phi = TRUE,
                           lambda = TRUE)
setup4 <- update_setup(setup1, obj1, spec4)
obj4 <- prepare_adfun(setup4$data, setup4$parameters, setup4$map, setup4$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fit4 <- optim(obj4$par, obj4$fn, obj4$gr, method = "BFGS")
hess4 <- hessian_spatq(obj4, fit4)
sdr4 <- sdreport_spatq(obj4, getJointPrecision = FALSE)


### Test differences with and without prior/penalty applied
prior <- grf_pcprior(rho0 = rep(5, 8),
                     sig0 = rep(0.1, 8),
                     alpha_rho = rep(0.05, 8),
                     alpha_sig = rep(0.05, 8),
                     use_prior = c(TRUE, TRUE))
## Default includes `use_prior = FALSE`
prior0 <- grf_pcprior()
prior1 <- grf_pcprior(rho0 = rep(5, 8),
                      sig0 = rep(0.1, 8),
                      alpha_rho = rep(0.05, 8),
                      alpha_sig = rep(0.05, 8),
                      use_prior = c(TRUE, FALSE))

spec <- specify_estimated(beta = TRUE, omega = TRUE, lambda = TRUE, phi = TRUE)

## limit to 3 years of data for faster running/testing
setup <- spatqsetup_sim(1, "spat", max_T = 3, spec_estd = spec, prior = prior)
setup0 <- spatqsetup_sim(1, "spat", max_T = 3, spec_estd = spec, prior = prior0)
setup1 <- spatqsetup_sim(1, "spat", max_T = 3, spec_estd = spec, prior = prior1)

objp <- prepare_adfun(setup$data, setup$parameters, setup$map, setup$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fitp <- optim(objp$par, objp$fn, objp$gr, method = "BFGS")
repp <- report_spatq(objp)
sdrp <- sdreport_spatq(objp)

obj0 <- prepare_adfun(setup0$data, setup0$parameters, setup0$map, setup0$random,
                      runSymbolicAnalysis = TRUE, normalize = TRUE)
fit0 <- optim(obj0$par, obj0$fn, obj0$gr, method = "BFGS")
rep0 <- report_spatq(obj0)
sdr0 <- sdreport_spatq(obj0)

## Different prior
setup1 <- spatqsetup_sim(1, "spat", max_T = 3, spec_estd = spec, prior = prior1)
objp1 <- prepare_adfun(setup1$data, setup1$parameters, setup1$map,
                       setup1$random, runSymbolicAnalysis = TRUE,
                       normalize = TRUE)
fitp1 <- optim(objp1$par, objp1$fn, objp1$gr, method = "BFGS")
repp1 <- report_spatq(objp1)
sdrp1 <- sdreport_spatq(objp1)
