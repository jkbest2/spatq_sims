## devtools::load_all("~/src/spatq", helpers = FALSE)
library(spatq)

spec <- specify_estimated(beta = TRUE, gamma = FALSE,
                          omega = list(omega_n = TRUE, omega_w = FALSE),
                          epsilon = FALSE,
                          lambda = TRUE, eta = FALSE,
                          phi = list(phi_n = TRUE, phi_w = FALSE),
                          psi = FALSE,
                          kappa_map = c(1, NA, NA, NA, 1, NA, NA, NA))

setup <- spatq_simsetup(1, "spat", max_T = 15, spec_estd = spec)




c_df <- read_catch(1, "spat")
cdata <- prepare_data()

estd <- specify_estimated(beta = TRUE, gamma = FALSE,
                          omega = FALSE, epsilon = FALSE,
                          lambda = TRUE, eta = FALSE,
                          phi = FALSE, psi = FALSE)

obj0 <- make_sim_adfun(1, "spat", max_T = 10, spec_estd = estd, runSymbolicAnalysis = FALSE)
fit0 <-

fit <- fit_spatq(obj)
rep <- report_spatq(obj)
hess <- hessian_spatq(obj, fit)
sdr <- sdreport_spatq(obj, getJointPrecision = FALSE)
