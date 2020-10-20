devtools::load_all("~/src/spatq", helpers = FALSE)

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
