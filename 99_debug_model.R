library(tidyverse)
library(TMB)
devtools::load_all("../spatq", helpers = FALSE, export_all = FALSE)

source("97_debug_fns.R")



compile("src/spatq_simplified.cpp")
dyn.load(dynlib("src/spatq_simplified"))

estd <- specify_estimated(omega = list(
                            omega_n = list(
                              log_kappa = TRUE,
                              log_tau = TRUE
                            ),
                            omega_w = list(
                              log_kappa = TRUE,
                              log_tau = TRUE
                            )
                          ), lambda = TRUE)

repl <- 1
scen <- "pref"
## Code from `make_sim_adfun`
catch_df <- read_catch(repl, scen)
truepop_df <- read_popstate(repl, scen)
index_df <- create_index_df(step = 5, T = attr(catch_df, "T"))

mesh <- generate_mesh()
fem <- generate_fem(mesh)

## prepare normal version
data <- prepare_data(catch_df, index_df, mesh, fem)
parameters <- prepare_pars(data, mesh)
map <- prepare_map(parameters, estd)
random <- prepare_random(map)

## simplify
data <- simplify_data(data)
parameters <- simplify_pars(parameters)
map <- simplify_map(map)
random <- simplify_random(random)

data$proc_switch <- simple_proc_switch(random)
data$norm_flag <- TRUE

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

fit3 <- fit_spatq(obj, fit)
fit3 <- fit_spatq(obj, fit3)
rep3 <- report_spatq(obj)
sdr3 <- sdreport_spatq(obj)

index_est <- rescale_index(sdr$value)
index_sd <- sdr$sd / attr(index_est, "scale")
index_df <- tibble(
  type = rep(c("true", "est"), each = 25),
  year = rep(1:25, 2),
  index = c(rescale_index(truepop_df$pop), index_est),
  low = index - c(rep(NA, 25), qnorm(0.975) * index_sd),
  high = index + c(rep(NA, 25), qnorm(0.975) * index_sd)
)

ggplot(index_df, aes(
  x = year, y = index,
  color = type, fill = type,
  ymin = low, ymax = high
)) +
  geom_ribbon(alpha = 0.5) +
  geom_point() +
  geom_line()

Ilog_n <- array(rep$Ilog_n, dim = c(20, 20, 25))
Ilog_w <- array(rep$Ilog_w, dim = c(20, 20, 25))
I_b <- exp(Ilog_n + Ilog_w)

image(Ilog_n[, , 25])
image(Ilog_w[, , 1])
image(I_b[, , 1])
