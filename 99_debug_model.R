library(tidyverse)
library(TMB)
devtools::load_all("../spatq", helpers = FALSE, export_all = FALSE)

source("97_debug_fns.R")



compile("src/spatq_simplified.cpp")
dyn.load(dynlib("src/spatq_simplified"))

## estd <- specify_estimated(omega = list(
##                             omega_n = list(
##                               log_kappa = TRUE,
##                               log_tau = TRUE
##                             ),
##                             omega_w = list(
##                               log_kappa = TRUE,
##                               log_tau = TRUE
##                             )
##                           ), lambda = TRUE)

estd <- specify_estimated(beta = TRUE,
                          gamma = FALSE,
                          omega = FALSE,
                          epsilon1 = list(epsilon1_n = TRUE,
                                          epsilon1_w = FALSE),
                          lambda = TRUE,
                          eta = FALSE,
                          phi = FALSE,
                          psi = FALSE)
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
names(map) <- gsub("epsilon1", "epsilon", names(map))
## ## map$epsilon_n <- factor(rep(NA, 404 * 25))
## map$epsilon_w <- factor(rep(NA, 404 * 25))
random <- simplify_random(random)

data$proc_switch <- simple_proc_switch(random)
data$norm_flag <- TRUE
data$incl_data <- FALSE

parameters$epsilon1_n <- matrix(rnorm(404 * 25, 0, 0.1), nrow = 404, ncol = 25)
parameters$epsilon1_w <- matrix(rnorm(404 * 25, 0, 0.1), nrow = 404, ncol = 25)
names(parameters) <- gsub("epsilon1", "epsilon", names(parameters))
## random[3:4] <- c("epsilon_n", "epsilon_w")
random <- "epsilon_n"

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

fit1 <- fit_spatq(obj)
fit1 <- fit_spatq(obj, fit1)
rep1 <- report_spatq(obj)
sdr1 <- sdreport_spatq(obj)

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
