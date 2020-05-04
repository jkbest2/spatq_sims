simple_data <- function() {
  c("catch_obs", "area_swept", "X_n", "X_w", "IX_n", "IX_w", "A_spat",
    "A_sptemp", "IA_spat", "IA_sptemp", "Ih", "R_n", "R_w", "spde",
    "proc_switch", "norm_flag", "incl_data")
}

simple_pars <- function() {
  c("beta_n", "beta_w", "omega_n", "omega_w", "epsilon1_n", "epsilon1_w",
    "lambda_n", "lambda_w", "log_kappa", "log_tau", "log_sigma")
}

## Define function to remove unnecessary pieces of model inputs
simplify_data <- function(data) {
  ## Note that proc_switch gets inserted later
  simple_data <- data[simple_data()]
  attr(simple_data, "T") <- attr(data, "T")
  return(simple_data)
}

simplify_pars <- function(pars) {
  pars <- pars[simple_pars()]
  pars$log_kappa <- pars$log_kappa[1:2]
  pars$log_tau <- pars$log_tau[1:2]
  return(pars)
}

simplify_map <- function(map) {
  ## Need to use %in% here because map won't have all of the parameter names
  simple_map <- map[names(map) %in% simple_pars()]
  if ("log_kappa" %in% names(simple_map)) {
    simple_map$log_kappa <- simple_map$log_kappa[1:4]
  }
  if ("log_tau" %in% names(simple_map)) {
    simple_map$log_tau <- simple_map$log_tau[1:4]
  }
  return(simple_map)
}

simplify_random <- function(random) {
  random[random %in% simple_pars()]
}

simple_proc_switch <- function(random) {
  c("omega_n", "omega_w", "epsilon1_n", "epsilon1_w") %in% random
}
