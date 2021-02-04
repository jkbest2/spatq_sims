library(spatq)
library(tidyverse)

### Useful functions
replace_with_sim <- function(sim, old) {
  for (nm in names(old)) {
    if (nm %in% names(sim)) {
      old[[nm]] <- sim[[nm]]
    }
  }
  old
}

objectify_sim <- function(sim, kappa_map, original_setup) {
  estd <- specify_estimated(beta = TRUE, gamma = FALSE,
                            omega = TRUE, epsilon = TRUE,
                            lambda = TRUE, eta = FALSE,
                            phi = TRUE, psi = TRUE,
                            kappa_map = kappa_map)

  setup <- original_setup
  setup$data <- replace_with_sim(sim, original_setup$data)
  setup <- update_setup(setup, sim, estd)

  spatq_obj(setup, silent = TRUE)
}

get_repl_str <- function(repl) {
  str_pad(repl, 2, pad = 0)
}

get_parval_csv <- function(repl, kappa_map, root = "misspec_kappa") {
  csvname <- paste0("parval", get_repl_str(repl), 
  file.path(root, "parval.csv")
}

write_parval <- function(df, repl, kappa_map, root = "misspec_kappa") {
  write_csv(df, get_parval_csv(root))
}

make_parval_df <- function(fit, kappa_map, repl, ref_df = ref_df) {
  pars <- gather_nvec(fit$par)
  kappa_est <- exp(pars$log_kappa)[kappa_map]
  tau_est <- exp(pars$log_tau)
  rho_est <- pars_rho(kappa_est)
  sigsp_est <- pars_sig(tau_est, kappa_est)

  tibble(repl = repl,
         nkappa = length(unique(kappa_map)),
         kappa_est = kappa_est,
         tau_est = tau_est,
         rho_est = rho_est,
         sigsp_est = sigsp_est,
         ref_df)
}

n_sims <- 50
max_T <- 5
kappa_maps <- list(twopar = c(1, 2, 1, 2, 1, 2, 1, 2),
                   fourpar = c(1, 2, 1, 2, 3, 4, 3, 4))

## Generative rho values
rho_gen <- c(70, 70, 35, 35, 60, 60, 30, 30)

### Read in one replicate
repl <- 1
sc <- "spat"
root_dir <- "."
sub_df <- data.frame(vessel_idx = 2, n = 1000)
max_T <- 10

estd_gen <- specify_estimated(beta = TRUE, gamma = FALSE,
                              omega = TRUE, epsilon = TRUE,
                              lambda = TRUE, eta = FALSE,
                              phi = TRUE, psi = TRUE,
                              kappa_map = c(1, 2, 3, 4, 5, 6, 7, 8))
setup_gen <- spatq_simsetup(1, "spat", max_T = 10,
                            spec_estd = estd_gen)
## Replace defaults with generative parameter values
setup_gen$parameters$beta_n <- c(0.5, rep(0, max_T - 1))
setup_gen$parameters$beta_w <- c(-10, rep(0, max_T - 1))
setup_gen$parameters$lambda_n <- 0.4
setup_gen$parameters$lambda_w <- 1.6
setup_gen$parameters$log_kappa <- log(pars_kappa(rho_gen))
setup_gen$parameters$log_tau <- log(pars_tau(0.1, rho_gen))
setup_gen$parameters$log_sigma <- log(0.25)
obj_gen <- spatq_obj(setup_gen,
                     runSymbolicAnalysis = TRUE,
                     normalize = TRUE,
                     silent = TRUE)
original_setup <- list(spec_estd = estd_gen,
                     data = setup_gen$data,
                     parameters = setup_gen$parameters,
                     map = setup_gen$map,
                     random = setup_gen$random)

ref_df <- tibble(par = c("omega_n", "omega_w",
                         "epsilon_n", "epsilon_w",
                         "phi_n", "phi_w",
                         "psi_n", "psi_w"),
                 kappa_true = exp(setup_gen$parameters$log_kappa),
                 tau_true = exp(setup_gen$parameters$log_tau),
                 rho_true = rho_gen,
                 spsig_true = pars_sig(tau_true, kappa_true))

sims <- replicate(n_sims, obj_gen$simulate(), simplify = FALSE)

## Run the simulations, write out to a csv file
root <- "misspec_kappa"
if (!dir.exists(root)) dir.create(root)
octl <- spatq_optcontrol(maxopts = 4)
for (repl in seq_along(sims)) {
  sim <- sims[[repl]]
  for (kmap in kappa_maps) {
    obj <- objectify_sim(sim, kmap, original_setup)
    fit <- tryCatch(
    fit_spatq(obj, optcontrol = octl),
    error = function(e) {
      l <- list()
      l$par <- rep(NA, 21 + length(unique(kmap)))
      names(l) <- c(rep("beta_n", max_T), rep("beta_w", max_T),
		    "lambda_n", "lambda_w",
		    rep("log_kappa", length(unique(kmap))),
		    rep("log_tau", 8),
		    "log_sigma")
      l
    })
    df <- make_parval_df(fit, kmap, repl, ref_df)
    append_parval(df, root)
  }
}
