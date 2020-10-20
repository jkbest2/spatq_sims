library(tidyverse)
library(INLA)
library(TMB)
devtools::load_all("~/dev/spatq", helpers = FALSE, compile = TRUE)
library(pbapply) # Include progress bar on sim fits

### Useful functions
replace_with_sim <- function(sim, old) {
  for (nm in names(old)) {
    if (nm %in% names(sim)) {
      old[[nm]] <- sim[[nm]]
    }
  }
  old
}

objectify_sim <- function(sim, original_sim, ..., replace_pars = FALSE, rand = NULL) {
    data <- replace_with_sim(sim, original_sim$data)
    ## data <- original_sim$data
    ## data$catch_obs <- sim$catch_obs
    pars <- original_sim$parameters
    if (replace_pars) {
      pars <- replace_with_sim(sim, pars)
    }
    map <- original_sim$map
    if (is.null(rand)) {
      rand <- original_sim$random
    }

    prepare_adfun(data = data,
                  parameters = pars,
                  map = map,
                  random = rand, ...)
}

### Read in one replicate
repl <- 1
sc <- "spat"
root_dir <- "."
sub_df <- data.frame(vessel_idx = 2, n = 1000)
max_T <- 5

estd_sim <- specify_estimated(beta = TRUE, gamma = FALSE,
                              omega = TRUE, epsilon = TRUE,
                              lambda = TRUE, eta = FALSE,
                              phi = TRUE, psi = FALSE,
                              kappa_map = factor(c(1, 2, 1, 2, 1, 2, NA, NA)))
prior <- grf_pcprior(rho0 = rep(10, 8),
                     alpha_rho = rep(0.05, 8),
                     use_prior = c(TRUE, FALSE))
setup_sim <- spatqsetup_sim(repl, sc, sub_df = sub_df,
                            max_T = max_T, index_step = 100,
                            spec_estd = estd_sim,
                            prior = prior)

## Initialize like this for Helmert contrasts
setup_sim$parameters$beta_n <- c(0.5, rep(0, max_T - 1))
setup_sim$parameters$beta_w <- c(-10, rep(0, max_T - 1))
setup_sim$parameters$lambda_n <- 0.4
setup_sim$parameters$lambda_w <- 1.6
## Give spatial parameters a 60-unit range and spatiotemporal a 30-unit range
setup_sim$parameters$log_kappa <- c(log(pars_kappa(60)), log(pars_kappa(60)),
                                    log(pars_kappa(30)), log(pars_kappa(30)),
                                    log(pars_kappa(60)), log(pars_kappa(60)),
                                    log(pars_kappa(30)), log(pars_kappa(30)))
setup_sim$parameters$log_tau <- rep(3, 8)
setup_sim$parameters$log_sigma <- 0

original_sim <- list(spec_estd = estd_sim,
                     data = setup_sim$data,
                     parameters = setup_sim$parameters,
                     map = setup_sim$map,
                     random = setup_sim$random)
obj_sim <- prepare_adfun(data = setup_sim$data,
                        parameters = setup_sim$parameters,
                        map = setup_sim$map,
                        random = setup_sim$random,
                        silent = FALSE,
                        runSymbolicAnalysis = TRUE,
                        normalize = FALSE)
gen <- obj_sim$simulate()

## setup2 <- setup_sim
## setup2$data$catch_obs <- gen$catch_obs
## obj2 <- prepare_adfun(setup2$data, setup2$parameters, setup2$map, setup2$random, normalize = TRUE)
## fit2 <- optim(obj2$par, obj2$fn, obj2$gr, method = "BFGS")
## rep2 <- report_spatq(obj2)
## sdr2 <- sdreport_spatq(obj2, bias.correct = FALSE, getJointPrecision = FALSE)

### Set up simulation study
n_repl <- 10
sims <- replicate(n_repl, obj_sim$simulate(), simplify = FALSE)
em_sims <- pblapply(sims, function(sim) {
  obj <- objectify_sim(sim, original_sim,
                       replace_pars = TRUE,
                       silent = TRUE,
                       runSymbolicAnalysis = TRUE,
                       normalize = TRUE)
  ## fit <- fit_spatq(obj)
  fit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")
  rep <- report_spatq(obj)
  sdr <- sdreport_spatq(obj)
  attr(sdr, "fit_mgc") <- attr(fit, "mgc")
  sdr
})

### Extract results of simulation study
pd_hess <- map_lgl(em_sims, pluck, "pdHess")
fix_pars <- map(em_sims, pluck, "par.fixed")
rand_pars <- map(em_sims, pluck, "par.random")

plot_sim <- function(repl, rand_pars, sims, par_regex) {
  plot_field(rand_pars[[repl]], sims[[repl]],
             par_regex = par_regex, colorbar = FALSE)
}

plot_all <- function(repl, rand_pars, sims) {
  for (par in c("omega", "epsilon", "phi", "^psi")) {
    dev.new()
    plot_sim(repl, rand_pars, sims, par)
  }
}

plot_all(5, rand_pars, sims)


## Can recover parameter values reasonably well if you start near them.
## Argument for phasing?
n_par <- length(em_sims[[1]]$par.fixed)
df_ref <- tibble(par = make.unique(names(obj_sim$par)),
                 val = obj_sim$par)
df_fp <- map(fix_pars, ~ tibble(par = make.unique(names(.x)), val = .x)) %>%
  bind_rows() %>%
  mutate(repl = rep(1:10, each = n_par)) %>%
  left_join(tibble(repl = 1:10, pd_hess = pd_hess), by = "repl")

df_fp %>%
  left_join(df_ref, by = "par", suffix = c("", "_ref")) %>%
  group_by(par) %>%
  summarize(absdiff = min(abs(val - val_ref)),
            nomov = any(absdiff == 0)) %>% View
            ## .groups = "keep")

df_fp %>%
  ggplot(aes(x = val, fill = pd_hess, color = pd_hess)) +
  facet_wrap(~ par, scales = 'free') +
  geom_dotplot(y = 0, width = 10, method = "histodot", stackgroups = TRUE) +
  geom_vline(aes(xintercept = val), data = df_ref) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

df_fp %>% mutate(repl = rep(1:10, each = n_par)) %>%
  pivot_wider(names_from = par, values_from = val) %>%
  select(repl, pd_hess, starts_with("log_kappa")) %>%
  View
df_fp %>% mutate(repl = rep(1:10, each = n_par)) %>%
  left_join(tibble(repl = 1:10, pd_hess = pd_hess), by = "repl") %>%
  pivot_wider(names_from = par, values_from = val) %>%
  select(repl, pd_hess, starts_with("log_tau")) %>%
  View


## First replicate; Hessian NOT PD
obj1 <- objectify_sim(sims[[1]], original_sim,
                      replace_pars = TRUE,
                      silent = TRUE,
                      runSymbolicAnalysis = TRUE,
                      normalize = TRUE)
obj1_he <- optimHess(em_sims[[1]]$par.fixed, obj1$fn, obj1$gr)
obj1_ne <- nulleigs(obj1_he)

tau_idx <- 21:28
kappa_idx <- 13:20
obj1_sig2 <- map2_dbl(exp(em_sims[[1]]$par.fixed[tau_idx]),
                      exp(em_sims[[1]]$par.fixed[kappa_idx]),
                      pars_sig2)
obj1_rho <- map_dbl(exp(em_sims[[1]]$par.fixed[kappa_idx]),
                    pars_rho)


## Second replicate; has PD Hessian
obj2 <- objectify_sim(sims[[2]], original_sim,
                      replace_pars = TRUE,
                      silent = TRUE,
                      runSymbolicAnalysis = TRUE,
                      normalize = TRUE)
obj2_he <- optimHess(em_sims[[2]]$par.fixed, obj2$fn, obj2$gr)
obj2_ne <- nulleigs(obj2_he)

tau_idx <- 21:28
kappa_idx <- 13:20
obj2_sig2 <- map2_dbl(exp(em_sims[[2]]$par.fixed[tau_idx]),
                      exp(em_sims[[2]]$par.fixed[kappa_idx]),
                      pars_sig2)
obj2_rho <- map_dbl(exp(em_sims[[2]]$par.fixed[kappa_idx]),
                    pars_rho)

sim_he <- map()

em_hess <- map2(sims, em_sims,
                function(sim, sdr) {
                  obj <- objectify_sim(sim, original_sim,
                                       replace_pars = TRUE,
                                       silent = TRUE,
                                       runSymbolicAnalysis = TRUE,
                                       normalize = TRUE)
                  obj2_he <- optimHess(sdr$par.fixed, obj$fn, obj$gr)
                })
em_ne <- map(em_hess, nulleigs, abstol = 0)
saveRDS(em_hess, "oct-update-hess.RData")

plot()
