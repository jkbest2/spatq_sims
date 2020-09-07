library(tidyverse)
library(INLA)
library(TMB)
devtools::load_all("~/src/spatq", helpers = FALSE, compile = TRUE)
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
max_T <- 10

## Copied from ~make_sim_adfun~
catch_df <- read_catch(repl, sc, root_dir)
if (!is.null(max_T)) {
  catch_df <- dplyr::filter(catch_df, time <= max_T)
  attr(catch_df, "T") <- max_T
}
## Subset observations
catch_df <- subsample_catch(catch_df, sub_df)

## Create index integration reference
index_step <- 10
index_df <- create_index_df(step = index_step, T = attr(catch_df, "T"))

## Discretize space
mesh <- generate_mesh()
fem <- generate_fem(mesh)

## Prepare model specification components and set parameter values to simulate
## from
spec_estd <- specify_estimated(beta = TRUE, gamma = FALSE,
                               omega = TRUE,
                               epsilon = list(epsilon_n = TRUE,
                                              epsilon_w = TRUE),
                               lambda = TRUE, eta = FALSE,
                               phi = FALSE, psi = FALSE)

data <- prepare_data(catch_df, index_df, mesh, fem)
parameters <- prepare_pars(data, mesh)
## parameters$beta_n <- rep(0.5, max_T)
## parameters$beta_w <- rep(-10, max_T)
parameters$beta_n <- c(0.5, rep(0, max_T - 1))
parameters$beta_w <- c(-10, rep(0, max_T - 1))
parameters$lambda_n <- 0.4
parameters$lambda_w <- 1.6
## Give spatial parameters a 60-unit range and spatiotemporal a 30-unit range
parameters$log_kappa <- c(log(pars_kappa(60)), log(pars_kappa(60)),
                        log(pars_kappa(30)), log(pars_kappa(30)),
                        log(pars_kappa(60)), log(pars_kappa(60)),
                        log(pars_kappa(30)), log(pars_kappa(30)))
parameters$log_tau <- rep(3, 8)
parameters$log_sigma <- 0
map <- prepare_map(parameters,
                    spec = spec_estd)
random <- prepare_random(map)
original_sim <- list(spec_estd = spec_estd,
                    data = data,
                    parameters = parameters,
                    map = map,
                    random = random)
obj_sim <- prepare_adfun(data = data,
                        parameters = parameters,
                        map = map,
                        random = random,
                        silent = FALSE,
                        runSymbolicAnalysis = TRUE,
                        normalize = FALSE)
gen <- obj_sim$simulate()

### Set up simulation study
n_repl <- 5
sims <- replicate(n_repl, obj_sim$simulate(), simplify = FALSE)
em_sims <- pblapply(sims, function(sim) {
  obj <- objectify_sim(sim, original_sim,
                       replace_pars = TRUE,
                       silent = TRUE,
                       runSymbolicAnalysis = TRUE,
                       normalize = TRUE)
  fit <- fit_spatq(obj)
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
  plot_field(rand_pars[[repl]], sims[[repl]], par_regex = par_regex, colorbar = FALSE)
}

plot_sim(1, rand_pars, sims, "omega")
dev.new()
plot_sim(2, rand_pars, sims, "omega")
dev.new()
plot_sim(3, rand_pars, sims, "omega")
dev.new()
plot_sim(4, rand_pars, sims, "omega")
dev.new()
plot_sim(5, rand_pars, sims, "omega")

plot_sim(1, rand_pars, sims, "epsilon")
dev.new()
plot_sim(2, rand_pars, sims, "epsilon")
dev.new()
plot_sim(3, rand_pars, sims, "epsilon")
dev.new()
plot_sim(4, rand_pars, sims, "epsilon")
dev.new()
plot_sim(5, rand_pars, sims, "epsilon")

##' Can recover parameter values reasonably well if you start near them.
##' Argument for phasing?
df_ref <- tibble(par = make.unique(names(obj_sim$par)),
                 val = obj_sim$par)
map(fix_pars, function(fp) tibble(par = make.unique(names(fp)), val = fp)) %>%
  bind_rows() %>%
  ggplot(aes(x = val)) +
  facet_wrap(~ par, scales = 'free') +
  geom_dotplot() +
  geom_vline(aes(xintercept = val), data = df_ref)
