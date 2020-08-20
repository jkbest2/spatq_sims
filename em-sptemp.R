library(tidyverse)
library(INLA)
library(TMB)
devtools::load_all("~/src/spatq", helpers = FALSE, compile = TRUE)
library(pbapply) # Include progress bar on sim fits

### Useful functions
replace_with_sim <- function(sim, old) {
  for (nm in names(old)) {
    if (nm %in% names(old)) {
      old[["nm"]] <- sim[["nm"]]
    }
  }
  old
}

objectify_sim <- function(sim, original_sim, ..., replace_pars = FALSE) {
    data <- replace_with_sim(sim, original_sim$data)
    data <- original_sim$data
    data$catch_obs <- sim$catch_obs
    pars <- original_sim$parameters
    if (replace_pars) {
      pars <- replace_with_sim(sim, pars)
    }
    map <- original_sim$map
    rand <- original_sim$random

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
                               omega =
                                 list(omega_n = list(log_kappa = TRUE,
                                                     log_tau = TRUE),
                                      omega_w = list(log_kappa = TRUE,
                                                     log_tau = TRUE)),
                               epsilon =
                                 list(epsilon_n = list(log_kappa = FALSE,
                                                       log_tau = TRUE),
                                      epsilon_w = list(log_kappa = FALSE,
                                                       log_tau = TRUE)),
                               lambda = TRUE, eta = FALSE,
                               phi = FALSE, psi = FALSE)
data <- prepare_data(catch_df, index_df, mesh, fem)
parameters <- prepare_pars(data, mesh)
parameters$beta_n <- rep(0.5, max_T)
parameters$beta_w <- rep(-10, max_T)
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
sim <- obj_sim$simulate()


obj3 <- objectify_sim(sim, original_sim,
                      silent = FALSE,
                      runSymbolicAnalysis = TRUE,
                      normalize = TRUE)
fit3 <- fit_spatq(obj3)
rep3 <- report_spatq(obj3)
sdr3 <- sdreport_spatq(obj3)

### Plot realization if desired
if (FALSE) {
  n_index_steps <- 100 / index_step
  st_n_sim <- array(data$IA_sptemp %*% as.vector(sim$epsilon_n),
                    dim = c(n_index_steps, n_index_steps, max_T))
  st_w_sim <- array(data$IA_sptemp %*% as.vector(sim$epsilon_w),
                    dim = c(n_index_steps, n_index_steps, max_T))

  ## dev.new()
  par(mfrow = c(2, max_T), mar = c(0, 1, 0, 1), oma = c(0, 2, 1, 0))
  apply(st_n_sim, 3, image, asp = 1, axes = FALSE)
  apply(st_w_sim, 3, image, asp = 1, axes = FALSE)
  mtext("Numbers density", side = 2, outer = TRUE, at = 0.75)
  mtext("Weight per group", side = 2, outer = TRUE, at = 0.25)
  ## mtext(1, side = 3, outer = TRUE, at = 0.1, padj = 1)
  ## mtext(2, side = 3, outer = TRUE, at = 0.3, padj = 1)
  ## mtext(3, side = 3, outer = TRUE, at = 0.5, padj = 1)
  ## mtext(4, side = 3, outer = TRUE, at = 0.7, padj = 1)
  ## mtext(5, side = 3, outer = TRUE, at = 0.9, padj = 1)
}
### Set up simulation study
n_repl <- 100
sims <- replicate(n_repl, obj_sim$simulate(), simplify = FALSE)
em_sims <- pblapply(sims, function(sim) {
  obj <- objectify_sim(sim, original_sim,
                       silent = TRUE,
                       runSymbolicAnalysis = TRUE,
                       normalize = TRUE)
  fit <- fit_spatq(obj)
  sdreport_spatq(obj)
})

### Extract results of simulation study
pd_hess <- map_lgl(em_sims, pluck, "pdHess")
fix_pars <- map(em_sims, pluck, "par.fixed")
fix_pars <- matrix(unlist(fix_pars), nrow = length(fix_pars[[1]]))
rand_pars <- map(em_sims, pluck, "par.random")
rand_pars <- matrix(unlist(rand_pars), nrow = length(rand_pars[[1]]))

parnames <- names(em_sims[[1]]$par.fixed)
## parnames[1:20] <- paste(parnames[1:20],
##                         rep(1:10, 2),
##                         sep = "_")
parnames[23:30] <- paste(parnames[23:30],
                         rep(c("omega_n", "omega_w", "phi_n", "phi_w"), 2),
                         sep = "_")

fixpar_df <- map(em_sims, ~ tibble(par = parnames, val = .$par.fixed)) %>%
  bind_rows()

randpar_df <- map2(em_sims, sims,
                   ~ tibble(par = names(.x$par.random),
                            est_val = .x$par.random,
                            gen_val = c(.y$omega_n,
                                        .y$omega_w,
                                        .y$phi_n,
                                        .y$phi_w))) %>%
  bind_rows()

project_randvec <- function(randvec, n = 404) {
  n_proc <- length(randvec) %/% n
  if (length(randvec) %% n != 0)
    stop("randvec length must be multiple of n")
  procs <- vapply(seq_len(n_proc) - 1, # Need zero-indexing here!
                  function(idx) {
                    as.vector(data$IA_spat %*% randvec[idx * n + seq_len(n)])
                  },
                  rep(0.0, nrow(data$IA_spat)))
  as.vector(procs)
}

proc_df <- map2(em_sims, sims,
                ~ tibble(proc = factor(
                           rep(c("spat_n", "spat_w", "catch_n", "catch_w"),
                               each = nrow(data$IA_spat)),
                          levels = c("spat_n", "spat_w", "catch_n", "catch_w")),
                         est_proc = project_randvec(.x$par.random),
                         gen_proc = c(as.vector(data$IA_spat %*% .y$omega_n),
                                      as.vector(data$IA_spat %*% .y$omega_w),
                                      as.vector(data$IA_spat %*% .y$phi_n),
                                      as.vector(data$IA_spat %*% .y$phi_w)))) %>%
  bind_rows()

fixtrue_df <- tibble(par = parnames,
                     val = c(parameters$beta_n,
                             parameters$beta_w,
                             parameters$lambda_n,
                             parameters$lambda_w,
                             parameters$log_kappa[spec_estd$log_kappa],
                             parameters$log_tau[spec_estd$log_tau],
                             parameters$log_sigma))

ggplot(fixpar_df, aes(x = val)) +
  geom_histogram(bins = 20) +
  geom_vline(data = fixtrue_df, aes(xintercept = val), linetype = "dashed") +
  facet_wrap(~ par, scales = "free") +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("fixpar_spatq.png", height = 7, width = 10)

ggplot(randpar_df, aes(x = gen_val, y = est_val)) +
  geom_hex(aes(fill = after_stat(density))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_smooth(method = lm, formula = y ~ x + 0) +
  facet_wrap(~ par, scales = "free") +
  scale_fill_viridis_c(option = "cividis") +
  labs(x = "Generative value",
       y = "Estimated value") +
  theme(aspect.ratio = 1)
ggsave("randpar_spatq.png", height = 7, width = 8)

randpar_df %>%
  mutate(diff = est_val - gen_val) %>%
  group_by(par) %>%
  summarize(diff_mean = mean(diff),
            diff_sd = sd(diff),
            corr = cor(est_val, gen_val)) %>%
  data.frame

ggplot(proc_df, aes(x = gen_proc, y = est_proc)) +
  geom_hex(aes(fill = after_stat(density))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_smooth(method = lm, formula = y ~ x + 0) +
  facet_wrap(~ proc, scales = "free") +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = "Generative value",
       y = "Estimated value") +
  theme(aspect.ratio = 1)
ggsave("proc_spatq.png", height = 7, width = 8)

proc_df %>%
  mutate(diff = est_proc - gen_proc) %>%
  group_by(proc) %>%
  summarize(diff_mean = mean(diff),
            diff_sd = sd(diff),
            corr = cor(est_proc, gen_proc)) %>%
  data.frame

plot_abund_procs <- function(repl, rand_pars, sims) {
  o_n <- rand_pars[seq_len(404), repl]
  o_n_sim <- sims[[repl]]$omega_n
  s_n <- matrix(data$IA_spat[1:10000, ] %*% o_n, nrow = 100)
  s_n_sim <- matrix(data$IA_spat[1:10000, ] %*% o_n_sim, nrow = 100)
  o_w <- rand_pars[404 + seq_len(404), repl]
  s_w <- matrix(data$IA_spat[1:10000, ] %*% o_w, nrow = 100)
  o_w_sim <- sims[[repl]]$omega_w
  s_w_sim <- matrix(data$IA_spat[1:10000, ] %*% o_w_sim, nrow = 100)

  s_n_breaks <- seq(min(c(s_n_sim, s_n)), max(c(s_n_sim, s_n)), length.out = 13)
  s_w_breaks <- seq(min(c(s_w_sim, s_w)), max(c(s_w_sim, s_w)), length.out = 13)

  par(mfrow = c(2, 2), mar = c(1, 0, 0, 1), oma = c(0, 3, 3, 0))
  image(s_n_sim, asp = 1, axes = FALSE, breaks = s_n_breaks)
  image(s_n, asp = 1, axes = FALSE, breaks = s_n_breaks)
  image(s_w_sim, asp = 1, axes = FALSE, breaks = s_w_breaks)
  image(s_w, asp = 1, axes = FALSE, breaks = s_w_breaks)
  title(paste("Abundance processes, replicate", repl), outer = TRUE)
  mtext("Generative", side = 3, outer = TRUE, at = 0.25)
  mtext("Estimated", side = 3, outer = TRUE, at = 0.75)
  mtext("Numbers density", side = 2, outer = TRUE, at = 0.75)
  mtext("Weight per group", side = 2, outer = TRUE, at = 0.25)
}

plot_cbility_procs <- function(repl, rand_pars, sims) {
  p_n <- rand_pars[2 * 404 + seq_len(404), repl]
  p_n_sim <- sims[[repl]]$phi_n
  q_n <- matrix(data$IA_spat[1:10000, ] %*% p_n, nrow = 100)
  q_n_sim <- matrix(data$IA_spat[1:10000, ] %*% p_n_sim, nrow = 100)
  p_w <- rand_pars[3 * 404 + seq_len(404), repl]
  q_w <- matrix(data$IA_spat[1:10000, ] %*% p_w, nrow = 100)
  p_w_sim <- sims[[repl]]$phi_w
  q_w_sim <- matrix(data$IA_spat[1:10000, ] %*% p_w_sim, nrow = 100)

  q_n_breaks <- seq(min(c(q_n_sim, q_n)), max(c(q_n_sim, q_n)), length.out = 13)
  q_w_breaks <- seq(min(c(q_w_sim, q_w)), max(c(q_w_sim, q_w)), length.out = 13)

  par(mfrow = c(2, 2), mar = c(1, 0, 0, 1), oma = c(0, 3, 3, 0))
  image(q_n_sim, asp = 1, axes = FALSE, breaks = q_n_breaks)
  image(q_n, asp = 1, axes = FALSE, breaks = q_n_breaks)
  image(q_w_sim, asp = 1, axes = FALSE, breaks = q_w_breaks)
  image(q_w, asp = 1, axes = FALSE, breaks = q_w_breaks)
  title(paste("Catchability processes, replicate", repl), outer = TRUE) #
  mtext("Generative", side = 3, outer = TRUE, at = 0.25)
  mtext("Estimated", side = 3, outer = TRUE, at = 0.75)
  mtext("Numbers density", side = 2, outer = TRUE, at = 0.75)
  mtext("Weight per group", side = 2, outer = TRUE, at = 0.25)
}

png("figs/abund_procs_%d.png", width = 480, height = 480)
plot_abund_procs(1, rand_pars, sims)
plot_abund_procs(25, rand_pars, sims)
plot_abund_procs(37, rand_pars, sims)
plot_abund_procs(97, rand_pars, sims)
dev.off()

png("figs/cbility_procs_%d.png", width = 480, height = 480)
plot_cbility_procs(1, rand_pars, sims)
plot_cbility_procs(25, rand_pars, sims)
plot_cbility_procs(37, rand_pars, sims)
plot_cbility_procs(97, rand_pars, sims)
dev.off()

index_df <- map2(em_sims, seq_len(n_repl),
                 function(.x, .y) {
                   index <- rescale_index(.x$value)
                   scale <- attr(index, "scale")
                   index_sd <- .x$sd / scale
                   tibble(repl = .y,
                          year = seq_len(max_T),
                          index = index,
                          sd = index_sd)}) %>%
  bind_rows()

hilite_idx <- sample(seq_len(n_repl), 5)
hilite <- rep(0.2, n_repl)
hilite[hilite_idx] <- 1
index_df %>%
  left_join(tibble(repl = seq_len(n_repl),
                   alpha = hilite), by = "repl") %>%
ggplot(aes(x = year, y = index, group = repl, alpha = alpha)) +
  geom_line() +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL, expand = c(0, 0) ) +
  geom_hline(yintercept = 1, linetype = 2) +
  guides(alpha = FALSE)
