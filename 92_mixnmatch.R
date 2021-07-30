if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
} else {
  library(spatq)
}
library(tidyverse)

study <- "prefintensity"
repl_years <- 1:15
opmod <- 1

extract_catch <- function(keep_year, study, repl, opmod, sub_df = NULL, root_dir = ".") {
  read_catch(study, repl, opmod, root_dir, feather = TRUE) %>%
    filter(year == keep_year) %>%
    subsample_catch(sub_df)
}

construct_catch <- function(study, repl_years, opmod, sub_df = NULL, root_dir = ".") {
  catch_df <- map_df(seq_along(repl_years),
                     ~ extract_catch(keep_year = .,
                                     study = study,
                                     repl = repl_years[.],
                                     opmod = opmod,
                                     sub_df = sub_df,
                                     root_dir = root_dir))
  attr(catch_df, "T") <- length(unique(catch_df$year))
  catch_df
}

sub_df <- data.frame(vessel_idx = 2, n = 0)


catch_df <- construct_catch(study, repl_years, opmod, sub_df)
spec <- specify_estimated(beta = TRUE, omega = TRUE, epsilon = TRUE, lambda = FALSE)
setup <- spatq_setup(catch_df, spec, index_step = 1)
obj <- spatq_obj(setup)
fit <- spatq_fit(obj)
fit <- spatq_fit(obj, fit, method = "BFGS")
fit <- spatq_fit(obj, fit)

sdr <- sdreport_spatq(obj)
rep <- report_spatq(obj)

spec1 <- specify_estimated(beta = TRUE, omega = TRUE, epsilon = TRUE, lambda = FALSE,
                           obs_lik = 1L)
setup1 <- spatq_setup(catch_df, spec1, index_step = 1)
obj1 <- spatq_obj(setup1)
fit1 <- spatq_fit(obj1)
fit1 <- spatq_fit(obj1, fit1, method = "BFGS")
fit1 <- spatq_fit(obj1, fit1)
sdr1 <- sdreport_spatq(obj1)
rep1 <- report_spatq(obj1)

spec2 <- specify_estimated(beta = TRUE, omega = TRUE, epsilon = TRUE, lambda = FALSE,
                           kappa_map = c(1, NA, 1, NA, NA, NA, NA, NA), obs_lik = 1L)
setup2 <- spatq_setup(catch_df, spec2, index_step = 1)
obj2 <- spatq_obj(setup2)
fit2 <- spatq_fit(obj2)
fit2 <- spatq_fit(obj2, fit2, method = "BFGS")
fit2 <- spatq_fit(obj2, fit2)
sdr2 <- sdreport_spatq(obj2)
rep2 <- report_spatq(obj2)

plot_field(gather_nvec(sdr1$par.random),
           gather_nvec(sdr2$par.random),
           par_regex = "omega",
           colorbar = TRUE)
dev.new()
plot_field(gather_nvec(sdr1$par.random),
           gather_nvec(sdr2$par.random),
           par_regex = "epsilon",
           colorbar = TRUE)

plot_field(sdr$par.random, par_regex = "omega", colorbar = TRUE)

catch_df %>%
  group_by(year) %>%
  summarize(design_index = mean(catch_biomass)) %>%
  mutate(index = rescale_index(design_index)$index)
