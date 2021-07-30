study_descr <- function(study) {
  switch(study,
         qdevscaling = "magnitude of catchability deviations. ",
         sharedq = "degree to which catchability deviations where shared. ",
         prefintensity = "intensity of the fishing location preference. ")
}

par_descr <- function(study) {
  switch(study,
         qdevscaling = "marginal standard deviation of the catchability deviations was ",
         sharedq = "catchability deviations of the survey fleet were equal to those of the commercial fleet, scaled by a factor of ",
         prefintensity = "preference was proportional to abundance to the power ")
}
om_val <- function(study, opmod) {
  qdevvals <- 10 ^ (-3:2)
  sharedqvals <- seq(0, 1, 0.2)
  prefintensvals <- c(0, 1, 2, 4, 8, 16)

  switch(study,
         qdevscaling = qdevvals[opmod],
         sharedq = sharedqvals[opmod],
         prefintensity = prefintensvals[opmod])
}

om_descr <- function(study, repl, opmod) {
  st <- paste0("This is replicate number ",
              repl, ". ",
              "This study varied the ",
              study_descr(study),
              "In this case the ",
              par_descr(study),
              om_val(study, opmod), ".")

}

om_pars <- function(study, opmod) {
  pars <- list(base_q = 0.01,
               logq_sd = 0.05,
               sharedq = 0,
               pref_power = 1,
               densdepq_mult = 0)

  qdevvals <- 10 ^ seq(-3, -0.5, 0.5)
  sharedqvals <- seq(0, 1, 0.2)
  prefintensvals <- c(0, 1, 2, 4, 8, 16)
  densdepqvals <- seq(0, 1.25, 0.25)

  if (study == "qdevscaling") {
    pars$logq_sd <- qdevvals[opmod]
  } else if (study == "sharedq") {
    pars$sharedq <- sharedqvals[opmod]
  } else if (study == "prefintensity") {
    pars$pref_power <- prefintensvals[opmod]
  } else if (study == "densdepq") {
    pars$logq_sd <- 0
    pars$densdepq_mult <- densdepqvals[opmod]
  }

  pars
}

comm_catchability <- function(logq_devs, pars) {
  pars$base_q * exp(logq_devs * pars$logq_sd - pars$logq_sd ^ 2 / 2)
}

survey_catchability <- function(logq_devs, pars) {
  logq_sd <- pars$logq_sd * pars$sharedq
  pars$base_q * exp(logq_devs * logq_sd - logq_sd ^ 2 / 2)
}

comm_target_prob <- function(abund, pars) {
  pref <- abund ^ pars$pref_power
  pref / sum(pref)
}

om_theme <- function() {
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        rect = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(100, "native"))
}
