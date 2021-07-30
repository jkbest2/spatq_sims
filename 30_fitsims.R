if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
  study <- "qdevscaling"
  repls <- 1:5
} else {
  library(spatq)
  args <- commandArgs(trailingOnly = TRUE)
  study <- args[1]
  repls <- args[2]:args[3]
}
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(RhpcBLASctl)

## Get number of cores via SLURM if available
ntotcores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(ntotcores))
  ntotcores <- detectCores()

nparfits <- as.numeric(Sys.getenv("SPATQ_PAR_FITS"))
if (is.na(nparfits))
  nparfits <- 7L

nblasthreads <- ntotcores %/% nparfits

opmods <- 1:6
estmods <- c("model",
             "survey",
             "spatial_ab",
             "spatial_q")
root_dir = "."

fitspec_ls <- cross(list(estmod = estmods,
                         opmod = opmods,
                         repl = repls,
                         study = study,
                         root_dir = root_dir))

fit <- function(fit_spec) {
  spec <- spatq_simstudyspec(fit_spec)
  max_T <- 15

  ## Suppress warnings about NaNs during fitting
  suppressWarnings({
    setup <- spatq_simsetup(study = spec$study,
                            repl = spec$repl,
                            opmod = spec$opmod,
                            estmod = spec$estmod,
                            max_T = max_T,
                            index_step = 1)
    obj <- spatq_obj(setup)
    fit <- spatq_fit(obj)
    fit <- spatq_fit(obj, fit = fit, method = "BFGS")
    fit <- spatq_fit(obj, fit = fit)
    sdr <- sdreport_spatq(obj)
    rep <- report_spatq(obj)
    lpb <- gather_nvec(obj$env$last.par.best)
  })

  save_fit(spec, fit, lpb, rep, sdr, root_dir = spec$root_dir)
  save_index(spec, sdr, max_T, feather = TRUE)
  ## Running this function for the side effects, just return a boolean
  TRUE
}

need_fit <- function(fitspec) {
  spec <- spatq_simstudyspec
  !file.exists(index_path(fitspec, "feather"))
}

create_res_dir(study = study,
               repl = repls)

blas_set_num_threads(nblasthreads)

cl <- makeForkCluster(nparfits)
registerDoParallel(cl)
foreach(fitspec = fitspec_ls,
        .combine = c,
        .inorder = FALSE,
        .packages = c("spatq")) %:%
  foreach::when(need_fit(fitspec)) %dopar% {
            fit(fitspec)
        }
stopCluster(cl)
