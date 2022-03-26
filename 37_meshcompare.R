if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
  ## study <- "qdevscaling"
  ## repls <- 1:5
} else {
  library(spatq)
  ## args <- commandArgs(trailingOnly = TRUE)
  ## study <- args[1]
  ## repls <- args[2]:args[3]
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
  nparfits <- 2L

nblasthreads <- ntotcores %/% nparfits

repls <- 1:5
studies <- c("qdevscaling",
             "habq")
opmods <- c(1, 3, 6)
estmods <- c("survey_spt",
             "sptemp_ab")
mesh_res <- c("coarse", "medium", "fine")
root_dir <- "./meshcompare"

fitspec_ls <- cross(list(estmod = estmods,
                         opmod = opmods,
                         repl = repls,
                         study = studies,
                         mesh_res = mesh_res)) %>%
  map(function(l) list_modify(l, root_dir = file.path(root_dir, l$mesh_res)))

fit <- function(fit_spec) {
  spec <- spatq_simstudyspec(fit_spec)
  max_T <- 15

  ## Suppress warnings about NaNs during fitting
  suppressWarnings({
    setup <- spatq_simsetup(study = fit_spec$study,
                            repl = fit_spec$repl,
                            opmod = fit_spec$opmod,
                            estmod = fit_spec$estmod,
                            max_T = max_T,
                            index_step = 1,
                            mesh_resolution = fit_spec$mesh_resolution)
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
