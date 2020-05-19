## library(spatq)
devtools::load_all("~/src/spatq", helpers = FALSE, export_all = FALSE, recompile = TRUE)
library(tidyverse)

repl <- 1
sc <- "pref"
estd <- specify_estimated(beta     = TRUE,
                          gamma    = FALSE,
                          omega    = TRUE,
                          epsilon1 = FALSE,
                          lambda   = TRUE,
                          eta      = FALSE,
                          phi      = TRUE,
                          psi1     = FALSE)

obj_1pref <- make_sim_adfun(repl = repl,
                            sc = sc,
                            sub_df = NULL,
                            root_dir = getwd(),
                            max_T = 10,
                            spec_estd = estd,
                            silent = FALSE,
                            runSymbolicAnalysis = TRUE,
                            normalize = TRUE)

fit_1pref <- fit_spatq(obj_1pref)
fit_1pref <- fit_spatq(obj_1pref, fit_1pref)
rep_1pref <- report_spatq(obj_1pref)
sdr_1pref <- sdreport_spatq(obj_1pref)
