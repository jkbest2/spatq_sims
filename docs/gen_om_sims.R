library(tidyverse)

if (!dir.exists("om_sims"))
  dir.create("om_sims")

param_list <- cross(list(study = c("qdevscaling",
                                   "sharedq",
                                   "prefintensity"),
                         repl = 1:5,
                         opmod = 1:5,
                         root_dir = "~/gscratch/spatq_sims"))

render_om <- function(pars) {
  outfile <- file.path("om_sims",
                       paste(pars$study,
                             pars$repl,
                             pars$opmod,
                             "opmod.html",
                             sep = "_"))
  rmarkdown::render("om_sims.Rmd",
                    params = pars,
                    output_file = outfile)
}

map(param_list, render_om)
