## if (interactive()) {
##   devtools::load_all("~/dev/spatq", helpers = FALSE)
## } else {
library(spatq)
## }
library(tidyverse)

fix_index <- function(spec) {
  path <- index_path(spec, "csv")
  try({
    index_df <- read_csv(path)
    if ("study.x" %in% names(index_df))
      index_df <- rename(index_df, study = study.x)
    if ("study.y" %in% names(index_df))
      index_df <- select(index_df, -study.y)
    write_csv(index_df, path)
  })
}

wl <- cross(list(study = c("qdevscaling",
                           "sharedq",
                           "prefintensity"),
                 repl = 1:25,
                 opmod = 1:6,
                 estmod = c("survey",
                            "spatial_ab",
                            "spatial_q")))

walk(wl, fix_index)
