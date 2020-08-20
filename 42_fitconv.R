library(tidyverse)
## library(spatq)
## devtools::load_all("~/src/spatq", helpers = FALSE)

list_Rdata <- function(root_dir = "results") {
  file.path(root_dir,
            list.files(root_dir, pattern = "\\.Rdata$"))
}

convg_df <- map_dfr(list_Rdata(), ~ readRDS(.)$spec) %>%
  select(estmod, opmod, repl, Rdata) %>%
  mutate(convcode = map_dbl(Rdata, ~ readRDS(.)$fit$convergence),
         outer_mgc = map_dbl(Rdata, ~ attr(readRDS(.)$fit, "mgc")),
         pdhess = map_lgl(Rdata, ~ readRDS(.)$sdr$pdHess))

