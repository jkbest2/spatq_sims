library(tidyverse)
## library(spatq)
## devtools::load_all("~/src/spatq", helpers = FALSE)

list_Rdata <- function(root_dir = "results") {
  file.path(root_dir,
            list.files(root_dir, pattern = "\\.Rdata$"))
}

is_pdhess <- function(rdata) {
  ispd <- readRDS(rdata)$sdr$pdHess
  if (is.null(ispd)) ispd <- FALSE
  return(ispd)
}

convg_df <- map_dfr(list_Rdata(), ~ readRDS(.)$spec) %>%
  select(estmod, opmod, repl, Rdata) %>%
  mutate(convcode = map_dbl(Rdata, ~ readRDS(.)$fit$convergence),
         outer_mgc = map_dbl(Rdata, ~ max(readRDS(.)$fit$grad)),
         pdhess = map_lgl(Rdata, is_pdhess))

convg_df %>%
  group_by(estmod) %>%
  summarize(n = n(),
            npd = sum(pdhess))
