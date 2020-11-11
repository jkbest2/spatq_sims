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
            pdhess = sum(pdhess))

nonpdhess_df <- convg_df %>%
  filter(!pdhess)

pdhess_df <- convg_df %>%
  filter(pdhess)

new_fixpar_df <- pdhess_df %>%
  group_by(estmod, opmod) %>%
  slice_min(n = 1, outer_mgc) %>%
  mutate(new_fixpar = map(Rdata, ~ readRDS(.)$fit$par)) %>%
  select(estmod, opmod, new_fixpar)


library(spatq)

res <- readRDS(nonpdhess_df$Rdata[1])

refit_spatq <- function(resfile, new_fixpar) {
  setup <- spatq_simsetup(repl = res$spec$repl,
                          sc = res$spec$opmod,
                          sub_df = res$spec$sub_df[[1]],
                          root_dir = ".",
                          max_T = 15,
                          index_step = 1,
                          spec_estd = res$spec$estd[[1]])
  setup <- update_parameters(setup, res$lpb)
  obj <- spatq_obj(setup)
  fit <- res$fit
  fit$par <- new_fixpar
  fit <- fit_spatq(obj, fit)
  sdr <- sdreport_spatq(obj)
  list(fit = fit, sdr = sdr)
}

refit_df <- nonpdhess_df %>%
  ## slice(n = 1) %>%
  left_join(new_fixpar_df, by = c("estmod", "opmod")) %>%
  mutate(refit = map2(Rdata, new_fixpar, refit_spatq))

for (idx in seq_len(nrow(refit_df))) {
  oldres <- readRDS(refit_df$Rdata)
  oldres
}
