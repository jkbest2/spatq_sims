if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
} else {
  library(spatq)
}
library(tidyverse)

get_est_pop <- function(spec) {
  rda <- read_rdata(spec)
  logpopmat <- rda$rep$Ilog_n
  nyr <- length(logpopmat) / 100 / 100
  logpop <- array(rda$rep$Ilog_n, dim = c(100, 100, nyr))
  exp(logpop)
}

## Get true population map and estimated maps for each `estmod`
get_all_pops <- function(spec, estmods = c("survey_spt", "sptemp_ab", "sptemp_q")) {
  pop_true <- list(true = read_popstate(spec$study,
                                        spec$repl,
                                        spec$opmod,
                                        filetype = "h5"))

  est_specs <- lapply(estmods,
                      function(em)
                        spatq_simstudyspec(spec$study,
                                           repl = spec$repl,
                                           opmod = spec$opmod,
                                           estmod = em,
                                           sub_df = NULL,
                                           estd = NULL,
                                           root_dir = "."))
  names(est_specs) <- estmods

  c(pop_true, lapply(est_specs, get_est_pop))
}

pop_comp_layout <- function(npops, nyears, colorbar = TRUE) {
  lay <- matrix(1:(npops * nyears),
                nrow = npops, ncol = nyears,
                byrow = TRUE)
  cbind(lay, max(lay) + 1)
}

rescale_pop <- function(pop) {
  index <- rescale_index(apply(pop, 3, sum))
  scale <- attr(index, "scale")
  pop / scale
}

plot_pop <- function(pop, years = 1:dim(pop)[3], ...) {
  walk(years,
       function(yr) {
         image(z = pop[, , yr],
               asp = 1,
               axes = FALSE,
               ...)
         tot <- sum(pop[, , yr])
         text(x = 1,
              y = -0.1,
              signif(tot, 4),
              pos = 2,
              cex = 1)
       })
}

label_pos <- function(nlabs) {
  offset <- 1 / nlabs / 2
  seq(offset, 1 - offset, length.out = nlabs)
}

plot_pop_comparison <- function(pop_list, years = NULL, rescale = TRUE) {
  if (is.null(years)) {
    year_tot <- map_int(pop_list,
                        function(pop) dim(pop)[3])
    years <- 1:min(year_tot)
  }
  nyears <- length(years)

  pop_names <- names(pop_list)
  npops <- length(pop_list)

  if (rescale) {
    pop_list <- lapply(pop_list, rescale_pop)
  }

  pop_range <- range(pop_list, na.rm = TRUE, finite = TRUE)

  lmat <- pop_comp_layout(npops, nyears, colorbar = TRUE)
  layout(lmat)
  par(oma = c(1, 2, 1, 1), mar = c(1, 1, 1, 1))

  col <- hcl.colors(12)

  walk(pop_list,
       plot_pop,
       years = years,
       zlim = pop_range,
       col = col)
  par(mar = c(3, 1, 3, 5))
  plot_colorbar(pop_range, col)
  mtext(rev(pop_names),
        side = 2,
        at = label_pos(length(pop_names)),
        outer = TRUE)
  mtext(c(paste("Year", years), ""),
        side = 3,
        at = label_pos(length(years) + 1),
        outer = TRUE,
        adj = 1, padj = 1)
}

studies <- c("qdevscaling",
             "prefintensity",
             "densdepq",
             "habq",
             "bycatch")
studies <- "habq"
repls <- 1
opmods <- 1:6
estmods <- c("survey_spt",
             "sptemp_ab",
             "sptemp_q")
spec_ls <- cross(list(study = studies,
                      repl = repls,
                      opmod = opmods))

ppc <- function(spec, estmods) {
  pops <- get_all_pops(spec, estmods)
  fn <- paste0(spec$study, "_",
               str_pad(spec$repl, 2), "_",
               str_pad(spec$opmod, 2), ".png")
  png(file.path("popfigs", fn),
      width = 1920, height = 1080)
  plot_pop_comparison(pops, years = seq(1, 15, 2))
  mtext(paste("Study:", spec$study,
              "Replicate:", spec$repl,
              "Sim value:", spec$opmod),
        side = 1,
        outer = TRUE,
        padj = 0)
  dev.off()
}

if (!dir.exists("popfigs")) {
  dir.create("popfigs")
}
walk(spec_ls, ppc, estmods = estmods)
