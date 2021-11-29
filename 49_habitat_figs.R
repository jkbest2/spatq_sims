if (interactive()) {
  devtools::load_all("~/dev/spatq", helpers = FALSE)
} else {
  library(spatq)
}
library(hdf5r)
library(tidyverse)
library(patchwork)

## qdevscaling example
qd_pref_fn <- function(h) exp(-(h + 5)^2 / 25)
h5 <- h5file("prep.h5", mode = "r")
qd_hab <- h5[["habitat"]][, , 1]
h5close(h5)
qd_pref <- qd_pref(qd_hab)
qd_pop <- read_popstate("qdevscaling", 1, 1, ".", "h5")[, , 1]
qd_df <- cross_df(list(x = seq(0.5, 99.5, 1),
                       y = seq(0.5, 99.5, 1))) %>%
  mutate(hab = as.vector(qd_hab),
         pref = qd_pref_fn(hab),
         pop = as.vector(qd_pop))

svg("om_results/qdevscaling_hab.svg",
    height = 12, width = 8)
hab_col = hcl.colors(12, "YlOrRd")
pref_col = hcl.colors(12, "cividis")
pop_col = hcl.colors(12, "viridis")
layout(matrix(1:6, byrow = TRUE, ncol = 2),
       widths = c(6, 1))
par(oma = rep(1, 4),
    mar = rep(1, 4))
image(qd_hab,
      main = "Habitat",
      asp = 1,
      xaxt = "n", yaxt = "n", bty = "n")
spatq:::plot_colorbar(c(min(qd_hab), max(qd_hab)),
                      hab_col)
image(qd_pref,
      main = "Preference",
      asp = 1,
      xaxt = "n", yaxt = "n", bty = "n",
      col = hcl.colors(12, "cividis"))
spatq:::plot_colorbar(c(min(qd_pref), max(qd_pref)),
                      pref_col)
image(qd_pop,
      main = "Abundance",
      asp = 1,
      xaxt = "n", yaxt = "n", bty = "n",
      col = hcl.colors(12, "viridis"))
spatq:::plot_colorbar(c(min(qd_pop), max(qd_pop)),
                      pop_col)
dev.off()

## habq example
habq_rep_fn <- function(h, p = 4) ifelse(h, p, 1)
h5 <- h5file("habq/repl_01/habq_prep.h5", mode = "r")
habq_hab <- h5[["rocky_hab"]]$read()
