library(tidyverse)
library(tweedie)
library(spatq)

## Tweedie parameters
p <- 1.84
phi <- 1.2

## base catchability
q <- 0.01

## Population map from example scenario
pop <- read_popstate("qdevscaling", 1, 1, ".", filetype = "h5")

## Calculate density for catches for mean abundance in a given cell
pop_bar <- 100 / 10000

plot_cpg <- function(mus, p, phi) {
  x_upper <- max(qtweedie(rep_along(mus, 0.95), p, mus, phi))
  x <- seq(0, x_upper, length.out = 1025)
  l <- layout(matrix(1:(2 * length(mus)), ncol = 2, byrow = TRUE), widths = c(1, 4))
  for (i in seq_along(mus)) {
    d <- dtweedie(x, p, mus[i], phi)
    barplot(d[1],
            xlab = "Catch = 0",
            ylim = c(0, 1),
            ylab = "Probability")
    mtext(paste0(letters[i], "."), adj = 0)
    plot(x[-1], d[-1], type = "l",
         xlab = "Catch > 0",
         xaxs = "i",
         ylab = "Density",
         yaxs = "i")
    abline(v = mus[i], lty = 2)
  }
}

qtiles <- c(0.05, 0.5, 0.95)
pop_qtile <- quantile(pop[, , 1], qtiles)
mus <- q * pop_qtile
svg(file.path("om_results", "cpg_plot.svg"),
    width = 8, height = 6)
par(oma = rep(0, 4),
    mar = c(4, 4, 1, 0))
plot_cpg(mus, p, phi)
dev.off()

