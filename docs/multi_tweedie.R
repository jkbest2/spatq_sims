library(tweedie)

## Remove fish in CPG-distributed increments for n years
nfish <- function(n, p0, q, xi = 1.84, phi = 1.2) {
  p <- rep(0, n)
  p[1] <- p0
  for (yr in 2:n) {
    ctch <- min(rtweedie(1, xi, p[yr - 1] * q, phi), p[yr - 1])
    p[yr] <- p[yr - 1] - ctch
    if (p[yr] == 0) {
      break
    }
  }
  p
}

f1 <- replicate(10000, nfish(50, 0.01, 0.01), simplify = FALSE)
hist(map_dbl(f1, ~ sum(. > 0)), breaks = 0:50)
