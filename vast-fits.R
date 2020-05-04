library(tidyverse)
library(VAST)
library(FishStatsUtils)
## devtools::load_all("~/src/FishStatsUtils", helpers = FALSE, export_all = FALSE)
## devtools::load_all("~/src/VAST", helpers = FALSE, export_all = FALSE)
## devtools::load_all("~/src/spatq", helpers = FALSE)

## coordref <- read_csv("coordref.csv",
##                      col_types = cols(
##                        loc_idx = col_integer(),
##                        s1 = col_double(),
##                        s2 = col_double()))
## catch_df <- read_csv("repl_01/catch_01_pref.csv",
##                      col_types = cols(
##                        time = col_integer(),
##                        vessel_idx = col_integer(),
##                        loc_idx = col_integer(),
##                        coordinates = col_character(),
##                        effort = col_double(),
##                        catch_biomass = col_double())) %>%
##   left_join(coordref, by = "loc_idx") %>%
##   select(-coordinates, -loc_idx)
## write_csv(catch_df, "catch_sim.csv")

catch_df <- read.csv("catch_sim.csv") %>%
  mutate(s1 = (s1 - 50) / 100,
         s2 = (s2 - 50) / 100)

extrap_grid <- cross_df(list(Lon = (seq(0.5, 99.5, 1) - 50) / 100,
                             Lat = (seq(0.5, 99.5, 1) - 50) / 100)) %>%
  mutate(Area_km2 = 1) %>%
  as.matrix

vast_extrap <- make_extrapolation_info(Region = "user",
                                       zone = NA,
                                       strata.limits =
                                         data.frame(STRATA = "All_areas"),
                                       input_grid = extrap_grid)

vast_spat <- make_spatial_info(n_x = 404,
                               Lon_i = catch_df$s1,
                               Lat_i = catch_df$s2,
                               Extrapolation_List = vast_extrap,
                               Method = "Mesh",
                               fine_scale = TRUE,
                               Save_Results = FALSE,
                               refine = TRUE)

vast_data <- make_data(b_i = catch_df$catch_biomass,
                       a_i = catch_df$effort,
                       ## c_iz needs to be zero-indexed, as
                       ## VAST/R/make_data.R:155 assigns `n_c = max(c_iz) + 1`
                       c_iz = matrix(rep_along(catch_df$catch_biomass, 0)),
                       t_iz = catch_df$time,
                       FieldConfig = c(Omega1 = 1, Epsilon1 = 0,
                                       Omega2 = 1, Epsilon2 = 0),
                       ObsModel_ez = c(PosDist = 1, Link = 1),
                       spatial_list = vast_spat,
                       Aniso = FALSE,
                       CheckForErrors = TRUE)

vast_model <- make_model(TmbData = vast_data,
                         Version = get_latest_version())

vast_fit <- optim(vast_model$Obj$par, vast_model$Obj$fn, vast_model$Obj$gr,
                  method = "L-BFGS-B")
vast_sdr <- sdreport(vast_model$Obj)

time_string <- format(Sys.time(), "%Y-%m-%dT%H:%M")
saveRDS(vast_sdr, paste0("vast_sdr_", time_string, ".RData"))
