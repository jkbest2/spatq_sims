using Plots
using TweedieDistributions
using HDF5
using SpatQSims
using Statistics

spec = QDevScalingSpec(1, 0.05, "prep.h5")
pop = load_spatial_eq("prep.h5", "qdevscaling", 1)

heatmap(pop.P, aspect_ratio = 1, cmap = :cividis)

qtiles = [0.05, 0.1, 0.5, 0.9, 0.95]
pop_qtile = quantile.(Ref(vec(pop.P)), qtiles)

# Tweedie parameters
ϕ = tweedie_dispersion(spec)
p = tweedie_shape(spec)

# Base catchability
q = 0.01
