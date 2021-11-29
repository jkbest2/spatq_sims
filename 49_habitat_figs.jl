using HDF5
using Plots
using SpatQSims

## qdevscaling example
qd_pref_fn(h) = exp(-(h + 5)^2 / 25)
qd_hab = h5read("prep.h5", "habitat")[:, :, 1]'
qd_pref = qd_pref_fn.(qd_hab)
qd_pop = h5read("prep.h5", "spatial_eq")[:, :, 1]'

qd_hab_plot = heatmap(qd_hab;
                      c = :thermal,
                      aspect_ratio = 1,
                      title = "Habitat",
                      showaxis = false, ticks = false)
qd_pref_plot = heatmap(qd_pref;
                       c = :cividis,
                       aspect_ratio = 1,
                       title = "Preference",
                       showaxis = false, ticks = false)
qd_pop_plot = heatmap(qd_pop;
                      c = :viridis, aspect_ratio = 1,
                      title = "Abundance",
                      showaxis = false, ticks = false)
plot(qd_hab_plot, qd_pref_plot, qd_pop_plot; size = (800, 800))
savefig("opmod_results/qdevscaling_hab.svg")

# Habq example
hq_parvals = 1.75 .^ (-2:3)
hq_pref_fn(h, p = 4) = ifelse(h, p, 1)
hq_hab = h5read("habq/repl_01/habq_prep.h5", "rocky_hab")'
hq_pref = hq_pref_fn.(hq_hab, hq_parvals[5])
hq_pop = h5read("habq/repl_01/habq_05_popstate.h5", "popstate")[:, :, 1]

hq_hab_plot = heatmap(hq_hab;
                      c = :thermal,
                      aspect_ratio = 1,
                      title = "Habitat",
                      showaxis = false, ticks = false)
hq_pref_plot = heatmap(hq_pref;
                       c = :cividis,
                       aspect_ratio = 1,
                       title = "Preference",
                       showaxis = false, ticks = false)
hq_pop_plot = heatmap(hq_pop;
                      c = :viridis, aspect_ratio = 1,
                      title = "Abundance",
                      showaxis = false, ticks = false)
plot(hq_hab_plot, hq_pref_plot, hq_pop_plot; size = (800, 800))
savefig("opmod_results/habq_hab.svg")
