using SpatQSims
using Plots

if !isdir("opmod_results")
    mkdir("opmod_results")
end

function plot_hab(spec::SpatQSimSpec)
    prep = load_prep(spec)
    hab = habitat(prep)
    preffn = HabitatPreference(spec)
    pref = preffn(hab)
    pop = init_pop(prep)
    mov = movement(prep)

    locs = [[25, 25],
            [25, 75],
            [50, 50],
            [75, 25],
            [75, 75]]
    mov_to = permutedims(mov.M[locidx, :])
    mov_from = mov.M[:, locidx]
    mov_net = mov_to .- mov_from
    mov_net_tot = sum(mov_net; dims = 2)
    mov_ex = extrema(mov_net)
    mov_clim = maximum(abs, mov_ex)

    x = y = 0.5:1.0:99.5

    hab_plot = heatmap(x, y, hab[1]';
                       c = :thermal,
                       aspect_ratio = 1,
                       title = "Habitat",
                       showaxis = false, ticks = false)
    pref_plot = heatmap(x, y, pref';
                        c = :cividis,
                        aspect_ratio = 1,
                        title = "Preference",
                        showaxis = false, ticks = false)
    pref_plot = scatter!(getindex.(locs, 1) .- 0.5, getindex.(locs, 2) .- 0.5,
                         label = nothing,
                         color = :black)
    pop_plot = heatmap(x, y, pop.P';
                       c = :viridis, aspect_ratio = 1,
                       title = "Abundance",
                       showaxis = false, ticks = false)

    mov_cmap = cgrad(:diverging_bwr_20_95_c54_n256; rev = true)
    mov_plot = heatmap(x, y, reshape(mov_net_tot, (100, 100))';
                       c = mov_cmap,
                       clims = (-mov_clim, mov_clim),
                       aspect_ratio = 1,
                       title = "Net movement",
                       showaxis = false,
                       ticks = false)
    mov_plot = scatter!(getindex.(locs, 1) .- 0.5, getindex.(locs, 2) .- 0.5,
                        label = nothing,
                        color = :black)

    plot(hab_plot, pref_plot,
         pop_plot, mov_plot;
         size = (1200, 900))
end

function plot_name(spec::SpatQSimSpec; fileext = "png")
    studytype = typeof(spec).name.wrapper
    study = lowercase(string(studytype))[begin:(end - 4)]
    repl = realization(spec)
    opmod = sim_value_idx(spec)
    study * "_" * string(repl, pad = 2) * "_" * string(opmod, pad = 2) * "." * fileext
end

repl = 1
specs = [QDevScalingSpec(repl, 0.1),
         PrefIntensitySpec(repl, 2),
         DensityDependentQSpec(repl, 1.0),
         HabQSpec(repl, 0.25),
         HabQSpec(repl, 4.0),
         BycatchSpec(repl, 0.4)]

for spec in specs
    plot_hab(spec)
    fn = plot_name(spec, fileext = "png")
    savefig(joinpath("opmod_results", fn))
end
