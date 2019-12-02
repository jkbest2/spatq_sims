pmap(1:100) do repl
    Pdict, Cdict = run_scenarios(repl, "prep.h5")
    save_scenarios(repl, Pdict, Cdict)
    for scen in keys(Pdict)
        gif_file = "repl_" * string(repl, pad = 2) * "/" *
                    string(scen) * ".gif"
        plot_gif(Pdict[scen], Cdict[scen], gif_file)
    end
end
