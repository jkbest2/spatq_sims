pmap(1:100) do repl
    Pdict, Cdict = run_scenarios(repl, "prep.h5")
    save_scenarios(repl, Pdict, Cdict)
end
