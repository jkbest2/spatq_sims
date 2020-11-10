using Distributed

n_sims = parse(Int, ARGS[1])
pmap(1:n_sims) do repl
    Pdict, Cdict = run_scenarios(repl, "prep.h5")
    save_scenarios(repl, Pdict, Cdict)
end
