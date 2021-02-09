using Distributed

n_sims = parse(Int, ARGS[1])

pmap(1:n_sims) do repl
    run_qdevscaling_sim(repl, "prep.h5")
end
