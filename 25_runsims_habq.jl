using Distributed

n_sims = parse(Int, ARGS[1])

pmap(1:n_sims) do rlz
    run_sims(HabQSpec, rlz)
end
