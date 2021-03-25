using Distributed

n_sims = parse(Int, ARGS[1])

run_sims(PrefIntensitySpec, n_sims)
