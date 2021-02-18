using Distributed

n_sims = parse(Int, ARGS[1])

# pmap(1:n_sims) do repl
#     run_sharedq_sim(repl, "prep.h5")
# end

for repl in 1:n_sims
    run_sharedq_sim(repl, "prep.h5")
end
