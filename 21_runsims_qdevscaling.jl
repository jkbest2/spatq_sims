using Distributed

n_sims = parse(Int, ARGS[1])

run_sims(QDevScalingSpec, n_sims)

# pmap(1:n_sims) do repl
#     run_qdevscaling_sim(repl, "prep.h5")
# end

# for repl in 1:n_sims
#     run_qdevscaling_sim(repl, "prep.h5")
# end
