@everywhere using Distributed
@everywhere using SpatQSims

simtypes = [
    QDevScalingSpec,
    # SharedQSpec,
    PrefIntensitySpec,
    DensityDependentQSpec,
    HabQSpec,
    BycatchSpec
]

repls = 1:100

specs  = Any[]
for st in simtypes
    svals = sim_values(st)
    for r in repls, sv in svals
        push!(specs, st(r, sv))
    end
end

pmap(specs) do sp
    run_sims([sp];
             checkpoint = true,
             csv = false,
             feather = true)
end
