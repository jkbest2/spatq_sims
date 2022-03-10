using SpatQSims

simtypes = [
    QDevScalingSpec,
    # SharedQSpec,
    PrefIntensitySpec,
    DensityDependentQSpec,
    HabQSpec,
    BycatchSpec
]

repls = 1:100

for simtype in simtypes
    run_sims(simtype,
             repls;
             checkpoint = true,
             csv = false,
             feather = true)
end
