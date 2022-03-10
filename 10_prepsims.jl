using SpatQSims

simtypes = [
    QDevScalingSpec,
    # SharedQSpec,
    PrefIntensitySpec,
    DensityDependentQSpec,
    HabQSpec,
    BycatchSpec
]

for simtype in simtypes
    prep_sims(simtype, 100)
end
