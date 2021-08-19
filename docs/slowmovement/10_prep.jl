# import Pkg
# Pkg.activate(".")

# using FisherySim, SpatQSims
using Distributions
# using HDF5
# import Random

struct PrepDistributions{D, H, Q}
    domain::D
    habitat_distribution::H
    log_catchability_devs_distr::Q

    function PrepDistributions(domain::D,
                               habitat_distribution::H,
                               log_catchability_devs_distr::Q) where {D<:AbstractFisheryDomain,
                                                                     H<:DomainDistribution,
                                                                     Q<:DomainDistribution}
        new{D, H, Q}(domain,
                     habitat_distribution,
                     log_catchability_devs_distr)
    end
end

domain = GriddedFisheryDomain()

#-Habitat-------------------------------------------------------------------
# Declare habitat distribution
habitat_kernel = Matérn32Cov(9.0, 20.0)
function hab_mean_fn(loc)
    hab_mean = ((-loc[1] + 50)^2 + (-loc[2] + 50)^2) / 500 - 5//2
    if hab_mean > 5
        hab_mean = 5one(hab_mean)
    end
    hab_mean
end
habitat_mean = hab_mean_fn.(domain.locs)
habitat_cov = cov(habitat_kernel, domain)
habitat_distribution = DomainDistribution(MvNormal(vec(habitat_mean),
                                                   habitat_cov), domain)

#-Catchability--------------------------------------------------------------
# Log-catchability deviations - rescalable option so that spatial structure
# of catchability deviations is constant, but has different magnitudes.
log_catchability_devs_kern = Matérn32Cov(1.0, 30.0)
log_catchability_devs_cov = cov(log_catchability_devs_kern, domain)
log_catchability_devs_distr = DomainDistribution(
    MvNormal(zeros(length(domain)),
             log_catchability_devs_cov),
    domain)

distrs = PrepDistributions(domain, habitat_distribution, log_catchability_devs_distr)

function prep_slowmove(rlz::Integer, distrs::PrepDistributions, K::Real = 100.0)
    habitat = rand(distrs.habitat_distribution)

    #-Movement------------------------------------------------------------------
    # Declare habitat preference function and plot realized preference
    hab_pref_fn(h) = exp(-(h + 5)^2 / 25)
    # Distance travelled in a single step decays exponentially
    distance_fn(d) = exp(-d ^ 2 / 5)
    # distance_fn(d) = 1
    movement = MovementModel(distrs.domain, habitat, hab_pref_fn, distance_fn)

    init_pop = eqdist(movement, K)

    # Generate *log*-catchability deviations that can be rescaled
    # depending on the exact OM specification
    log_qdev = rand(distrs.log_catchability_devs_distr)

    SpatQSimPrep(rlz, movement, init_pop, log_qdev)
end
