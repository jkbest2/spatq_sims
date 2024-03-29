import Pkg
Pkg.activate(".")

using FisherySim, SpatQSims
using ProgressMeter
using DataFrames, CategoricalArrays
using Statistics
using AlgebraOfGraphics, GLMakie
using Arrow
using GLM

include("10_prep.jl")

n_repl = 50

survey_stations = survey_targeting(PrefIntensitySpec(1, 1)).locations

survey_vec = Vector{DataFrame}()
pop_vec = Vector{DataFrame}()
@showprogress for rlz in 1:n_repl
    prep = prep_slowmove(rlz, distrs)
    spec = PrefIntensitySpec(rlz, 4)
    setup = SpatQSimSetup(spec, prep)
    result = simulate(setup)
    survey_pop = [(repl = rlz,
                   year = y,
                   loc_idx = loc,
                   pop = result.pop_state[y][loc])
                  for y in 1:20,
                      loc in survey_stations] |> vec |> DataFrame
    survey_catch = filter(:vessel_idx => ==(1),
                          DataFrame(result.catch_record))
    rename!(survey_catch, :time => :year)
    survey_df = leftjoin(survey_catch, survey_pop, on = [:year, :loc_idx])

    tot_pop = [(repl = rlz,
                year = y,
                tot_pop = sum(result.pop_state[y])) for y in 1:20] |> vec |> DataFrame
    pop_df = leftjoin(combine(groupby(survey_df, :year), :catch_biomass => sum => :tot_catch, :pop => sum => :survey_pop),
                      tot_pop, on = :year)

    push!(survey_vec, survey_df)
    push!(pop_vec, pop_df)
end

rescale(pop::Vector) = pop ./ exp(mean(log.(pop)))

survey_df = reduce(vcat, survey_vec)
disallowmissing!(survey_df)
Arrow.write("survey_df.feather", survey_df)

for rlz in 1:n_repl
    transform!(pop_vec[rlz],
               :repl => categorical => :repl,
               :tot_pop => rescale => :index_pop,
               :tot_catch => rescale => :index_catch,
               :survey_pop => rescale => :index_survey)
end
pop_df = reduce(vcat, pop_vec)
disallowmissing!(pop_df)
Arrow.write("pop_df.feather", pop_df)

data(pop_df) *
    mapping(:index_pop, :index_survey, color = :repl) +
    data(DataFrame(index_pop = [0, 3],
                   index_survey = [0, 3])) *
    mapping(:index_pop, :index_survey) *
    visual(Lines) |> draw

mod_st = lm(@formula(log(survey_pop) ~ -1 + log(tot_pop) + repl), pop_df)
mod_ct = lm(@formula(log(tot_catch) ~ -1 + log(tot_pop) + repl), pop_df)

function bias_metric(formula, df::AbstractDataFrame)
    mod = lm(formula, df)
    coef(mod)[1]
end

form_st = @formula(log(pop) ~ -1 + log(tot_pop) + repl)
form_ct = @formula(log(catch_biomass) ~ -1 + log(tot_pop) + repl)

survey2_df = leftjoin(survey_df, pop_df[!, [:repl, :year, :tot_pop]], on = [:repl, :year])
survey_gdf = groupby(survey2_df, [:repl, :loc_idx])

delta_st = Vector{Float64}()
# delta_ct = Vector{Float64}()
for sdf in survey_gdf
    push!(delta_st, bias_metric(form_st, sdf))
    # push!(delta_ct, bias_metric(form_ct, sdf))
end

survey3_df = combine(survey_gdf,
                     :repl => unique => :repl,
                     :loc_idx => unique => :loc_idx,
                     :pop => maximum => :pop_max)
survey3_df.delta = delta_st

## Shows some systematic bias/index response depending on initial population
## under pref. sampling. Why does this only show up in spatial estimation
## models?
data(survey3_df) *
    mapping(:pop_max, :delta) |> draw
