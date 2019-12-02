using SpatQSims
using CSV

pmap(1:100) do repl
  for scen in ["naive", "simple", "scaled", "shared"]
      dirnm = "repl_" * string(repl, pad = 2)
      ctchflnm =  "catch_" * string(repl, pad = 2) * scen * ".csv"
      popflnm = "pop_" * string(repl, pad = 2) * ".h5"

      P = Array{Float64, 3}(undef, 100, 100, 25)
      h5read(joinpath(dirnm, popflnm)) do fid
          P .= fid[sc]["popstate"]


