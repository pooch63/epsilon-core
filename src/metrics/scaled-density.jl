#=
Scaled bipartite density of an induced subgraph.

  include(joinpath(@__DIR__, "scaled-density.jl"))

  scaled_density(g, SubGraph(Set(left), Set(right)))

Scaled density is |E| / sqrt(|U| |V|). The denominator is the geometric
mean of the two side sizes. A side with no vertices scores 0. g is not
mutated.
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "..", "core.jl"))

"""
    scaled_density(g, sg::SubGraph) -> Float64

`|E| / sqrt(|U| |V|)` on the subgraph `sg` induces in `g`. `|E|` counts
edges of `g` between `sg.U` and `sg.V`. Either side empty scores `0.0`.
"""
function scaled_density(g::Union{BipartiteGraph,FrozenBipartite}, sg::SubGraph)
    n_u = length(sg.U)
    n_v = length(sg.V)
    (n_u == 0 || n_v == 0) && return 0.0
    return Subgraph.edge_count(g, sg) / sqrt(Float64(n_u) * Float64(n_v))
end
