#=
Bipartite density of an induced subgraph.

  include(joinpath(@__DIR__, "density.jl"))

  density(g, SubGraph(Set(left), Set(right)))

Density is |E| / (|U| |V|). |U| and |V| are the two sides of the subgraph.
|E| counts edges of g with one end in each side. A side with no vertices
scores 0. g is not mutated.
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "..", "core.jl"))

"""
    density(g, sg::SubGraph) -> Float64

`|E| / (|U| |V|)` on the subgraph `sg` induces in `g`. `|E|` counts edges
of `g` between `sg.U` and `sg.V`. Either side empty scores `0.0`.
"""
function density(g::Union{BipartiteGraph,FrozenBipartite}, sg::SubGraph)
    n_u = length(sg.U)
    n_v = length(sg.V)
    (n_u == 0 || n_v == 0) && return 0.0
    return Subgraph.edge_count(g, sg) / (Float64(n_u) * Float64(n_v))
end
