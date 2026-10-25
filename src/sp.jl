#=
Shortest-path distance between right-side vertices of a bipartite graph.

  include(joinpath(@__DIR__, "sp.jl"))

  shortest_path(g, r1, r2, 2)       # 2 when they share a neighbor, else 3
  d = shortest_path_metric(g, 2)    # same rule, one search per source
  h = map_graph(g, 2, d)

The length is the number of edges on a shortest path in `g`. Right-side
vertices are an even number of hops apart, so two items that share a buyer
are at distance 2. The walk stops before taking a hop past `epsilon`. A pair
not reached by then — or a pair with a missing vertex — is reported as
`epsilon + 1`, which is not the true length.
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

function _sp_require_epsilon(epsilon::Real)
    if !(epsilon >= 0 && isfinite(epsilon))
        throw(ArgumentError("epsilon must be finite and ≥ 0"))
    end
    return epsilon
end

# Ids are not unique across the bipartition, so callers have to say which side.
_has_right(g::BipartiteGraph, v::Int) = haskey(g.adjV, v)
_has_right(g::FrozenBipartite, v::Int) = haskey(g.v_index, v)

# Left-side neighbors of a right vertex, and right-side neighbors of a left vertex.
@inline _left_neighbors(g::BipartiteGraph, v::Int) = g.adjV[v]
@inline _right_neighbors(g::BipartiteGraph, u::Int) = g.adjU[u]
@inline _left_neighbors(g::FrozenBipartite, v::Int) = neighbors_v(g, v)
@inline _right_neighbors(g::FrozenBipartite, u::Int) = neighbors_u(g, u)

# Largest hop count that can still be ≤ epsilon. Capped so `limit - d` cannot overflow.
function _hop_limit(epsilon::Real)
    cap = typemax(Int) - 2
    return epsilon >= cap ? cap : floor(Int, epsilon)
end

"""
Distances from right-side `source` to other right-side vertices, up to the
largest even hop count that is still ≤ `epsilon`.

Each iteration steps source → left neighbor → right neighbor, adding 2.
`stop_at` returns as soon as that vertex is reached; the dict is then only
complete for vertices discovered before it.
"""
function _bfs_ball(g, source::Int, epsilon::Real; stop_at::Union{Int,Nothing}=nothing)
    limit = _hop_limit(epsilon)
    dist = Dict{Int,Int}()
    # Separate sets: the same Int can name both a left vertex and a right vertex.
    seen_right = Set{Int}((source,))
    seen_left = Set{Int}()
    frontier = [source]
    d = 0
    while !isempty(frontier) && limit - d >= 2
        lefts = Int[]
        for r in frontier
            for u in _left_neighbors(g, r)
                u in seen_left && continue
                push!(seen_left, u)
                push!(lefts, u)
            end
        end
        d += 2
        next = Int[]
        for u in lefts
            for r in _right_neighbors(g, u)
                r in seen_right && continue
                push!(seen_right, r)
                dist[r] = d
                if stop_at !== nothing && r == stop_at
                    return dist
                end
                push!(next, r)
            end
        end
        frontier = next
    end
    return dist
end

"""
    shortest_path(g, r1, r2, epsilon) -> Real

Hop length of a shortest path between right-side vertices `r1` and `r2`.
The search does not continue past `epsilon`. If every path is longer, or
either id is not a right-side vertex, the result is `epsilon + 1`.
`r1 == r2` is `0`.
"""
function shortest_path(g::Union{BipartiteGraph,FrozenBipartite}, r1::Int, r2::Int, epsilon::Real)
    _sp_require_epsilon(epsilon)
    r1 == r2 && return 0
    if !_has_right(g, r1) || !_has_right(g, r2)
        return epsilon + 1
    end
    dist = _bfs_ball(g, r1, epsilon; stop_at=r2)
    return get(dist, r2, epsilon + 1)
end

"""
    shortest_path_metric(g, epsilon) -> Function

`(r1, r2) ->` hop length capped at `epsilon + 1`, suitable for `map_graph`.
The ball around a source is searched once and reused for later targets.
Use the same `epsilon` `map_graph` thresholds on, so a pair the search
stopped on is a pair the map leaves unconnected.
"""
function shortest_path_metric(g::Union{BipartiteGraph,FrozenBipartite}, epsilon::Real)
    _sp_require_epsilon(epsilon)
    sentinel = epsilon + 1
    cache = Dict{Int,Dict{Int,Int}}()
    function metric(r1::Int, r2::Int)
        r1 == r2 && return 0
        if !_has_right(g, r1) || !_has_right(g, r2)
            return sentinel
        end
        if !haskey(cache, r1)
            if haskey(cache, r2)
                return get(cache[r2], r1, sentinel)
            end
            cache[r1] = _bfs_ball(g, r1, epsilon)
        end
        return get(cache[r1], r2, sentinel)
    end
    return metric
end
