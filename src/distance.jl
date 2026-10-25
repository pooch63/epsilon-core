#=
Neighborhood distances between right-side vertices.

  include(joinpath(@__DIR__, "distance.jl"))

  d = overlap_metric(g)
  d(r1, r2)                          # 0 when one neighborhood contains the other
  h = map_graph(g, 0.5, d)

`N(r)` is the left-side neighborhood of a right-side vertex. Each metric is
that intersection count turned into a distance in [0, 1]:

  overlap   1 - |N(u) ∩ N(v)| / min(|N(u)|, |N(v)|)
  jaccard   1 - |N(u) ∩ N(v)| / |N(u) ∪ N(v)|
  dice      1 - 2|N(u) ∩ N(v)| / (|N(u)| + |N(v)|)
  cosine    1 - |N(u) ∩ N(v)| / sqrt(|N(u)| · |N(v)|)

`u == v` is 0. A vertex that is not on the right side, or a denominator of
0, is distance 1. The `*_metric` constructors cache each neighborhood the
first time it is read; pass those to `map_graph`. The `*_distance` functions
answer one pair and do not cache.

Diffusion distance is not a neighborhood coefficient. `diffusion_metric`
embeds every right vertex from the spectrum of a random walk and returns
Euclidean distance in that embedding. It is not scaled into [0, 1].
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

using LinearAlgebra

const _DistanceGraph = Union{BipartiteGraph,FrozenBipartite}

@inline _right_present(g::BipartiteGraph, v::Int) = haskey(g.adjV, v)
@inline _right_present(g::FrozenBipartite, v::Int) = haskey(g.v_index, v)

@inline _right_neighborhood(g::BipartiteGraph, v::Int) = g.adjV[v]
@inline _right_neighborhood(g::FrozenBipartite, v::Int) = neighbors_v(g, v)

function _neighborhood_set(g::_DistanceGraph, r::Int)
    _right_present(g, r) || return nothing
    return Set{Int}(_right_neighborhood(g, r))
end

function _count_intersect(a::Set{Int}, b::Set{Int})
    if length(a) > length(b)
        a, b = b, a
    end
    n = 0
    for x in a
        x in b && (n += 1)
    end
    return n
end

function _overlap(a::Set{Int}, b::Set{Int})
    denom = min(length(a), length(b))
    denom == 0 && return 1.0
    return 1 - _count_intersect(a, b) / denom
end

function _jaccard(a::Set{Int}, b::Set{Int})
    inter = _count_intersect(a, b)
    denom = length(a) + length(b) - inter
    denom == 0 && return 1.0
    return 1 - inter / denom
end

function _dice(a::Set{Int}, b::Set{Int})
    denom = length(a) + length(b)
    denom == 0 && return 1.0
    return 1 - 2 * _count_intersect(a, b) / denom
end

function _cosine(a::Set{Int}, b::Set{Int})
    n1, n2 = length(a), length(b)
    (n1 == 0 || n2 == 0) && return 1.0
    return 1 - _count_intersect(a, b) / sqrt(n1 * n2)
end

function _apply_distance(g::_DistanceGraph, r1::Int, r2::Int, distance)
    r1 == r2 && return 0.0
    a = _neighborhood_set(g, r1)
    b = _neighborhood_set(g, r2)
    (a === nothing || b === nothing) && return 1.0
    return distance(a, b)
end

function _neighborhood_metric(g::_DistanceGraph, distance)
    cache = Dict{Int,Union{Nothing,Set{Int}}}()
    function metric(r1::Int, r2::Int)
        r1 == r2 && return 0.0
        a = get!(cache, r1) do
            _neighborhood_set(g, r1)
        end
        b = get!(cache, r2) do
            _neighborhood_set(g, r2)
        end
        (a === nothing || b === nothing) && return 1.0
        return distance(a, b)
    end
    return metric
end

"""
    overlap_distance(g, r1, r2) -> Float64

Overlap-coefficient distance, `1 - |N(r1) ∩ N(r2)| / min(|N(r1)|, |N(r2)|)`.
This is 0 when one neighborhood contains the other. `r1 == r2` is 0. A
missing right-side vertex, or a neighborhood of size 0, yields 1.
"""
overlap_distance(g::_DistanceGraph, r1::Int, r2::Int) = _apply_distance(g, r1, r2, _overlap)

"""
    overlap_metric(g) -> Function

`(r1, r2) ->` [`overlap_distance`](@ref) for `map_graph`. Neighborhoods are
cached on first read.
"""
overlap_metric(g::_DistanceGraph) = _neighborhood_metric(g, _overlap)

"""
    jaccard_distance(g, r1, r2) -> Float64

Jaccard distance, `1 - |N(r1) ∩ N(r2)| / |N(r1) ∪ N(r2)|`. `r1 == r2` is 0.
A missing right-side vertex, or two empty neighborhoods, yields 1.
"""
jaccard_distance(g::_DistanceGraph, r1::Int, r2::Int) = _apply_distance(g, r1, r2, _jaccard)

"""
    jaccard_metric(g) -> Function

`(r1, r2) ->` [`jaccard_distance`](@ref) for `map_graph`. Neighborhoods are
cached on first read.
"""
jaccard_metric(g::_DistanceGraph) = _neighborhood_metric(g, _jaccard)

"""
    dice_distance(g, r1, r2) -> Float64

Sørensen–Dice distance, `1 - 2|N(r1) ∩ N(r2)| / (|N(r1)| + |N(r2)|)`.
`r1 == r2` is 0. A missing right-side vertex, or two empty neighborhoods,
yields 1.
"""
dice_distance(g::_DistanceGraph, r1::Int, r2::Int) = _apply_distance(g, r1, r2, _dice)

"""
    dice_metric(g) -> Function

`(r1, r2) ->` [`dice_distance`](@ref) for `map_graph`. Neighborhoods are
cached on first read.
"""
dice_metric(g::_DistanceGraph) = _neighborhood_metric(g, _dice)

"""
    cosine_distance(g, r1, r2) -> Float64

Cosine (Ochiai) distance, `1 - |N(r1) ∩ N(r2)| / sqrt(|N(r1)| · |N(r2)|)`.
`r1 == r2` is 0. A missing right-side vertex, or an empty neighborhood,
yields 1.
"""
cosine_distance(g::_DistanceGraph, r1::Int, r2::Int) = _apply_distance(g, r1, r2, _cosine)

"""
    cosine_metric(g) -> Function

`(r1, r2) ->` [`cosine_distance`](@ref) for `map_graph`. Neighborhoods are
cached on first read.
"""
cosine_metric(g::_DistanceGraph) = _neighborhood_metric(g, _cosine)

# --- lazy and two-step diffusion distance on the right side ----------------
#
# Lazy walk: P_L = (I + D^{-1} A) / 2, eigenvalues μ = (1+λ)/2.
# Two-step walk on V: P_V = D_V^{-1} Bᵀ D_U^{-1} B, eigenvalues σ(M)^2.
# Both come from the SVD of M = D_U^{-1/2} B D_V^{-1/2}. A singular value σ
# with right vector b gives
#   lazy ±σ:  μ = (1±σ)/2,  φ_V = D_V^{-1/2} (±b) / √2
#   lazy kernel on V (Mb = 0): μ = 1/2, φ_V = D_V^{-1/2} b
#   two-step: λ = σ^2, φ_V = D_V^{-1/2} b
# D_t(i, j) = √vol ‖μ^t φ(i) − μ^t φ(j)‖, which equals
# Σ_k (P^t[i,k] − P^t[j,k])^2 / π_k. vol is the volume of the chain
# (both sides for the lazy walk, the right side for the two-step walk).

_as_weight(::Nothing) = 1.0
_as_weight(w::Real) = Float64(w)
_as_weight(w) = throw(ArgumentError("diffusion distance needs numeric edge weights, got $(typeof(w))"))

function _require_weight(w::Float64)
    if !(w > 0 && isfinite(w))
        throw(ArgumentError("diffusion distance needs positive finite edge weights"))
    end
    return w
end

function _each_edge(f, g::BipartiteGraph{Nothing})
    for (u, nbrs) in g.adjU
        for v in nbrs
            f(u, v, 1.0)
        end
    end
    return nothing
end

function _each_edge(f, g::BipartiteGraph)
    for (u, nbrs) in g.adjU
        for v in nbrs
            w = _require_weight(_as_weight(g.edge_data[(u, v)]))
            f(u, v, w)
        end
    end
    return nothing
end

function _each_edge(f, g::FrozenBipartite)
    for ui in eachindex(g.u_ids)
        u = g.u_ids[ui]
        for k in neighbor_range_u(g, ui)
            w = _require_weight(_as_weight(g.edge_data[k]))
            f(u, g.v_ids[g.v_adj[k]], w)
        end
    end
    return nothing
end

"""Biadjacency of the edges, with isolated vertices left out. `nothing` if there are none."""
function _biadjacency(g::_DistanceGraph)
    us = Int[]
    vs = Int[]
    ws = Float64[]
    _each_edge(g) do u, v, w
        push!(us, u)
        push!(vs, v)
        push!(ws, w)
    end
    isempty(ws) && return nothing
    u_ids = sort!(unique(us))
    v_ids = sort!(unique(vs))
    u_of = Dict(id => i for (i, id) in enumerate(u_ids))
    v_of = Dict(id => i for (i, id) in enumerate(v_ids))
    B = zeros(Float64, length(u_ids), length(v_ids))
    for i in eachindex(ws)
        B[u_of[us[i]], v_of[vs[i]]] = ws[i]
    end
    return (; B, v_ids, dv=vec(sum(B; dims=1)))
end

_normalize_walk(walk::Symbol) =
    walk in (:lazy, :step) ? walk :
    throw(ArgumentError("diffusion walk must be lazy or step"))
_normalize_walk(walk::AbstractString) = _normalize_walk(Symbol(walk))

function _diffusion_settings(time::Integer, rank::Union{Integer,Nothing}, walk)
    if time < 0
        throw(ArgumentError("diffusion time must be ≥ 0"))
    end
    if rank !== nothing && rank < 1
        throw(ArgumentError("diffusion rank must be ≥ 1"))
    end
    return Int(time), rank === nothing ? nothing : Int(rank), _normalize_walk(walk)
end

"""True when `phi` is not constant on the right side. Constant modes cancel in distances."""
function _varies_on_right(phi::Vector{Float64})
    lo, hi = extrema(phi)
    span = hi - lo
    span > 0 || return false
    return span > 1e-9 * max(abs(lo), abs(hi), 1e-12)
end

function _diffusion_modes(B::Matrix{Float64}, dv::Vector{Float64}, walk::Symbol)
    du = vec(sum(B; dims=2))
    M = Diagonal(1 ./ sqrt.(du)) * B * Diagonal(1 ./ sqrt.(dv))
    F = svd(M; full=true)
    inv_sqrt = 1 ./ sqrt.(dv)
    modes = Tuple{Float64,Vector{Float64}}[]
    nσ = length(F.S)
    if walk === :lazy
        for k in 1:size(F.V, 2)
            phi = inv_sqrt .* F.V[:, k]
            if k <= nσ && F.S[k] > 0
                σ = F.S[k] > 1 ? 1.0 : F.S[k]
                phi = phi ./ sqrt(2)
                push!(modes, ((1 + σ) / 2, phi))
                push!(modes, ((1 - σ) / 2, -phi))
            else
                push!(modes, (0.5, phi))
            end
        end
        vol = sum(du) + sum(dv)
    else
        for k in 1:size(F.V, 2)
            σ = k <= nσ ? (F.S[k] > 1 ? 1.0 : F.S[k]) : 0.0
            push!(modes, (σ^2, inv_sqrt .* F.V[:, k]))
        end
        vol = sum(dv)
    end
    return modes, vol
end

function _select_modes(modes, time::Int, rank::Union{Int,Nothing})
    kept = Tuple{Float64,Vector{Float64}}[]
    for (μ, phi) in modes
        time > 0 && μ == 0 && continue
        _varies_on_right(phi) || continue
        push!(kept, (μ, phi))
    end
    if rank !== nothing
        sort!(kept; by=first, rev=true)
        rank < length(kept) && resize!(kept, rank)
    end
    return kept
end

function _diffusion_embedding(g::_DistanceGraph, time::Int, rank::Union{Int,Nothing}, walk::Symbol)
    built = _biadjacency(g)
    index = Dict{Int,Int}()
    built === nothing && return index, zeros(0, 0)
    modes, vol = _diffusion_modes(built.B, built.dv, walk)
    kept = _select_modes(modes, time, rank)
    nV = length(built.v_ids)
    E = zeros(nV, length(kept))
    scale = sqrt(vol)
    for (col, (μ, phi)) in enumerate(kept)
        weight = scale * μ^time
        weight == 0 && continue
        E[:, col] = weight .* phi
    end
    for (i, id) in enumerate(built.v_ids)
        index[id] = i
    end
    return index, E
end

"""
    diffusion_metric(g; time=1, rank=nothing, walk=:lazy) -> Function

`(r1, r2) ->` diffusion distance at time `time` between right-side vertices.

`walk=:lazy` (the default) is the lazy random walk `P_L = (I + D^{-1} A) / 2`.
`walk=:step` is the two-step walk on the right side,
`P_V = D_V^{-1} Bᵀ D_U^{-1} B`, so one step goes right → left → right.
`D_t(i, j)^2 = Σ_k (P^t[i,k] − P^t[j,k])^2 / π_k` with `π_k = d_k / vol`
on the vertices of that chain. The value is exact when `rank` is `nothing`.
`rank` keeps the `rank` slowest modes (largest eigenvalue). The error from
dropping the rest is controlled by the next eigenvalue raised to `time`.

`time` is an integer ≥ 0. On a connected graph the distance tends to 0 as
`time` grows. Under the lazy walk a vertex still carries its own leftover
mass, so two right vertices with the same neighbors are not at distance 0
for small `time`. Under the two-step walk they are.

`r1 == r2` is 0. A right vertex with no incident edge, or an id that is not
on the right side, is at distance `Inf`. `Nothing` edge data has weight 1.
Real edge data is the weight and must be positive and finite. The map is one
dense SVD of the `|U| × |V|` normalized biadjacency.
"""
function diffusion_metric(g::_DistanceGraph; time::Integer=1, rank::Union{Integer,Nothing}=nothing, walk=:lazy)
    time, rank, walk = _diffusion_settings(time, rank, walk)
    index, E = _diffusion_embedding(g, time, rank, walk)
    function metric(r1::Int, r2::Int)
        r1 == r2 && return 0.0
        i = get(index, r1, 0)
        j = get(index, r2, 0)
        (i == 0 || j == 0) && return Inf
        acc = 0.0
        @inbounds for col in axes(E, 2)
            delta = E[i, col] - E[j, col]
            acc += delta * delta
        end
        # Roundoff on a symmetric pair sits well below any distance the
        # walk can mean. Returning 0 keeps epsilon = 0 meaningful.
        dist = sqrt(acc)
        return dist <= 1e-12 ? 0.0 : dist
    end
    return metric
end

"""
    diffusion_distance(g, r1, r2; time=1, rank=nothing, walk=:lazy) -> Float64

[`diffusion_metric`](@ref) for one pair. The spectrum is rebuilt on every
call; use `diffusion_metric` when many pairs share a graph.
"""
function diffusion_distance(g::_DistanceGraph, r1::Int, r2::Int; time::Integer=1, rank::Union{Integer,Nothing}=nothing, walk=:lazy)
    return diffusion_metric(g; time, rank, walk)(r1, r2)
end
