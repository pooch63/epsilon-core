#=
Maximum-edge biclique by the BSC sampler of Shaham, Yu, and Li,
SDM 2016, Algorithm 1.

  include(joinpath(@__DIR__, "max-edge-biclique.jl"))

  found = max_edge_biclique(g)
  found.left
  found.right
  found.edges

  score_biclique(g, found.left, found.right, 2, 0.5, overlap_metric(g))

Repeat: draw a seed P on one side, keep the common neighborhood on the
other side, then keep every seed-side vertex that contains that
neighborhood. The returned pair is the sample with the most edges.
`g` is not mutated. The search does not read θ, ε, or a distance.

  julia --project=. test/max-edge-biclique.jl
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

using Random

"""
    is_biclique(g, left, right) -> Bool

`left` and `right` are duplicate-free, lie on those sides of `g`, and every
listed left vertex is adjacent to every listed right vertex. Either side
being empty is rejected.
"""
function is_biclique(g::BipartiteGraph, left, right)
    L = collect(Int, left)
    R = collect(Int, right)
    (isempty(L) || isempty(R)) && return false
    length(L) == length(Set{Int}(L)) || return false
    length(R) == length(Set{Int}(R)) || return false
    for v in R
        haskey(g.adjV, v) || return false
    end
    for u in L
        haskey(g.adjU, u) || return false
        nbrs = g.adjU[u]
        for v in R
            v in nbrs || return false
        end
    end
    return true
end

"""Vertices of `seeds` all touch. The smallest neighborhood is the seed of the intersection."""
function _meb_common(seeds::Vector{Int}, adj::Dict{Int,Set{Int}})
    isempty(seeds) && return Int[]
    smallest = seeds[1]
    smallest_n = length(adj[smallest])
    for s in seeds
        n = length(adj[s])
        if n < smallest_n
            smallest = s
            smallest_n = n
        end
    end
    acc = copy(adj[smallest])
    for s in seeds
        s == smallest && continue
        intersect!(acc, adj[s])
        isempty(acc) && return Int[]
    end
    return sort!(collect(acc))
end

"""Uniform sample of `p` ids. `pool` stays in its original order."""
function _meb_sample(rng::AbstractRNG, pool::Vector{Int}, p::Int)
    n = length(pool)
    p >= n && return copy(pool)
    drawn = copy(pool)
    for i in 1:p
        j = rand(rng, i:n)
        drawn[i], drawn[j] = drawn[j], drawn[i]
    end
    return drawn[1:p]
end

function _meb_keep(best, left::Vector{Int}, right::Vector{Int})
    (isempty(left) || isempty(right)) && return best
    edges = Int128(length(left)) * length(right)
    key = (edges, length(left) + length(right), length(left))
    if best === nothing || key > best.key
        return (; key, left, right, edges)
    end
    return best
end

function _meb_density(adjU::Dict{Int,Set{Int}}, m::Int, n::Int)
    edges = 0
    for nbrs in values(adjU)
        edges += length(nbrs)
    end
    return edges / (float(m) * float(n))
end

"""
Shaham, Yu, and Li fit `|P| ≈ (0.6 − 0.8β) e^{(4+3β)ρ}` on planted bicliques.
The product is clamped to at least 1. β > 0.75 makes the linear factor
negative; that still clamps to 1.
"""
function _meb_formula(rho::Float64, beta::Float64)
    p = (0.6 - 0.8 * beta) * exp((4 + 3 * beta) * rho)
    return max(1, round(Int, p))
end

function _meb_require_beta(beta::Real)
    if !(0 <= beta <= 1 && isfinite(beta))
        throw(ArgumentError("beta must be finite and in [0, 1]"))
    end
    return Float64(beta)
end

"""
    max_edge_biclique(g; beta=0.5, samples=nothing, seed_size=nothing,
                      max_samples=2_000, both_sides=true, seed=1) -> NamedTuple

A large edge biclique `(left, right)` of `g`. `edges` is `|left| * |right|`.

Each sample draws `seed_size` vertices on one side, keeps their common
neighborhood, and then keeps every vertex on the seed side whose
neighborhood contains that set. The sample with the most edges is returned.
Ties keep the pair with more vertices, then the pair with more left
vertices. `g` is not mutated. An empty side, or a run that never builds a
pair, returns empty vectors and `edges == 0`.

`both_sides` alternates the seed side, starting with the left. Shaham, Yu,
and Li note that the construction is symmetric and that alternating the
seed side helps when the biclique is unbalanced. `both_sides=false` seeds
the left side only, which is Algorithm 1 as printed.

`seed_size` defaults to the planted-biclique fit above, with `beta` the
guess for the fraction of the opposite side the biclique covers. `samples`
defaults to 1% of `|L| * |R|`, the paper's easy setting, and is capped by
`max_samples` so a large graph does not draw hundreds of millions of seeds.
Pass `samples` to use an exact budget; that value is not capped.

The optimality argument in the paper is for the ranked objective
`|I| · (1/γ)^|J|`, and it needs the seeded side of an optimal biclique to
have at least `α|L| + seed_size` vertices. The edge count `|I| * |J|` is
what this function maximizes, and it is what the paper compares with iMBEA.
A seed of size 1 misses a biclique in which every member has a private
neighbor. Raise `seed_size`, or seed both sides, when that happens. The
paper's formulas are a starting point to tune.

The search does not read θ, ε, or a distance. [`score_biclique`](@ref)
applies those after the fact.
"""
function max_edge_biclique(g::BipartiteGraph;
                           beta::Real=0.5,
                           samples::Union{Nothing,Integer}=nothing,
                           seed_size::Union{Nothing,Integer}=nothing,
                           max_samples::Integer=2_000,
                           both_sides::Bool=true,
                           seed::Integer=1)
    max_samples < 1 && throw(ArgumentError("max_samples must be ≥ 1"))
    beta64 = _meb_require_beta(beta)
    if samples !== nothing && samples < 1
        throw(ArgumentError("samples must be ≥ 1"))
    end
    if seed_size !== nothing && seed_size < 1
        throw(ArgumentError("seed_size must be ≥ 1"))
    end
    m = length(g.adjU)
    n = length(g.adjV)
    if m == 0 || n == 0
        return (; left=Int[], right=Int[], edges=0)
    end
    if samples === nothing
        target = max(1, round(Int, 0.01 * float(m) * float(n)))
        n_samples = min(target, Int(max_samples))
    else
        n_samples = Int(samples)
    end
    p = seed_size === nothing ? _meb_formula(_meb_density(g.adjU, m, n), beta64) : Int(seed_size)
    u_ids = sort!(collect(keys(g.adjU)))
    v_ids = sort!(collect(keys(g.adjV)))
    rng = Xoshiro(Int(seed))
    best = nothing
    for t in 1:n_samples
        if !both_sides || isodd(t)
            picked = _meb_sample(rng, u_ids, min(p, m))
            right = _meb_common(picked, g.adjU)
            left = _meb_common(right, g.adjV)
        else
            picked = _meb_sample(rng, v_ids, min(p, n))
            left = _meb_common(picked, g.adjV)
            right = _meb_common(left, g.adjU)
        end
        best = _meb_keep(best, left, right)
    end
    best === nothing && return (; left=Int[], right=Int[], edges=0)
    edge_count = best.edges
    if edge_count > typemax(Int)
        throw(OverflowError("biclique has too many edges to store in Int"))
    end
    return (; left=best.left, right=best.right, edges=Int(edge_count))
end

function _bcs_epsilon(epsilon::Real)
    if !(epsilon >= 0 && isfinite(epsilon))
        throw(ArgumentError("epsilon must be finite and ≥ 0"))
    end
    return epsilon
end

@inline function _bcs_dist(metric, a::Int, b::Int)
    return a <= b ? metric(a, b) : metric(b, a)
end

"""
`|ids| ≤ 1` has diameter 0, which meets `M < ε` only when `ε > 0`.
A larger set meets it when every pair is strictly inside `ε`.
"""
function _bcs_close(ids::Vector{Int}, epsilon::Real, metric)
    n = length(ids)
    n == 0 && return false
    n == 1 && return epsilon > 0
    for i in 1:(n - 1)
        a = ids[i]
        for j in (i + 1):n
            _bcs_dist(metric, a, ids[j]) < epsilon || return false
        end
    end
    return true
end

function _bcs_colors(cands::Vector{Int}, adj::Vector{BitSet})
    col = Dict{Int,Int}()
    ncolor = 0
    for v in cands
        used = BitSet()
        for u in adj[v]
            c = get(col, u, 0)
            c == 0 || push!(used, c)
        end
        c = 1
        while c in used
            c += 1
        end
        col[v] = c
        ncolor = max(ncolor, c)
    end
    return ncolor
end

function _bcs_expand!(best::Vector{Int}, clique::Vector{Int}, cands::Vector{Int}, adj::Vector{BitSet}, allow_one::Bool)
    if length(clique) > length(best) && (length(clique) >= 2 || allow_one)
        empty!(best)
        append!(best, clique)
    end
    isempty(cands) && return nothing
    length(clique) + _bcs_colors(cands, adj) <= length(best) && return nothing
    for t in eachindex(cands)
        length(clique) + (length(cands) - t + 1) <= length(best) && return nothing
        v = cands[t]
        push!(clique, v)
        nxt = Int[]
        nbrs = adj[v]
        for j in (t + 1):length(cands)
            cands[j] in nbrs && push!(nxt, cands[j])
        end
        _bcs_expand!(best, clique, nxt, adj, allow_one)
        pop!(clique)
    end
    return nothing
end

"""
Largest subset whose diameter is strictly below `epsilon`. Indices are
positions into `ids`. A one-vertex subset is allowed only when `epsilon > 0`,
because a single right vertex has diameter 0.
"""
function _bcs_largest(ids::Vector{Int}, epsilon::Real, metric)
    n = length(ids)
    n == 0 && return Int[]
    allow_one = epsilon > 0
    n == 1 && return allow_one ? copy(ids) : Int[]
    adj = [BitSet() for _ in 1:n]
    for i in 1:(n - 1)
        for j in (i + 1):n
            if _bcs_dist(metric, ids[i], ids[j]) < epsilon
                push!(adj[i], j)
                push!(adj[j], i)
            end
        end
    end
    best = Int[]
    _bcs_expand!(best, Int[], collect(1:n), adj, allow_one)
    return sort!(ids[best])
end

"""Left vertices of `left` adjacent to every vertex of `right`."""
function _bcs_anchor(g::BipartiteGraph, left::Vector{Int}, right::Vector{Int})
    kept = Int[]
    for u in left
        haskey(g.adjU, u) || continue
        nbrs = g.adjU[u]
        ok = true
        for v in right
            v in nbrs || (ok = false; break)
        end
        ok && push!(kept, u)
    end
    return kept
end

function _bcs_side(ids, present::Dict{Int,Set{Int}})
    out = Int[]
    seen = Set{Int}()
    for id in ids
        x = Int(id)
        haskey(present, x) || return nothing
        x in seen && return nothing
        push!(seen, x)
        push!(out, x)
    end
    return out
end

"""
    score_biclique(g, left, right, theta, epsilon, metric) -> NamedTuple

How a biclique `(left, right)` reads as an ε-cohesive core. The biclique
search that produced the pair does not see `theta`, `epsilon`, or `metric`.

A biclique whose right side has at least `theta` vertices and diameter
strictly below `epsilon` is a core: every left vertex can use that same
right side as its witness. The core optimum is therefore at least `|left|`
whenever the raw test passes.

`raw_ok` is that test on the pair as returned. `raw_size` is `|left|` when
it passes and 0 otherwise. A single right vertex has diameter 0, so it
passes only when `epsilon > 0`.

`repaired_right` is a largest subset of `right` that still has at least
`theta` vertices and diameter strictly below `epsilon`. `repaired_left` is
every vertex of `left` adjacent to all of that subset, so the repaired pair
is still a biclique. `repaired_size` is `|repaired_left|` when such a subset
exists and 0 otherwise. On a biclique the repaired left side is the original
left side, because every left vertex touches every right vertex.

The diameter test is strict, matching `M < ε`. `metric(r1, r2)` is read in
identifier order, the same convention as the cohesive-core search.
"""
function score_biclique(g::BipartiteGraph, left, right, theta::Integer, epsilon::Real, metric)
    theta < 1 && throw(ArgumentError("theta must be ≥ 1"))
    eps = _bcs_epsilon(epsilon)
    L = _bcs_side(left, g.adjU)
    R = _bcs_side(right, g.adjV)
    if L === nothing || R === nothing
        return (;
            raw_ok=false, raw_size=0,
            repaired_ok=false, repaired_size=0,
            repaired_left=Int[], repaired_right=Int[],
        )
    end
    raw_ok = !isempty(L) && _bcs_close(R, eps, metric) && length(R) >= theta && is_biclique(g, L, R)
    raw_size = raw_ok ? length(L) : 0
    repaired_right = _bcs_largest(R, eps, metric)
    if length(repaired_right) < theta
        return (;
            raw_ok, raw_size,
            repaired_ok=false, repaired_size=0,
            repaired_left=Int[], repaired_right=Int[],
        )
    end
    repaired_left = sort!(_bcs_anchor(g, L, repaired_right))
    repaired_ok = !isempty(repaired_left)
    return (;
        raw_ok, raw_size,
        repaired_ok,
        repaired_size=repaired_ok ? length(repaired_left) : 0,
        repaired_left=repaired_ok ? repaired_left : Int[],
        repaired_right=repaired_ok ? repaired_right : Int[],
    )
end
