#=
Greedy search for a large ε-cohesive core.

  include(joinpath(@__DIR__, "heuristic.jl"))
  include(joinpath(@__DIR__, "distance.jl"))

  core = cohesive_core(g, 2, 0.5, overlap_metric(g))
  core = cohesive_core(g, 2, 0.5, diffusion_metric(g; time=2))
  core.left     # left vertices with ≥ theta neighbors in the group
  core.right    # all-pairs-similar right vertices that witness the core

Two right vertices are similar when `metric` returns a distance ≤ `epsilon`,
the same test `map_graph` uses. The metric is read on the graph passed in.
Peeling runs on a copy, so a metric closed over `g` (such as `overlap_metric`)
keeps the original neighborhoods.

  julia --project=. bin/heuristic.jl konect-small/euroroad --theta=2 --epsilon=0.5
  julia --project=. test/heuristic.jl
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

function _require_epsilon(epsilon::Real)
    if !(epsilon >= 0 && isfinite(epsilon))
        throw(ArgumentError("epsilon must be finite and ≥ 0"))
    end
    return epsilon
end

"""Structural copy. The original stays intact for a metric closed over it."""
function _copy_graph(g::BipartiteGraph{Nothing})
    h = BipartiteGraph{Nothing}()
    for (u, nbrs) in g.adjU
        h.adjU[u] = copy(nbrs)
    end
    for (v, nbrs) in g.adjV
        h.adjV[v] = copy(nbrs)
    end
    return h
end

@inline function _pair_distance(metric, r1::Int, r2::Int)
    return r1 <= r2 ? metric(r1, r2) : metric(r2, r1)
end

"""
    peel_core!(g, theta) -> NamedTuple

Drop every left vertex of degree `< theta`. Those vertices can never belong
to a core. Then drop every right vertex that has no left neighbor left.
Repeat until a pass removes nothing. Mutates `g`.
"""
function peel_core!(g::BipartiteGraph, theta::Int)
    theta < 1 && throw(ArgumentError("theta must be ≥ 1"))
    removed_left = 0
    removed_right = 0
    sweeps = 0
    while true
        changed = false
        drop_u = Int[]
        for (u, nbrs) in g.adjU
            length(nbrs) < theta && push!(drop_u, u)
        end
        for u in drop_u
            rem_u!(g, u) && (removed_left += 1; changed = true)
        end
        drop_v = Int[]
        for (v, nbrs) in g.adjV
            isempty(nbrs) && push!(drop_v, v)
        end
        for v in drop_v
            rem_v!(g, v) && (removed_right += 1; changed = true)
        end
        changed || return (; removed_left, removed_right, sweeps)
        sweeps += 1
    end
end

"""Right vertices joined when `metric ≤ epsilon`. No self-loops."""
function _similar_adj(ids::Vector{Int}, epsilon::Real, metric)
    adj = Dict{Int,Set{Int}}(id => Set{Int}() for id in ids)
    n = length(ids)
    for i in 1:n
        r1 = ids[i]
        for j in (i + 1):n
            r2 = ids[j]
            if _pair_distance(metric, r1, r2) <= epsilon
                push!(adj[r1], r2)
                push!(adj[r2], r1)
            end
        end
    end
    return adj
end

"""
Score of adding `c`. A left neighbor that already has `hits` group-neighbors,
and is still short of `theta`, contributes `hits + 1`. Vertices already at
`theta` contribute nothing: another neighbor does not grow the core. A vertex
one step from `theta` therefore outweighs one that has none.
"""
function _candidate_score(g::BipartiteGraph, c::Int, hits::Dict{Int,Int}, theta::Int)
    score = 0
    for u in g.adjV[c]
        h = get(hits, u, 0)
        h >= theta && continue
        score += h + 1
    end
    return score
end

"""Highest score, then more left neighbors, then the smaller id."""
function _choose_candidate(g::BipartiteGraph, cands::Set{Int}, hits::Dict{Int,Int}, theta::Int)
    best = nothing
    best_key = nothing
    for c in cands
        score = _candidate_score(g, c, hits, theta)
        score == 0 && continue
        key = (score, length(g.adjV[c]), -c)
        if best_key === nothing || key > best_key
            best_key = key
            best = c
        end
    end
    return best
end

function _alive_count(hits::Dict{Int,Int}, pool::Dict{Int,Int}, theta::Int)
    n = 0
    seen = Set{Int}()
    for (u, h) in hits
        push!(seen, u)
        h + get(pool, u, 0) >= theta && (n += 1)
    end
    for (u, p) in pool
        u in seen && continue
        p >= theta && (n += 1)
    end
    return n
end

function _move_into_group!(g::BipartiteGraph, c::Int, cands::Set{Int}, hits::Dict{Int,Int}, pool::Dict{Int,Int})
    delete!(cands, c)
    for u in g.adjV[c]
        hits[u] = get(hits, u, 0) + 1
        pool[u] -= 1
    end
    return nothing
end

function _shrink_pool!(g::BipartiteGraph, cands::Set{Int}, keep::Set{Int}, pool::Dict{Int,Int})
    drop = Int[]
    for c in cands
        c in keep || push!(drop, c)
    end
    for c in drop
        delete!(cands, c)
        for u in g.adjV[c]
            pool[u] -= 1
        end
    end
    return nothing
end

function _core_vertices(hits::Dict{Int,Int}, theta::Int)
    core = Int[]
    for (u, h) in hits
        h >= theta && push!(core, u)
    end
    sort!(core)
    return core
end

"""
Grow one all-pairs-similar group from `seed`. Returns `nothing` when the
number of left vertices that could still reach `theta` is no better than
`best_size`.
"""
function _grow_from(g::BipartiteGraph, seed::Int, theta::Int, similar::Dict{Int,Set{Int}}, best_size::Int)
    group = Int[seed]
    hits = Dict{Int,Int}()
    for u in g.adjV[seed]
        hits[u] = 1
    end
    cands = copy(similar[seed])
    pool = Dict{Int,Int}()
    for c in cands
        for u in g.adjV[c]
            pool[u] = get(pool, u, 0) + 1
        end
    end
    _alive_count(hits, pool, theta) <= best_size && return nothing
    while true
        pick = _choose_candidate(g, cands, hits, theta)
        pick === nothing && break
        push!(group, pick)
        _move_into_group!(g, pick, cands, hits, pool)
        _shrink_pool!(g, cands, similar[pick], pool)
        _alive_count(hits, pool, theta) <= best_size && return nothing
    end
    return (; left=_core_vertices(hits, theta), right=sort!(group))
end

"""
    cohesive_core(g, theta, epsilon, metric) -> NamedTuple

A large ε-cohesive core of `g`, witnessed by an all-pairs-similar right group.

Seeds are the right vertices that survive [`peel_core!`](@ref), tried
highest left-degree first. The group starts as the seed. The candidate pool
starts as every right vertex similar to the seed, and after each addition
shrinks to the vertices similar to every member, so the group stays
all-pairs similar. The vertex added next is the one that helps the most left
vertices still short of `theta`, weighting a vertex with `hits` group
neighbors by `hits + 1`.

A left vertex is alive while its neighbors in the group plus its neighbors
still in the pool could reach `theta`. That count is an upper bound on the
core this seed can produce. A seed whose bound is no larger than the best
core so far is abandoned.

`left` is every surviving left vertex with at least `theta` neighbors in
`right`. `pruned` is how many seeds were abandoned by the bound. `g` is not
mutated. `metric(r1, r2) ≤ epsilon` means similar; pass `overlap_metric(g)`
to use the overlap coefficient on the original neighborhoods.
"""
function cohesive_core(g::BipartiteGraph{Nothing}, theta::Int, epsilon::Real, metric)
    theta < 1 && throw(ArgumentError("theta must be ≥ 1"))
    _require_epsilon(epsilon)
    work = _copy_graph(g)
    peel_core!(work, theta)
    ids = collect(keys(work.adjV))
    sort!(ids, by=v -> (-length(work.adjV[v]), v))
    similar = _similar_adj(ids, epsilon, metric)
    best_left = Int[]
    best_right = Int[]
    pruned = 0
    for seed in ids
        grown = _grow_from(work, seed, theta, similar, length(best_left))
        if grown === nothing
            pruned += 1
            continue
        end
        if length(grown.left) > length(best_left)
            best_left = grown.left
            best_right = grown.right
        end
    end
    return (; left=best_left, right=best_right, pruned)
end

"""
    is_cohesive_core(g, left, right, theta, epsilon, metric) -> Bool

`right` is pairwise at distance ≤ `epsilon`, and `left` is exactly the set of
left vertices with at least `theta` neighbors in `right`.
"""
function is_cohesive_core(g::BipartiteGraph, left, right, theta::Int, epsilon::Real, metric)
    theta < 1 && return false
    _require_epsilon(epsilon)
    ids = collect(Int, right)
    for i in eachindex(ids)
        r1 = ids[i]
        haskey(g.adjV, r1) || return false
        for j in (i + 1):length(ids)
            _pair_distance(metric, r1, ids[j]) <= epsilon || return false
        end
    end
    wanted = Set{Int}(ids)
    found = Int[]
    for (u, nbrs) in g.adjU
        n = 0
        for v in nbrs
            v in wanted && (n += 1)
        end
        n >= theta && push!(found, u)
    end
    sort!(found)
    return found == sort!(collect(Int, left))
end
