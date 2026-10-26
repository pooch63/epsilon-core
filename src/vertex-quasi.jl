#=
Greedy search for a large vertex δ-quasi-biclique.

  include(joinpath(@__DIR__, "vertex-quasi.jl"))
  qb = vertex_quasi_biclique(g, 0.2; tau=2)
  qb.left
  qb.right

Liu, Li, and Wang, as outlined by Ignatov, Ivanova, and Zamaletdinova
(arXiv:2002.09880, Definition 7 and Section 2.3). A pair of sets is a
δ-quasi-biclique when every vertex keeps at least a (1-δ) fraction of the
opposite side as neighbors. δ = 0 is a biclique. The paper's density
parameter is γ = 1-δ; the vertex rule is the stronger of the two.

One pass grows the left side from empty. The right side starts as every
right vertex and is pruned after each addition. The pass is repeated with
the sides swapped, and the larger pair is kept. `g` is not mutated.

  julia --project=. test/vertex-quasi.jl
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

function _require_delta(delta::Real)
    if !(delta >= 0 && delta <= 0.5 && isfinite(delta))
        throw(ArgumentError("delta must be finite and in [0, 0.5]"))
    end
    return Float64(delta)
end

"""`deg` neighbors out of `n` miss at most a `delta` fraction. An empty
opposite side holds vacuously: the threshold is 0."""
@inline function _fraction_ok(deg::Int, n::Int, delta::Float64)
    n == 0 && return true
    return (n - deg) <= delta * n + 1e-8
end

function _count_in(nbrs::Set{Int}, chosen::Set{Int})
    if length(chosen) < length(nbrs)
        n = 0
        for x in chosen
            x in nbrs && (n += 1)
        end
        return n
    end
    n = 0
    for x in nbrs
        x in chosen && (n += 1)
    end
    return n
end

function _pair_ok(grown::Set{Int}, kept::Set{Int}, grown_deg, kept_deg, delta::Float64)
    ng, nk = length(grown), length(kept)
    (ng == 0 || nk == 0) && return false
    for u in grown
        _fraction_ok(grown_deg[u], nk, delta) || return false
    end
    for v in kept
        _fraction_ok(kept_deg[v], ng, delta) || return false
    end
    return true
end

"""Highest degree into the current opposite side, then the smaller id."""
function _pick_degree(ids, skip::Set{Int}, deg)
    best = nothing
    best_deg = -1
    for u in ids
        u in skip && continue
        d = deg[u]
        if best === nothing || d > best_deg || (d == best_deg && u < best)
            best = u
            best_deg = d
        end
    end
    return best
end

function _add!(u::Int, own::Set{Int}, other_deg, nbrs)
    push!(own, u)
    for v in nbrs
        other_deg[v] += 1
    end
    return nothing
end

function _drop!(victims, own::Set{Int}, other_deg, adj)
    for u in victims
        delete!(own, u)
        for v in adj[u]
            other_deg[v] -= 1
        end
    end
    return nothing
end

"""Drop kept vertices that miss the degree rule against the grown side."""
function _prune_kept!(grown, kept, grown_deg, kept_deg, adj_kept, delta)
    ng = length(grown)
    drop = Int[]
    for v in kept
        _fraction_ok(kept_deg[v], ng, delta) || push!(drop, v)
    end
    _drop!(drop, kept, grown_deg, adj_kept)
    return nothing
end

"""Shed violators until both sides meet the rule. The grown side goes first,
so a vertex added too eagerly is removed before the opposite side is cut
back to match it."""
function _close!(grown, kept, grown_deg, kept_deg, adj_grown, adj_kept, delta)
    while true
        nk = length(kept)
        drop_g = Int[]
        for u in grown
            _fraction_ok(grown_deg[u], nk, delta) || push!(drop_g, u)
        end
        if !isempty(drop_g)
            _drop!(drop_g, grown, kept_deg, adj_grown)
            continue
        end
        ng = length(grown)
        drop_k = Int[]
        for v in kept
            _fraction_ok(kept_deg[v], ng, delta) || push!(drop_k, v)
        end
        isempty(drop_k) && return nothing
        _drop!(drop_k, kept, grown_deg, adj_kept)
    end
end

"""`u` may join `own` when the induced subgraph on the enlarged sets still
meets the degree rule. The opposite side is not pruned."""
function _can_add(nbrs, own_n::Int, other::Set{Int}, own_deg::Int, other_deg, delta)
    other_n = length(other)
    (own_n == 0 || other_n == 0) && return false
    _fraction_ok(own_deg, other_n, delta) || return false
    need_n = own_n + 1
    for v in other
        inc = v in nbrs ? 1 : 0
        _fraction_ok(other_deg[v] + inc, need_n, delta) || return false
    end
    return true
end

function _pick_addable(ids, own, other, own_deg, other_deg, adj_own, delta)
    best = nothing
    best_deg = -1
    own_n = length(own)
    for u in ids
        u in own && continue
        d = own_deg[u]
        if best !== nothing && (d < best_deg || (d == best_deg && u > best))
            continue
        end
        _can_add(adj_own[u], own_n, other, d, other_deg, delta) || continue
        best = u
        best_deg = d
    end
    return best
end

function _expand!(grown_ids, kept_ids, grown, kept, grown_deg, kept_deg, adj_grown, adj_kept, delta)
    while true
        added = false
        u = _pick_addable(grown_ids, grown, kept, grown_deg, kept_deg, adj_grown, delta)
        if u !== nothing
            _add!(u, grown, kept_deg, adj_grown[u])
            added = true
        end
        v = _pick_addable(kept_ids, kept, grown, kept_deg, grown_deg, adj_kept, delta)
        if v !== nothing
            _add!(v, kept, grown_deg, adj_kept[v])
            added = true
        end
        added || return nothing
    end
end

function _recount(grown, kept, adj_grown, adj_kept)
    grown_deg = Dict{Int,Int}(u => _count_in(nbrs, kept) for (u, nbrs) in adj_grown)
    kept_deg = Dict{Int,Int}(v => _count_in(nbrs, grown) for (v, nbrs) in adj_kept)
    return grown_deg, kept_deg
end

"""Remember `grown`/`kept` when it meets the rule and beats the stored pair.
`equal_ok` keeps the later pair on a tie."""
function _consider(best_g, best_k, best_n, grown, kept, grown_deg, kept_deg, delta, equal_ok::Bool)
    _pair_ok(grown, kept, grown_deg, kept_deg, delta) || return best_g, best_k, best_n
    n = length(grown) + length(kept)
    better = equal_ok ? n >= best_n : n > best_n
    better || return best_g, best_k, best_n
    return copy(grown), copy(kept), n
end

"""
One orientation. `adj_grown` starts empty and receives up to `tau` vertices;
`adj_kept` starts as the whole part and is pruned. Returns the two sets in
that order.
"""
function _quasi_pass(adj_grown::Dict{Int,Set{Int}}, adj_kept::Dict{Int,Set{Int}}, delta::Float64, tau::Int)
    if isempty(adj_grown) || isempty(adj_kept)
        return Int[], Int[]
    end
    grown_ids = collect(keys(adj_grown))
    kept_ids = collect(keys(adj_kept))
    grown = Set{Int}()
    kept = Set{Int}(kept_ids)
    grown_deg = Dict{Int,Int}(u => length(nbrs) for (u, nbrs) in adj_grown)
    kept_deg = Dict{Int,Int}(v => 0 for v in kept_ids)

    best_g = Set{Int}()
    best_k = Set{Int}()
    best_n = 0
    # An addition that breaks the rule stays in the working pair. The next
    # one can prune the opposite side and leave a larger valid pair than the
    # best saved so far, so the search does not roll back until the end.
    while length(grown) < tau
        isempty(kept) && break
        u = _pick_degree(grown_ids, grown, grown_deg)
        u === nothing && break
        _add!(u, grown, kept_deg, adj_grown[u])
        _prune_kept!(grown, kept, grown_deg, kept_deg, adj_kept, delta)
        best_g, best_k, best_n = _consider(best_g, best_k, best_n, grown, kept, grown_deg, kept_deg, delta, true)
    end
    _close!(grown, kept, grown_deg, kept_deg, adj_grown, adj_kept, delta)
    best_g, best_k, best_n = _consider(best_g, best_k, best_n, grown, kept, grown_deg, kept_deg, delta, false)
    best_n == 0 && return Int[], Int[]

    grown_deg, kept_deg = _recount(best_g, best_k, adj_grown, adj_kept)
    _expand!(grown_ids, kept_ids, best_g, best_k, grown_deg, kept_deg, adj_grown, adj_kept, delta)
    return sort!(collect(best_g)), sort!(collect(best_k))
end

"""
    vertex_quasi_biclique(g, delta; tau=1) -> NamedTuple

A large δ-quasi-biclique `(left, right)` of `g`.

Every returned left vertex has at least a `(1 - delta)` fraction of `right`
as neighbors, and every returned right vertex has at least a `(1 - delta)`
fraction of `left`. `delta` is in `[0, 0.5]`. `tau ≥ 1` is how many vertices
the first phase places on the side that starts empty.

That phase follows the outline in Section 2.3. The side starts empty, the
other starts full. Until the empty side has `tau` vertices, add the vertex
with the most neighbors on the current opposite side (ties: smaller id) and
drop every opposite vertex below the degree rule. The largest pair that
still meets the rule is kept, and a tie keeps the later pair. Vertices that
still miss the rule are then removed, grown side first; that repaired pair
replaces the saved one only when it is strictly larger. The second phase
adds a vertex of either side, grown side first, when the enlarged pair still
meets the rule.

The pass runs once in each direction. The pair with more vertices is
returned; a tie keeps the pass that grew the left side. `left` and `right`
are sorted. `g` is not mutated. A graph with no usable pair returns empty
sets.
"""
function vertex_quasi_biclique(g::BipartiteGraph, delta::Real; tau::Integer=1)
    tau < 1 && throw(ArgumentError("tau must be ≥ 1"))
    d = _require_delta(delta)
    t = Int(tau)
    left, right = _quasi_pass(g.adjU, g.adjV, d, t)
    rev_right, rev_left = _quasi_pass(g.adjV, g.adjU, d, t)
    if length(rev_left) + length(rev_right) > length(left) + length(right)
        return (; left=rev_left, right=rev_right)
    end
    return (; left, right)
end

"""
    is_vertex_quasi_biclique(g, left, right, delta) -> Bool

Both sides are nonempty, every listed vertex is on that side of `g`, and
each one meets the δ degree rule against the other side.
"""
function is_vertex_quasi_biclique(g::BipartiteGraph, left, right, delta::Real)
    d = _require_delta(delta)
    L = Set{Int}(collect(Int, left))
    R = Set{Int}(collect(Int, right))
    (isempty(L) || isempty(R)) && return false
    for u in L
        haskey(g.adjU, u) || return false
        _fraction_ok(_count_in(g.adjU[u], R), length(R), d) || return false
    end
    for v in R
        haskey(g.adjV, v) || return false
        _fraction_ok(_count_in(g.adjV[v], L), length(L), d) || return false
    end
    return true
end
