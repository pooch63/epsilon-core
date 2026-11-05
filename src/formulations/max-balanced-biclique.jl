#=
Maximum balanced biclique by TSGR-MBBP (Zhou and Hao, EAAI 2019).

  include(joinpath(@__DIR__, "max-balanced-biclique.jl"))

  found = max_balanced_biclique(g)
  found.left
  found.right     # |left| == |right|
  found.balanced
  found.proven

The search grows slightly unbalanced bicliques with a tabu "push", then
deletes vertices of degree at most the best balanced size and solves the
small components that remain. `g` is not mutated. θ, ε, and any distance
are applied afterwards by `score_biclique` in max-edge-biclique.jl.

Li, Hao, and Wu (Computers & Operations Research, 2020) report better
best-known values on 19 instances with a k-swap search. The authors post
C++ for that solver. This file is the 2019 baseline, which is specified
in full and runs on the bipartite graph directly.

  julia --project=. test/max-balanced-biclique.jl
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

using Random

"""
    is_balanced_biclique(g, left, right) -> Bool

`left` and `right` have the same positive size, and every listed left vertex
is adjacent to every listed right vertex.
"""
function is_balanced_biclique(g::BipartiteGraph, left, right)
    L = collect(Int, left)
    R = collect(Int, right)
    length(L) == length(R) || return false
    isempty(L) && return false
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

function _mbb_copy(g::BipartiteGraph)
    adjU = Dict{Int,Set{Int}}(u => copy(nbrs) for (u, nbrs) in g.adjU)
    adjV = Dict{Int,Set{Int}}(v => copy(nbrs) for (v, nbrs) in g.adjV)
    return adjU, adjV
end

function _mbb_drop_u!(adjU, adjV, u::Int)
    nbrs = pop!(adjU, u, nothing)
    nbrs === nothing && return nothing
    for v in nbrs
        haskey(adjV, v) && delete!(adjV[v], u)
    end
    return nothing
end

function _mbb_drop_v!(adjU, adjV, v::Int)
    nbrs = pop!(adjV, v, nothing)
    nbrs === nothing && return nothing
    for u in nbrs
        haskey(adjU, u) && delete!(adjU[u], v)
    end
    return nothing
end

"""Delete every vertex of degree `≤ omega`, then repeat on the degrees that fall."""
function _mbb_peel!(adjU, adjV, omega::Int)
    while true
        drop_u = Int[]
        for (u, nbrs) in adjU
            length(nbrs) <= omega && push!(drop_u, u)
        end
        drop_v = Int[]
        for (v, nbrs) in adjV
            length(nbrs) <= omega && push!(drop_v, v)
        end
        isempty(drop_u) && isempty(drop_v) && return nothing
        for u in drop_u
            _mbb_drop_u!(adjU, adjV, u)
        end
        for v in drop_v
            _mbb_drop_v!(adjU, adjV, v)
        end
    end
end

function _mbb_min_degree(adjU, adjV)
    m = typemax(Int)
    for nbrs in values(adjU)
        m = min(m, length(nbrs))
    end
    for nbrs in values(adjV)
        m = min(m, length(nbrs))
    end
    return m
end

function _mbb_common(owners, adj::Dict{Int,Set{Int}}, exclude::Set{Int})
    smallest = nothing
    smallest_n = typemax(Int)
    for o in owners
        n = length(adj[o])
        if smallest === nothing || n < smallest_n || (n == smallest_n && o < smallest)
            smallest = o
            smallest_n = n
        end
    end
    smallest === nothing && return Int[]
    acc = copy(adj[smallest])
    for o in owners
        o == smallest && continue
        intersect!(acc, adj[o])
        isempty(acc) && return Int[]
    end
    for x in exclude
        delete!(acc, x)
    end
    return sort!(collect(acc))
end

"""
Alternate one vertex at a time, starting from a random vertex. The paper's
construction stops when the side whose turn it is has no common neighbor
left, so the two sides differ by at most one vertex.
"""
function _mbb_init(rng::AbstractRNG, adjU, adjV)
    verts = Tuple{Int,Bool}[]
    for u in sort!(collect(keys(adjU)))
        push!(verts, (u, true))
    end
    for v in sort!(collect(keys(adjV)))
        push!(verts, (v, false))
    end
    isempty(verts) && return Set{Int}(), Set{Int}()
    id, on_left = verts[rand(rng, eachindex(verts))]
    X = Set{Int}()
    Y = Set{Int}()
    if on_left
        push!(X, id)
        next_right = true
    else
        push!(Y, id)
        next_right = false
    end
    while true
        cands = next_right ? _mbb_common(X, adjU, Y) : _mbb_common(Y, adjV, X)
        isempty(cands) && break
        v = cands[rand(rng, eachindex(cands))]
        if next_right
            push!(Y, v)
        else
            push!(X, v)
        end
        next_right = !next_right
    end
    return X, Y
end

"""`tt(α, l) = max(7, α · uniform{0,…,l})`, truncated toward zero."""
function _mbb_tenure(rng::AbstractRNG, alpha::Float64, nside::Int)
    span = max(nside, 0)
    return max(7, floor(Int, alpha * rand(rng, 0:span)))
end

"""Balanced-size change of adding one vertex and dropping `missing` on the other side."""
@inline function _mbb_delta(to_n::Int, from_n::Int, missing::Int)
    return min(to_n + 1, from_n - missing) - min(to_n, from_n)
end

function _mbb_offer!(expand, plateau, id::Int, on_left::Bool, missing::Int, to_n::Int, from_n::Int,
                     iter::Int, tabu, best_bal::Int)
    delta = _mbb_delta(to_n, from_n, missing)
    free = get(tabu, (id, on_left), 0) <= iter
    # Equation (6) prints δ > 1. One push raises the balanced size by at most 1,
    # and Section 3.4.2 describes the improving set as δ > 0.
    if delta > 0 && (free || min(to_n, from_n) + 1 > best_bal)
        push!(expand, (id, on_left))
    elseif delta == 0 && free
        push!(plateau, (id, on_left))
    end
    return nothing
end

"""
Candidates that miss at most one vertex of the opposite side. `on_left` is
the side the candidate would join.
"""
function _mbb_scan!(expand, plateau, opposite::Set{Int}, opp_adj, own::Set{Int}, on_left::Bool,
                    iter::Int, tabu, best_bal::Int)
    from_n = length(opposite)
    from_n == 0 && return nothing
    to_n = length(own)
    counts = Dict{Int,Int}()
    for y in opposite
        for u in opp_adj[y]
            u in own && continue
            counts[u] = get(counts, u, 0) + 1
        end
    end
    need = from_n - 1
    for (u, c) in counts
        c >= need || continue
        _mbb_offer!(expand, plateau, u, on_left, from_n - c, to_n, from_n, iter, tabu, best_bal)
    end
    return nothing
end

function _mbb_candidates(adjU, adjV, X, Y, iter, tabu, best_bal)
    expand = Tuple{Int,Bool}[]
    plateau = Tuple{Int,Bool}[]
    _mbb_scan!(expand, plateau, Y, adjV, X, true, iter, tabu, best_bal)
    _mbb_scan!(expand, plateau, X, adjU, Y, false, iter, tabu, best_bal)
    sort!(expand)
    sort!(plateau)
    return expand, plateau
end

function _mbb_rand_id(rng::AbstractRNG, side::Set{Int})
    ids = sort!(collect(side))
    return ids[rand(rng, eachindex(ids))]
end

function _mbb_push!(X, Y, v::Int, on_left::Bool, adjU, adjV)
    expelled = Int[]
    if on_left
        push!(X, v)
        nbrs = adjU[v]
        for y in Y
            y in nbrs || push!(expelled, y)
        end
        for y in expelled
            delete!(Y, y)
        end
    else
        push!(Y, v)
        nbrs = adjV[v]
        for x in X
            x in nbrs || push!(expelled, x)
        end
        for x in expelled
            delete!(X, x)
        end
    end
    return expelled
end

"""
Constraint-based tabu search. Pushes add a vertex and drop at most the one
non-neighbor the candidate is allowed to miss. When the two sides differ by
more than 2, vertices are dropped from the larger side until the sides match.
Tabu keys carry the side because the two parts may reuse integer ids.
"""
function _mbb_cbts(rng::AbstractRNG, adjU, adjV, X::Set{Int}, Y::Set{Int}, depth::Int, alpha::Float64)
    bestX = copy(X)
    bestY = copy(Y)
    best_bal = min(length(X), length(Y))
    tabu = Dict{Tuple{Int,Bool},Int}()
    for iter in 0:depth
        expand, plateau = _mbb_candidates(adjU, adjV, X, Y, iter, tabu, best_bal)
        moved = nothing
        if !isempty(expand)
            moved = expand[rand(rng, eachindex(expand))]
        elseif !isempty(plateau)
            moved = plateau[rand(rng, eachindex(plateau))]
        end
        if moved !== nothing
            v, on_left = moved
            expelled = _mbb_push!(X, Y, v, on_left, adjU, adjV)
            side_n = on_left ? length(Y) : length(X)
            for u in expelled
                tabu[(u, !on_left)] = iter + _mbb_tenure(rng, alpha, side_n)
            end
        end
        if abs(length(X) - length(Y)) > 2
            while length(X) > length(Y)
                u = _mbb_rand_id(rng, X)
                delete!(X, u)
                tabu[(u, true)] = iter + _mbb_tenure(rng, alpha, length(X))
            end
            while length(Y) > length(X)
                u = _mbb_rand_id(rng, Y)
                delete!(Y, u)
                tabu[(u, false)] = iter + _mbb_tenure(rng, alpha, length(Y))
            end
        end
        bal = min(length(X), length(Y))
        if bal > best_bal
            best_bal = bal
            bestX = copy(X)
            bestY = copy(Y)
        end
    end
    return bestX, bestY
end

function _mbb_components(adjU, adjV)
    seen_u = Set{Int}()
    seen_v = Set{Int}()
    comps = Tuple{Vector{Int},Vector{Int}}[]
    function grow(start::Int, start_left::Bool)
        cu = Int[]
        cv = Int[]
        stack = Tuple{Int,Bool}[(start, start_left)]
        if start_left
            push!(seen_u, start)
        else
            push!(seen_v, start)
        end
        while !isempty(stack)
            id, on_left = pop!(stack)
            if on_left
                push!(cu, id)
                for v in adjU[id]
                    v in seen_v && continue
                    push!(seen_v, v)
                    push!(stack, (v, false))
                end
            else
                push!(cv, id)
                for u in adjV[id]
                    u in seen_u && continue
                    push!(seen_u, u)
                    push!(stack, (u, true))
                end
            end
        end
        return sort!(cu), sort!(cv)
    end
    for u in sort!(collect(keys(adjU)))
        u in seen_u && continue
        push!(comps, grow(u, true))
    end
    for v in sort!(collect(keys(adjV)))
        v in seen_v && continue
        push!(comps, grow(v, false))
    end
    sort!(comps, by=c -> (length(c[1]) + length(c[2]), isempty(c[1]) ? typemax(Int) : c[1][1]))
    return comps
end

function _mbb_cap(a::Vector{Int}, b::Vector{Int})
    out = Int[]
    i = 1
    j = 1
    while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            push!(out, a[i])
            i += 1
            j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
    return out
end

mutable struct _MbbExact
    adj_u::Vector{Vector{Int}}
    adj_v::Vector{Vector{Int}}
    u_ids::Vector{Int}
    v_ids::Vector{Int}
    lb::Int
    best_u::Vector{Int}
    best_v::Vector{Int}
    deadline::Float64
    timed_out::Bool
end

"""
Carraghan–Pardalos enumeration adapted to alternate sides (Appendix A).
`len_a` is the balanced size whenever the side about to be extended has no
candidates left: that side is the smaller one, or the two sides are equal.
Vertices on a side are tried in increasing position, and a later call only
receives the positions that follow the one just chosen.
"""
function _mbb_bb!(st::_MbbExact, len_a::Int, len_b::Int, Ca::Vector{Int}, Cb::Vector{Int}, a_is_u::Bool, Au::Vector{Int}, Bu::Vector{Int})
    st.timed_out && return nothing
    if time() > st.deadline
        st.timed_out = true
        return nothing
    end
    if isempty(Ca)
        if len_a > st.lb
            st.lb = len_a
            st.best_u = copy(Au)
            st.best_v = copy(Bu)
        end
        return nothing
    end
    for t in eachindex(Ca)
        if len_a + (length(Ca) - t + 1) <= st.lb
            return nothing
        end
        i = Ca[t]
        new_cb = Ca[(t + 1):end]
        if a_is_u
            push!(Au, st.u_ids[i])
            new_ca = _mbb_cap(Cb, st.adj_u[i])
            _mbb_bb!(st, len_b, len_a + 1, new_ca, new_cb, false, Au, Bu)
            pop!(Au)
        else
            push!(Bu, st.v_ids[i])
            new_ca = _mbb_cap(Cb, st.adj_v[i])
            _mbb_bb!(st, len_b, len_a + 1, new_ca, new_cb, true, Au, Bu)
            pop!(Bu)
        end
        st.timed_out && return nothing
    end
    return nothing
end

"""
Exact balanced biclique of one component, or any biclique better than
`omega` found before `seconds` elapses. The third value is false when the
clock cuts the search off; the component is then not safe to delete.
"""
function _mbb_exact(adjU, adjV, cu::Vector{Int}, cv::Vector{Int}, omega::Int, seconds::Float64)
    (isempty(cu) || isempty(cv)) && return Int[], Int[], true
    u_index = Dict{Int,Int}(id => i for (i, id) in enumerate(cu))
    v_index = Dict{Int,Int}(id => i for (i, id) in enumerate(cv))
    adj_u = [Int[] for _ in cu]
    adj_v = [Int[] for _ in cv]
    for (i, u) in enumerate(cu)
        for v in adjU[u]
            j = get(v_index, v, 0)
            j == 0 && continue
            push!(adj_u[i], j)
        end
        sort!(adj_u[i])
    end
    for (j, v) in enumerate(cv)
        for u in adjV[v]
            i = get(u_index, u, 0)
            i == 0 && continue
            push!(adj_v[j], i)
        end
        sort!(adj_v[j])
    end
    st = _MbbExact(adj_u, adj_v, cu, cv, omega, Int[], Int[], time() + seconds, false)
    _mbb_bb!(st, 0, 0, collect(1:length(cu)), collect(1:length(cv)), true, Int[], Int[])
    return st.best_u, st.best_v, !st.timed_out
end

"""Keep the smallest ids. At most two vertices are excess for a pair this search stores."""
function _mbb_balance(left, right)
    L = sort!(collect(Int, left))
    R = sort!(collect(Int, right))
    n = min(length(L), length(R))
    return L[1:n], R[1:n]
end

function _mbb_preset(preset::Symbol)
    if preset == :sparse
        return 100, 1.74, 500
    elseif preset == :dense
        return 1000, 0.30, 100
    end
    throw(ArgumentError("preset must be :sparse or :dense"))
end

"""
    max_balanced_biclique(g; preset=:sparse, restarts=20, depth=nothing,
                          alpha=nothing, component_limit=nothing,
                          exact_seconds=10, seconds=nothing, seed=1) -> NamedTuple

A balanced biclique `(left, right)` of `g`, with `|left| == |right| == balanced`.
`proven` is true when a part of the reduced graph is no larger than that
size, which is an upper bound, so the size is optimal.

Each restart builds a biclique by alternating random additions, then runs
`depth + 1` tabu pushes. A push adds a vertex that is adjacent to every
vertex of the opposite side, or to all but one of them, and drops that
missing vertex. Improving pushes are preferred. The other side's expelled
vertex is tabu for `max(7, floor(alpha * U{0,…,|side|}))` iterations, and a
tabu vertex may still be pushed when the push would beat the best biclique
of the current restart. If the sides differ by more than 2, the larger side
is thinned until they match.

Whenever the best balanced size is at least the minimum degree, vertices at
or below that degree are deleted. A remaining component of at most
`component_limit` vertices is solved by the branch-and-bound in Appendix A,
then deleted. A component that exceeds `exact_seconds` stays in the graph.
`preset=:sparse` uses the KONECT setting `(depth=100, alpha=1.74,
component_limit=500)`. `preset=:dense` uses `(1000, 0.30, 100)`. An explicit
`depth`, `alpha`, or `component_limit` replaces the preset value.

`seconds` is an optional wall-clock budget for the whole run. The search
stops early once the size is proven. The final equal-size pair keeps the
smallest identifiers on the side that was larger. `g` is not mutated.

Zhou and Hao allow the working pair to be short two vertices on one side.
The vertex-biclique optimum has no such balance requirement, so a comparison
against this baseline should expect a different shape. Neither baseline reads
θ, ε, or a distance; score the returned pair with `score_biclique`.
"""
function max_balanced_biclique(g::BipartiteGraph;
                               preset::Symbol=:sparse,
                               restarts::Integer=20,
                               depth::Union{Nothing,Integer}=nothing,
                               alpha::Union{Nothing,Real}=nothing,
                               component_limit::Union{Nothing,Integer}=nothing,
                               exact_seconds::Real=10,
                               seconds::Union{Nothing,Real}=nothing,
                               seed::Integer=1)
    restarts < 0 && throw(ArgumentError("restarts must be ≥ 0"))
    exact_seconds < 0 && throw(ArgumentError("exact_seconds must be ≥ 0"))
    if seconds !== nothing && seconds < 0
        throw(ArgumentError("seconds must be ≥ 0"))
    end
    default_depth, default_alpha, default_limit = _mbb_preset(preset)
    depth_i = depth === nothing ? default_depth : Int(depth)
    alpha_f = alpha === nothing ? default_alpha : Float64(alpha)
    limit_i = component_limit === nothing ? default_limit : Int(component_limit)
    depth_i < 0 && throw(ArgumentError("depth must be ≥ 0"))
    alpha_f < 0 && throw(ArgumentError("alpha must be ≥ 0"))
    limit_i < 0 && throw(ArgumentError("component_limit must be ≥ 0"))
    adjU, adjV = _mbb_copy(g)
    best_left = Int[]
    best_right = Int[]
    omega = 0
    proven = isempty(adjU) || isempty(adjV)
    rng = Xoshiro(Int(seed))
    started = time()
    budget = seconds === nothing ? Inf : Float64(seconds)
    for _ in 1:Int(restarts)
        proven && break
        time() - started >= budget && break
        if isempty(adjU) || isempty(adjV)
            proven = true
            break
        end
        X, Y = _mbb_init(rng, adjU, adjV)
        X, Y = _mbb_cbts(rng, adjU, adjV, X, Y, depth_i, alpha_f)
        bal = min(length(X), length(Y))
        if bal > omega
            omega = bal
            best_left = collect(X)
            best_right = collect(Y)
        end
        while !isempty(adjU) && !isempty(adjV) && omega >= _mbb_min_degree(adjU, adjV)
            _mbb_peel!(adjU, adjV, omega)
            for (cu, cv) in _mbb_components(adjU, adjV)
                length(cu) + length(cv) <= limit_i || continue
                found_u, found_v, done = _mbb_exact(adjU, adjV, cu, cv, omega, Float64(exact_seconds))
                found_bal = min(length(found_u), length(found_v))
                if found_bal > omega
                    omega = found_bal
                    best_left = found_u
                    best_right = found_v
                end
                done && (_mbb_delete!(adjU, adjV, cu, cv))
            end
            if length(adjU) <= omega || length(adjV) <= omega
                proven = true
                break
            end
        end
    end
    if isempty(adjU) || isempty(adjV) || length(adjU) <= omega || length(adjV) <= omega
        proven = true
    end
    left, right = _mbb_balance(best_left, best_right)
    return (; left, right, balanced=length(left), proven)
end

function _mbb_delete!(adjU, adjV, cu, cv)
    for u in cu
        _mbb_drop_u!(adjU, adjV, u)
    end
    for v in cv
        _mbb_drop_v!(adjU, adjV, v)
    end
    return nothing
end
