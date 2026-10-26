#=
A large k-biplex, by the deletion-and-expansion heuristic of Cui, Li, Dai,
and Wang, arXiv:2607.07419, Algorithm 3.

  include(joinpath(@__DIR__, "k-biplex.jl"))

  found = k_biplex(g, 1, 3)
  found.left    # one side of the biplex
  found.right   # the other side
  found.edges   # edges induced between them

  julia --project=. test/k-biplex.jl

A k-biplex is an induced biclique with at most k missing edges at every
vertex. Both sides must have at least theta vertices. The search does not
mutate g.
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

"""
    is_k_biplex(g, left, right, k) -> Bool

Every vertex of `left` has at most `k` non-neighbors in `right`, and every
vertex of `right` has at most `k` non-neighbors in `left`. Both collections
are free of duplicates and use only vertices of `g`.
"""
function is_k_biplex(g::BipartiteGraph, left, right, k::Integer)
    k < 0 && return false
    left_ids = collect(Int, left)
    right_ids = collect(Int, right)
    left_set = Set{Int}(left_ids)
    right_set = Set{Int}(right_ids)
    length(left_ids) == length(left_set) || return false
    length(right_ids) == length(right_set) || return false
    n_right = length(right_set)
    n_left = length(left_set)
    for u in left_ids
        haskey(g.adjU, u) || return false
        hits = 0
        for v in g.adjU[u]
            v in right_set && (hits += 1)
        end
        n_right - hits > k && return false
    end
    for v in right_ids
        haskey(g.adjV, v) || return false
        hits = 0
        for u in g.adjV[v]
            u in left_set && (hits += 1)
        end
        n_left - hits > k && return false
    end
    return true
end

# Min-heap of (degree, id, generation). Generation drops stale degree updates.
@inline function _kbp_before(a::NTuple{3,Int}, b::NTuple{3,Int})
    return a[1] < b[1] || (a[1] == b[1] && a[2] <= b[2])
end

function _kbp_hpush!(h::Vector{NTuple{3,Int}}, item::NTuple{3,Int})
    push!(h, item)
    i = length(h)
    while i > 1
        p = i >> 1
        _kbp_before(h[p], h[i]) && break
        h[p], h[i] = h[i], h[p]
        i = p
    end
    return h
end

function _kbp_hpop!(h::Vector{NTuple{3,Int}})
    top = h[1]
    last = pop!(h)
    if !isempty(h)
        h[1] = last
        i = 1
        n = length(h)
        while true
            l = i << 1
            l > n && break
            child = l
            r = l + 1
            if r <= n && _kbp_before(h[r], h[l])
                child = r
            end
            _kbp_before(h[i], h[child]) && break
            h[i], h[child] = h[child], h[i]
            i = child
        end
    end
    return top
end

function _kbp_peek!(h::Vector{NTuple{3,Int}}, gen::Dict{Int,Int})
    while !isempty(h)
        degree, id, generation = h[1]
        get(gen, id, -1) == generation && return degree
        _kbp_hpop!(h)
    end
    return 0
end

function _kbp_extract!(h::Vector{NTuple{3,Int}}, gen::Dict{Int,Int})
    while !isempty(h)
        _, id, generation = _kbp_hpop!(h)
        get(gen, id, -1) == generation && return id
    end
    error("k-biplex peel ran out of vertices")
end

function _kbp_drop_side!(id::Int, neighbors::Set{Int}, own_deg::Dict{Int,Int}, own_gen::Dict{Int,Int},
                         opp_deg::Dict{Int,Int}, opp_gen::Dict{Int,Int}, opp_heap::Vector{NTuple{3,Int}})
    delete!(own_deg, id)
    delete!(own_gen, id)
    for other in neighbors
        haskey(opp_deg, other) || continue
        opp_deg[other] -= 1
        opp_gen[other] += 1
        _kbp_hpush!(opp_heap, (opp_deg[other], other, opp_gen[other]))
    end
    return nothing
end

"""
Balanced greedy deletion. Remove a minimum-degree vertex of the larger side
(the right side when the sides have equal size) until the induced subgraph
is a k-biplex. Degree ties break toward the smaller id.
"""
function _kbp_peel(g::BipartiteGraph, k::Int)
    deg_u = Dict{Int,Int}()
    deg_v = Dict{Int,Int}()
    gen_u = Dict{Int,Int}()
    gen_v = Dict{Int,Int}()
    heap_u = NTuple{3,Int}[]
    heap_v = NTuple{3,Int}[]
    sizehint!(heap_u, length(g.adjU))
    sizehint!(heap_v, length(g.adjV))
    for (u, nbrs) in g.adjU
        degree = length(nbrs)
        deg_u[u] = degree
        gen_u[u] = 1
        _kbp_hpush!(heap_u, (degree, u, 1))
    end
    for (v, nbrs) in g.adjV
        degree = length(nbrs)
        deg_v[v] = degree
        gen_v[v] = 1
        _kbp_hpush!(heap_v, (degree, v, 1))
    end
    n_u = length(deg_u)
    n_v = length(deg_v)
    while n_u > 0 && n_v > 0
        min_u = _kbp_peek!(heap_u, gen_u)
        min_v = _kbp_peek!(heap_v, gen_v)
        min_u >= n_v - k && min_v >= n_u - k && break
        if n_u > n_v
            u = _kbp_extract!(heap_u, gen_u)
            _kbp_drop_side!(u, g.adjU[u], deg_u, gen_u, deg_v, gen_v, heap_v)
            n_u -= 1
        else
            v = _kbp_extract!(heap_v, gen_v)
            _kbp_drop_side!(v, g.adjV[v], deg_v, gen_v, deg_u, gen_u, heap_u)
            n_v -= 1
        end
    end
    edges = 0
    for degree in values(deg_u)
        edges += degree
    end
    return collect(keys(deg_u)), collect(keys(deg_v)), edges
end

mutable struct _KbpInst
    u_ids::Vector{Int}
    v_ids::Vector{Int}
    adj_u::Vector{Vector{Int}}
    adj_v::Vector{Vector{Int}}
    in_s_u::BitVector
    in_c_u::BitVector
    in_s_v::BitVector
    in_c_v::BitVector
    in_sc_u::BitVector
    in_sc_v::BitVector
    deg_sc_u::Vector{Int}
    deg_sc_v::Vector{Int}
    deg_s_u::Vector{Int}
    deg_s_v::Vector{Int}
    s_list_u::Vector{Int}
    s_list_v::Vector{Int}
    n_s_u::Int
    n_s_v::Int
    n_c_u::Int
    n_c_v::Int
    k::Int
    theta::Int
    mark_u::Vector{Int}
    mark_v::Vector{Int}
    mark_gen::Int
end

function _kbp_discard_u!(st::_KbpInst, i::Int)
    st.in_c_u[i] || return false
    st.in_c_u[i] = false
    st.in_sc_u[i] = false
    st.n_c_u -= 1
    for j in st.adj_u[i]
        st.in_sc_v[j] || continue
        st.deg_sc_v[j] -= 1
    end
    return true
end

function _kbp_discard_v!(st::_KbpInst, j::Int)
    st.in_c_v[j] || return false
    st.in_c_v[j] = false
    st.in_sc_v[j] = false
    st.n_c_v -= 1
    for i in st.adj_v[j]
        st.in_sc_u[i] || continue
        st.deg_sc_u[i] -= 1
    end
    return true
end

function _kbp_accept_u!(st::_KbpInst, i::Int)
    st.in_c_u[i] || return false
    st.in_c_u[i] = false
    st.in_s_u[i] = true
    st.n_c_u -= 1
    st.n_s_u += 1
    push!(st.s_list_u, i)
    for j in st.adj_u[i]
        st.in_sc_v[j] || continue
        st.deg_s_v[j] += 1
    end
    return true
end

function _kbp_accept_v!(st::_KbpInst, j::Int)
    st.in_c_v[j] || return false
    st.in_c_v[j] = false
    st.in_s_v[j] = true
    st.n_c_v -= 1
    st.n_s_v += 1
    push!(st.s_list_v, j)
    for i in st.adj_v[j]
        st.in_sc_u[i] || continue
        st.deg_s_u[i] += 1
    end
    return true
end

function _kbp_fresh!(st::_KbpInst)
    gen = st.mark_gen + 1
    if gen == typemax(Int)
        fill!(st.mark_u, 0)
        fill!(st.mark_v, 0)
        gen = 1
    end
    st.mark_gen = gen
    return gen
end

"""Adding `i` keeps the committed left-to-right biplex within k missing edges."""
function _kbp_can_add_u(st::_KbpInst, i::Int)
    st.n_s_v - st.deg_s_u[i] > st.k && return false
    gen = _kbp_fresh!(st)
    for j in st.adj_u[i]
        st.mark_v[j] = gen
    end
    for j in st.s_list_v
        st.mark_v[j] == gen && continue
        st.n_s_u - st.deg_s_v[j] + 1 > st.k && return false
    end
    return true
end

function _kbp_can_add_v(st::_KbpInst, j::Int)
    st.n_s_u - st.deg_s_v[j] > st.k && return false
    gen = _kbp_fresh!(st)
    for i in st.adj_v[j]
        st.mark_u[i] = gen
    end
    for i in st.s_list_u
        st.mark_u[i] == gen && continue
        st.n_s_v - st.deg_s_u[i] + 1 > st.k && return false
    end
    return true
end

"""Lemma 4: a candidate with fewer than theta − k neighbors in the pool can never meet theta."""
function _kbp_lemma4!(st::_KbpInst)
    need = st.theta - st.k
    need <= 0 && return false
    queue = Tuple{Bool,Int}[]
    for i in eachindex(st.u_ids)
        st.in_c_u[i] && st.deg_sc_u[i] < need && push!(queue, (true, i))
    end
    for j in eachindex(st.v_ids)
        st.in_c_v[j] && st.deg_sc_v[j] < need && push!(queue, (false, j))
    end
    changed = false
    head = 1
    while head <= length(queue)
        is_u, idx = queue[head]
        head += 1
        if is_u
            (st.in_c_u[idx] && st.deg_sc_u[idx] < need) || continue
            _kbp_discard_u!(st, idx)
            changed = true
            for j in st.adj_u[idx]
                st.in_c_v[j] && st.deg_sc_v[j] < need && push!(queue, (false, j))
            end
        else
            (st.in_c_v[idx] && st.deg_sc_v[idx] < need) || continue
            _kbp_discard_v!(st, idx)
            changed = true
            for i in st.adj_v[idx]
                st.in_c_u[i] && st.deg_sc_u[i] < need && push!(queue, (true, i))
            end
        end
    end
    return changed
end

"""Lemma 3: drop a candidate that would put the committed biplex over k missing edges."""
function _kbp_lemma3!(st::_KbpInst)
    drop_u = Int[]
    drop_v = Int[]
    for i in eachindex(st.u_ids)
        st.in_c_u[i] && !_kbp_can_add_u(st, i) && push!(drop_u, i)
    end
    for j in eachindex(st.v_ids)
        st.in_c_v[j] && !_kbp_can_add_v(st, j) && push!(drop_v, j)
    end
    for i in drop_u
        _kbp_discard_u!(st, i)
    end
    for j in drop_v
        _kbp_discard_v!(st, j)
    end
    return !isempty(drop_u) || !isempty(drop_v)
end

"""
Lemma 5: keep a candidate whose non-neighbors in the whole pool, and the
candidate itself, already have at most k non-neighbors there. Moving it into
the committed set cannot break any extension.
"""
function _kbp_lemma5!(st::_KbpInst)
    n_u = st.n_s_u + st.n_c_u
    n_v = st.n_s_v + st.n_c_v
    k = st.k
    bad_u = Int[]
    bad_v = Int[]
    for i in eachindex(st.u_ids)
        st.in_sc_u[i] && n_v - st.deg_sc_u[i] > k && push!(bad_u, i)
    end
    for j in eachindex(st.v_ids)
        st.in_sc_v[j] && n_u - st.deg_sc_v[j] > k && push!(bad_v, j)
    end
    gen_v = _kbp_fresh!(st)
    for j in bad_v
        st.mark_v[j] = gen_v
    end
    gen_u = _kbp_fresh!(st)
    for i in bad_u
        st.mark_u[i] = gen_u
    end
    add_u = Int[]
    add_v = Int[]
    n_bad_v = length(bad_v)
    n_bad_u = length(bad_u)
    for i in eachindex(st.u_ids)
        st.in_c_u[i] || continue
        n_v - st.deg_sc_u[i] > k && continue
        if n_bad_v > 0
            hits = 0
            for j in st.adj_u[i]
                st.mark_v[j] == gen_v && (hits += 1)
            end
            hits == n_bad_v || continue
        end
        push!(add_u, i)
    end
    for j in eachindex(st.v_ids)
        st.in_c_v[j] || continue
        n_u - st.deg_sc_v[j] > k && continue
        if n_bad_u > 0
            hits = 0
            for i in st.adj_v[j]
                st.mark_u[i] == gen_u && (hits += 1)
            end
            hits == n_bad_u || continue
        end
        push!(add_v, j)
    end
    for i in add_u
        _kbp_accept_u!(st, i)
    end
    for j in add_v
        _kbp_accept_v!(st, j)
    end
    return !isempty(add_u) || !isempty(add_v)
end

function _kbp_reduce!(st::_KbpInst)
    while true
        changed = false
        _kbp_lemma4!(st) && (changed = true)
        _kbp_lemma3!(st) && (changed = true)
        _kbp_lemma5!(st) && (changed = true)
        changed || return nothing
    end
end

"""Highest induced degree, then a left vertex, then the smaller id."""
function _kbp_pick(st::_KbpInst)
    best_key = nothing
    best_u = true
    best_idx = 0
    for i in eachindex(st.u_ids)
        st.in_c_u[i] || continue
        key = (st.deg_sc_u[i], 1, -st.u_ids[i])
        if best_key === nothing || key > best_key
            best_key = key
            best_u = true
            best_idx = i
        end
    end
    for j in eachindex(st.v_ids)
        st.in_c_v[j] || continue
        key = (st.deg_sc_v[j], 0, -st.v_ids[j])
        if best_key === nothing || key > best_key
            best_key = key
            best_u = false
            best_idx = j
        end
    end
    best_key === nothing && return nothing
    return best_u, best_idx
end

function _kbp_count_edges(st::_KbpInst)
    edges = 0
    for i in eachindex(st.u_ids)
        st.in_s_u[i] || continue
        for j in st.adj_u[i]
            st.in_s_v[j] && (edges += 1)
        end
    end
    return edges
end

function _kbp_members(st::_KbpInst)
    left = Int[]
    right = Int[]
    for i in eachindex(st.u_ids)
        st.in_s_u[i] && push!(left, st.u_ids[i])
    end
    for j in eachindex(st.v_ids)
        st.in_s_v[j] && push!(right, st.v_ids[j])
    end
    return left, right
end

function _kbp_local(g::BipartiteGraph, ui::Int, later::Set{Int}, k::Int, theta::Int)
    theta - k > 0 && length(g.adjU[ui]) < theta - k && return nothing
    length(later) + 1 < theta && return nothing
    hop2 = Set{Int}((ui,))
    for v in g.adjU[ui]
        for u in g.adjV[v]
            push!(hop2, u)
        end
    end
    candidates_u = Int[]
    for u in hop2
        u == ui && continue
        u in later && push!(candidates_u, u)
    end
    sort!(candidates_u)
    length(candidates_u) + 1 < theta && return nothing
    candidates_v = Set{Int}()
    for u in hop2
        for v in g.adjU[u]
            push!(candidates_v, v)
        end
    end
    length(candidates_v) < theta && return nothing
    v_ids = sort!(collect(candidates_v))
    u_ids = Vector{Int}(undef, 1 + length(candidates_u))
    u_ids[1] = ui
    for (t, u) in enumerate(candidates_u)
        u_ids[t + 1] = u
    end
    v_pos = Dict{Int,Int}(v => j for (j, v) in enumerate(v_ids))
    n_u = length(u_ids)
    n_v = length(v_ids)
    adj_u = [Int[] for _ in 1:n_u]
    adj_v = [Int[] for _ in 1:n_v]
    for (i, u) in enumerate(u_ids)
        for v in g.adjU[u]
            j = get(v_pos, v, 0)
            j == 0 && continue
            push!(adj_u[i], j)
            push!(adj_v[j], i)
        end
    end
    in_s_u = falses(n_u)
    in_c_u = falses(n_u)
    in_s_v = falses(n_v)
    in_c_v = trues(n_v)
    in_s_u[1] = true
    for i in 2:n_u
        in_c_u[i] = true
    end
    deg_s_v = zeros(Int, n_v)
    for j in adj_u[1]
        deg_s_v[j] += 1
    end
    return _KbpInst(
        u_ids, v_ids, adj_u, adj_v,
        in_s_u, in_c_u, in_s_v, in_c_v,
        trues(n_u), trues(n_v),
        [length(nbrs) for nbrs in adj_u],
        [length(nbrs) for nbrs in adj_v],
        zeros(Int, n_u), deg_s_v,
        Int[1], Int[],
        1, 0, n_u - 1, n_v,
        k, theta,
        zeros(Int, n_u), zeros(Int, n_v), 0,
    )
end

function _kbp_expand(g::BipartiteGraph, ui::Int, later::Set{Int}, k::Int, theta::Int, best_edges::Int)
    st = _kbp_local(g, ui, later, k, theta)
    st === nothing && return nothing
    while st.n_c_u > 0 || st.n_c_v > 0
        st.n_s_u + st.n_c_u < theta && return nothing
        st.n_s_v + st.n_c_v < theta && return nothing
        _kbp_reduce!(st)
        st.n_s_u + st.n_c_u < theta && return nothing
        st.n_s_v + st.n_c_v < theta && return nothing
        st.n_c_u == 0 && st.n_c_v == 0 && break
        picked = _kbp_pick(st)
        picked === nothing && break
        is_u, idx = picked
        if is_u
            _kbp_can_add_u(st, idx) ? _kbp_accept_u!(st, idx) : _kbp_discard_u!(st, idx)
        else
            _kbp_can_add_v(st, idx) ? _kbp_accept_v!(st, idx) : _kbp_discard_v!(st, idx)
        end
    end
    (st.n_s_u < theta || st.n_s_v < theta) && return nothing
    edges = _kbp_count_edges(st)
    edges > best_edges || return nothing
    left, right = _kbp_members(st)
    return left, right, edges
end

"""
    k_biplex(g, k, theta) -> NamedTuple

A large k-biplex of `g`. Every returned vertex has at most `k` non-neighbors
on the opposite side, and both `left` and `right` have at least `theta`
vertices. `edges` counts the edges induced between the two sides. `g` is not
mutated. Empty sides and `edges == 0` mean no such biplex was found.

The procedure is Algorithm 3 of Cui, Li, Dai, and Wang (arXiv:2607.07419).
Deletion peels a minimum-degree vertex off the larger side until the
remainder is a k-biplex; the right side is peeled when the two sides have
the same size, and degree ties keep the smaller id. Expansion then rebuilds
a biplex around each left vertex in non-decreasing degree order. The
candidates are the vertices within distance 3, excluding left vertices that
come earlier in that order. Inside that pool, a vertex with fewer than
`theta - k` neighbors is dropped, a vertex that would push a committed
vertex past `k` non-neighbors is dropped, and a vertex that is already safe
for every subset of the pool is kept at once. The vertex of largest induced
degree is tried next.

Any k-biplex whose two sides both have at least `2k + 1` vertices has
diameter at most 3, which is the reason the pool is a 3-hop ball. Use
`theta ≥ 2k + 1` for that guarantee.
"""
function k_biplex(g::BipartiteGraph, k::Int, theta::Int)
    k < 0 && throw(ArgumentError("k must be ≥ 0"))
    theta < 1 && throw(ArgumentError("theta must be ≥ 1"))
    best_left = Int[]
    best_right = Int[]
    best_edges = -1
    if length(g.adjU) >= theta && length(g.adjV) >= theta
        left, right, edges = _kbp_peel(g, k)
        if length(left) >= theta && length(right) >= theta
            best_left = left
            best_right = right
            best_edges = edges
        end
        order = collect(keys(g.adjU))
        sort!(order, by = u -> (length(g.adjU[u]), u))
        later = Set{Int}(order)
        for u in order
            delete!(later, u)
            grown = _kbp_expand(g, u, later, k, theta, best_edges)
            grown === nothing && continue
            best_left, best_right, best_edges = grown
        end
    end
    sort!(best_left)
    sort!(best_right)
    return (; left=best_left, right=best_right, edges=best_edges < 0 ? 0 : best_edges)
end
