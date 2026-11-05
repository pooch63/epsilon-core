#=
Exact maximum vertex biclique: a biclique maximizing |left| + |right|.

  include(joinpath(@__DIR__, "max-vertex.jl"))

  found = max_vertex_biclique(g)
  found.left
  found.right

  julia --project=. test/max-vertex.jl

A biclique is a pair of sets in which every left vertex is adjacent to every
right vertex. Either side may be empty; that pair is an independent set of
the bipartite complement and is optimal when no biclique is larger. The
search does not mutate g.

Yannakakis observed that this is polynomial. The bipartite complement has an
edge exactly where g does not, and a biclique of g is an independent set of
that complement. A maximum independent set is the complement of a minimum
vertex cover, and König's theorem builds the cover from a maximum matching.
Hopcroft–Karp finds the matching. The complement is not stored: its edges are
the gaps in the sorted original neighborhoods.
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

"""
    is_vertex_biclique(g, left, right) -> Bool

`left` and `right` are duplicate-free vertices of their own side of `g`, and
every left vertex is adjacent to every right vertex. Either side may be
empty: there is then no pair to check.
"""
function is_vertex_biclique(g::BipartiteGraph, left, right)
    left_ids = collect(Int, left)
    right_ids = collect(Int, right)
    length(left_ids) == length(Set{Int}(left_ids)) || return false
    length(right_ids) == length(Set{Int}(right_ids)) || return false
    for u in left_ids
        haskey(g.adjU, u) || return false
    end
    for v in right_ids
        haskey(g.adjV, v) || return false
    end
    for u in left_ids
        nbrs = g.adjU[u]
        for v in right_ids
            v in nbrs || return false
        end
    end
    return true
end

"""Sorted original neighborhoods as dense right indices. Rows are strictly
increasing, so a complement neighbor is a gap in that row."""
function _mv_original(g::BipartiteGraph, u_ids::Vector{Int}, v_ids::Vector{Int})
    v_index = Dict{Int,Int}(v => j for (j, v) in enumerate(v_ids))
    adj = Vector{Vector{Int}}(undef, length(u_ids))
    n_edges = 0
    for (i, u) in enumerate(u_ids)
        row = Int[]
        sizehint!(row, length(g.adjU[u]))
        for v in g.adjU[u]
            j = get(v_index, v, 0)
            j == 0 && continue
            push!(row, j)
        end
        sort!(row)
        unique!(row)
        adj[i] = row
        n_edges += length(row)
    end
    return adj, n_edges
end

"""Next complement neighbor at or after dense index `j`. `p` points at the
first original neighbor that might still meet `j`. Returns `(0, p)` past
the last right vertex."""
@inline function _mv_next_gap(nbrs::Vector{Int}, nR::Int, j::Int, p::Int)
    nn = length(nbrs)
    while j <= nR
        if p <= nn && nbrs[p] < j
            p += 1
            continue
        end
        if p <= nn && nbrs[p] == j
            j += 1
            p += 1
            continue
        end
        return j, p
    end
    return 0, p
end

function _mv_bfs!(adj, pairU, pairV, dist, nil_dist::Ref{Int}, nL::Int, nR::Int)
    inf = nL + 1
    nil_dist[] = inf
    queue = Int[]
    sizehint!(queue, nL)
    for u in 1:nL
        if pairU[u] == 0
            dist[u] = 0
            push!(queue, u)
        else
            dist[u] = inf
        end
    end
    head = 1
    while head <= length(queue)
        u = queue[head]
        head += 1
        # Deeper than the shortest augmenting path cannot improve this phase.
        dist[u] >= nil_dist[] && continue
        nbrs = adj[u]
        j = 1
        p = 1
        while true
            v, p = _mv_next_gap(nbrs, nR, j, p)
            v == 0 && break
            j = v + 1
            w = pairV[v]
            if w == 0
                nil_dist[] == inf && (nil_dist[] = dist[u] + 1)
            elseif dist[w] == inf
                dist[w] = dist[u] + 1
                push!(queue, w)
            end
        end
    end
    return nil_dist[] != inf
end

function _mv_dfs!(adj, pairU, pairV, dist, nil_dist::Ref{Int}, u::Int, nR::Int)
    inf = length(dist) + 1
    nbrs = adj[u]
    j = 1
    p = 1
    while true
        v, p = _mv_next_gap(nbrs, nR, j, p)
        v == 0 && break
        j = v + 1
        w = pairV[v]
        if w == 0
            nil_dist[] == dist[u] + 1 || continue
            pairV[v] = u
            pairU[u] = v
            return true
        end
        dist[w] == dist[u] + 1 || continue
        _mv_dfs!(adj, pairU, pairV, dist, nil_dist, w, nR) || continue
        pairV[v] = u
        pairU[u] = v
        return true
    end
    dist[u] = inf
    return false
end

function _mv_matching(adj, nL::Int, nR::Int)
    pairU = zeros(Int, nL)
    pairV = zeros(Int, nR)
    dist = zeros(Int, nL)
    nil_dist = Ref(0)
    while _mv_bfs!(adj, pairU, pairV, dist, nil_dist, nL, nR)
        for u in 1:nL
            pairU[u] == 0 && _mv_dfs!(adj, pairU, pairV, dist, nil_dist, u, nR)
        end
    end
    return pairU, pairV
end

"""
Independent set from König's construction. Unmatched left vertices of the
complement start an alternating search: non-matching edges run left to right,
matching edges run right to left. The reached left vertices and the unreached
right vertices are a maximum independent set.
"""
function _mv_konig(adj, pairU, pairV, nL::Int, nR::Int)
    seenL = falses(nL)
    seenR = falses(nR)
    queue = Int[]
    for u in 1:nL
        pairU[u] == 0 || continue
        seenL[u] = true
        push!(queue, u)
    end
    head = 1
    while head <= length(queue)
        u = queue[head]
        head += 1
        matched = pairU[u]
        nbrs = adj[u]
        j = 1
        p = 1
        while true
            v, p = _mv_next_gap(nbrs, nR, j, p)
            v == 0 && break
            j = v + 1
            v == matched && continue
            seenR[v] && continue
            seenR[v] = true
            w = pairV[v]
            if w != 0 && !seenL[w]
                seenL[w] = true
                push!(queue, w)
            end
        end
    end
    inR = falses(nR)
    for v in 1:nR
        seenR[v] || (inR[v] = true)
    end
    n_match = 0
    for u in 1:nL
        pairU[u] != 0 && (n_match += 1)
    end
    n_mis = count(seenL) + count(inR)
    n_mis == nL + nR - n_match ||
        error("König construction disagrees with the maximum matching")
    return seenL, inR
end

function _mv_keep(ids::Vector{Int}, keep::BitVector)
    out = Int[]
    for i in eachindex(ids)
        keep[i] && push!(out, ids[i])
    end
    return out
end

"""
    max_vertex_biclique(g) -> NamedTuple

A biclique `(left, right)` of `g` with the largest value of
`|left| + |right|`.

Every returned left vertex is adjacent to every returned right vertex. Either
side may be empty. `left` and `right` are sorted. `g` is not mutated.

The complement of `g` has an edge wherever `g` does not. Hopcroft–Karp finds
a maximum matching of that complement, and König's theorem turns the matching
into a minimum vertex cover. The vertices outside the cover are the biclique.
Each phase scans every missing edge, so a graph that is neither complete nor
edgeless takes `O(√(|L| + |R|) · |L| · |R|)` time. An edgeless graph skips the
matching: the complement is complete, the biclique is the larger side, and a
tie returns the right side. A complete bipartite graph is already the unique
maximum and is returned whole.
"""
function max_vertex_biclique(g::BipartiteGraph)
    u_ids = sort!(collect(keys(g.adjU)))
    v_ids = sort!(collect(keys(g.adjV)))
    nL = length(u_ids)
    nR = length(v_ids)
    if nL == 0 || nR == 0
        return (; left=u_ids, right=v_ids)
    end
    adj, n_edges = _mv_original(g, u_ids, v_ids)
    if n_edges == 0
        # The complement is complete. When |left| ≤ |right| a perfect matching
        # covers every left vertex, the cover is the left side, and the
        # biclique is the right side. When the left side is larger, the
        # biclique is the left side.
        if nL > nR
            return (; left=u_ids, right=Int[])
        end
        return (; left=Int[], right=v_ids)
    end
    complete = true
    for row in adj
        if length(row) != nR
            complete = false
            break
        end
    end
    if complete
        return (; left=u_ids, right=v_ids)
    end
    pairU, pairV = _mv_matching(adj, nL, nR)
    inL, inR = _mv_konig(adj, pairU, pairV, nL, nR)
    return (; left=_mv_keep(u_ids, inL), right=_mv_keep(v_ids, inR))
end
