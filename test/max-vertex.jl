#=
Maximum vertex biclique against the definition and exhaustive search.

  julia --project=. test/max-vertex.jl

The search is exact. On every graph small enough to enumerate, the returned
pair is a biclique and its number of vertices equals the best common
neighborhood.
=#

include(joinpath(@__DIR__, "..", "src", "max-vertex.jl"))

using Random
using Test

function _graph(edges)
    g = BipartiteGraph{Nothing}()
    for (u, v) in edges
        add_edge!(g, u, v, nothing)
    end
    return g
end

function _snapshot(g)
    return (
        Dict(u => copy(nbrs) for (u, nbrs) in g.adjU),
        Dict(v => copy(nbrs) for (v, nbrs) in g.adjV),
    )
end

"""Largest |U| + |V|. The right side of a left set is its common neighborhood.
The empty left set takes every right vertex."""
function _brute_size(g)
    us = sort!(collect(keys(g.adjU)))
    vs = collect(keys(g.adjV))
    n = length(us)
    n > 12 && error("exhaustive check is only for graphs with ≤ 12 left vertices")
    best = 0
    for mask in 0:((1 << n) - 1)
        common = nothing
        k = 0
        for i in 1:n
            mask & (1 << (i - 1)) == 0 && continue
            k += 1
            nbrs = g.adjU[us[i]]
            if common === nothing
                common = copy(nbrs)
            else
                intersect!(common, nbrs)
            end
        end
        size = k == 0 ? length(vs) : k + length(common)
        size > best && (best = size)
    end
    return best
end

function _gaps(nbrs, nR)
    out = Int[]
    j = 1
    p = 1
    while true
        v, p = _mv_next_gap(nbrs, nR, j, p)
        v == 0 && break
        push!(out, v)
        j = v + 1
    end
    return out
end

function _check(g; left=nothing, right=nothing)
    before = _snapshot(g)
    found = max_vertex_biclique(g)
    @test _snapshot(g) == before
    @test issorted(found.left)
    @test issorted(found.right)
    @test is_vertex_biclique(g, found.left, found.right)
    @test length(found.left) + length(found.right) == _brute_size(g)
    if left !== nothing
        @test found.left == left
        @test found.right == right
    end
    again = max_vertex_biclique(g)
    @test again.left == found.left
    @test again.right == found.right
    return found
end

@testset "complement gaps" begin
    @test _gaps(Int[], 3) == [1, 2, 3]
    @test _gaps([1, 2, 3], 3) == Int[]
    @test _gaps([2], 4) == [1, 3, 4]
    @test _gaps([1, 3], 3) == [2]
    @test _gaps([2, 4], 5) == [1, 3, 5]
end

@testset "definition" begin
    g = _graph([(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)])
    @test is_vertex_biclique(g, [1, 2], [1, 2])
    @test is_vertex_biclique(g, Int[], Int[])
    @test is_vertex_biclique(g, Int[], [1, 2])
    @test is_vertex_biclique(g, [1], Int[])
    @test !is_vertex_biclique(g, [1, 2, 3], [1, 2])
    @test !is_vertex_biclique(g, [1, 1], [1])
    @test !is_vertex_biclique(g, [1], [1, 1])
    @test !is_vertex_biclique(g, [9], [1])
    @test !is_vertex_biclique(g, [1], [9])
end

@testset "empty and one side" begin
    bare = BipartiteGraph{Nothing}()
    _check(bare; left=Int[], right=Int[])

    left_only = BipartiteGraph{Nothing}()
    add_u!(left_only, 4)
    add_u!(left_only, 1)
    _check(left_only; left=[1, 4], right=Int[])

    right_only = BipartiteGraph{Nothing}()
    add_v!(right_only, 8)
    add_v!(right_only, 2)
    _check(right_only; left=Int[], right=[2, 8])
end

@testset "edgeless complement is complete" begin
    # |left| > |right|: the cover is the right side.
    wide = BipartiteGraph{Nothing}()
    for u in 1:3
        add_u!(wide, u)
    end
    for v in 1:2
        add_v!(wide, v)
    end
    _check(wide; left=[1, 2, 3], right=Int[])

    # |right| > |left|, and a tie: a perfect matching leaves no free left
    # vertex, so the biclique is the right side.
    tall = BipartiteGraph{Nothing}()
    for u in 1:2
        add_u!(tall, u)
    end
    for v in 10:12
        add_v!(tall, v)
    end
    _check(tall; left=Int[], right=[10, 11, 12])

    square = BipartiteGraph{Nothing}()
    for u in 1:2
        add_u!(square, u)
    end
    for v in 1:2
        add_v!(square, v)
    end
    _check(square; left=Int[], right=[1, 2])
end

@testset "complete and near-complete" begin
    complete = _graph([(u, v) for u in 1:2 for v in 10:12])
    _check(complete; left=[1, 2], right=[10, 11, 12])

    # One forbidden complement edge. Dropping either endpoint leaves five
    # vertices; König drops the left endpoint of that edge.
    hole = _graph([(u, v) for u in 1:2 for v in 1:2 if (u, v) != (2, 2)])
    _check(hole; left=[1], right=[1, 2])

    single = _graph([(0, -1)])
    _check(single; left=[0], right=[-1])
end

@testset "unique proper biclique" begin
    # Five left vertices share one right vertex. The other right vertices are
    # isolated, so the only optimum is the star: size 6, against 5 and 4.
    star = _graph([(u, 1) for u in 1:5])
    add_v!(star, 2)
    add_v!(star, 3)
    add_v!(star, 4)
    _check(star; left=[1, 2, 3, 4, 5], right=[1])

    # Left 1 meets every right vertex. Left 2 meets only the first.
    universal = _graph([(1, 10), (1, 11), (1, 12), (2, 10)])
    _check(universal; left=[1], right=[10, 11, 12])

    # K_{3,2} plus an isolated left vertex. The isolated vertex is the only
    # complement edge's left side, and it is the vertex König deletes.
    block = _graph([(u, v) for u in 1:3 for v in 1:2])
    add_u!(block, 4)
    _check(block; left=[1, 2, 3], right=[1, 2])
end

@testset "several optima" begin
    # Each edge, either full side, all have two vertices. The size is what
    # the optimum fixes; König picks one of them.
    disjoint = _graph([(1, 1), (2, 2)])
    found = _check(disjoint)
    @test length(found.left) + length(found.right) == 2

    cycle = _graph([(1, 1), (2, 1), (2, 2), (3, 2), (3, 3), (1, 3)])
    _check(cycle)
end

@testset "matches exhaustive search" begin
    rng = MersenneTwister(1)
    for _ in 1:40
        nL = rand(rng, 0:8)
        nR = rand(rng, 0:8)
        p = rand(rng, (0.0, 0.25, 0.5, 1.0))
        edges = Tuple{Int,Int}[]
        for u in 1:nL, v in 1:nR
            rand(rng) < p && push!(edges, (u + 3, -(v + 1)))
        end
        g = _graph(edges)
        # Keep an isolated vertex when the draw produced none of that side
        # but the range still asked for it.
        if nL > 0 && isempty(g.adjU)
            add_u!(g, 4)
        end
        if nR > 0 && isempty(g.adjV)
            add_v!(g, -2)
        end
        _check(g)
    end
end

@testset "edgeless shortcut agrees with König" begin
    for (nL, nR) in ((1, 1), (3, 2), (2, 3), (4, 4), (5, 1))
        adj = [Int[] for _ in 1:nL]
        pairU, pairV = _mv_matching(adj, nL, nR)
        inL, inR = _mv_konig(adj, pairU, pairV, nL, nR)
        g = BipartiteGraph{Nothing}()
        for u in 1:nL
            add_u!(g, u)
        end
        for v in 1:nR
            add_v!(g, v)
        end
        found = max_vertex_biclique(g)
        @test found.left == [u for u in 1:nL if inL[u]]
        @test found.right == [v for v in 1:nR if inR[v]]
    end
end
