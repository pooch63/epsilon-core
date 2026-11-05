#=
Maximum-edge biclique (BSC) and the core score applied to a biclique.

  julia --project=. test/max-edge-biclique.jl

The sampler is not an exact maximum-edge solver. These checks pin the
expansion step, the both-sides seed, and the raw / repaired core numbers.
=#

include(joinpath(@__DIR__, "..", "src", "max-edge-biclique.jl"))

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

function _brute_edges(g)
    us = sort!(collect(keys(g.adjU)))
    vs = sort!(collect(keys(g.adjV)))
    best = 0
    nu, nv = length(us), length(vs)
    for mu in 0:((1 << nu) - 1), mv in 0:((1 << nv) - 1)
        lu = [us[i] for i in 1:nu if mu & (1 << (i - 1)) != 0]
        lv = [vs[i] for i in 1:nv if mv & (1 << (i - 1)) != 0]
        (isempty(lu) || isempty(lv)) && continue
        edges = length(lu) * length(lv)
        edges <= best && continue
        ok = true
        for u in lu
            nbrs = g.adjU[u]
            for v in lv
                v in nbrs || (ok = false; break)
            end
            ok || break
        end
        ok && (best = edges)
    end
    return best
end

"""Largest core under the strict diameter test, by every right subset."""
function _brute_core(g, theta, epsilon, metric)
    rights = sort!(collect(keys(g.adjV)))
    n = length(rights)
    best = 0
    for mask in 0:((1 << n) - 1)
        group = Int[]
        for i in 1:n
            mask & (1 << (i - 1)) == 0 && continue
            push!(group, rights[i])
        end
        length(group) >= theta || continue
        ok = true
        if length(group) == 1
            ok = epsilon > 0
        else
            for a in 1:length(group), b in (a + 1):length(group)
                d = group[a] <= group[b] ? metric(group[a], group[b]) : metric(group[b], group[a])
                d < epsilon || (ok = false; break)
            end
        end
        ok || continue
        wanted = Set(group)
        covered = 0
        for (_, nbrs) in g.adjU
            hits = 0
            for v in nbrs
                v in wanted && (hits += 1)
            end
            hits >= theta && (covered += 1)
        end
        best = max(best, covered)
    end
    return best
end

@testset "biclique predicate" begin
    g = _graph([(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)])
    @test is_biclique(g, [1, 2], [1, 2])
    @test is_biclique(g, [2, 1], [2, 1])
    @test !is_biclique(g, [1, 2, 3], [1, 2])
    @test !is_biclique(g, Int[], [1])
    @test !is_biclique(g, [1], Int[])
    @test !is_biclique(g, [1, 1], [1])
    @test !is_biclique(g, [9], [1])
end

@testset "expansion finds the planted edge biclique" begin
    # Shaham, Yu, and Li, Figure 1. ({3,4,6}, {3,5}) has 6 edges.
    paper = _graph([
        (1, 2),
        (3, 2), (3, 3), (3, 5),
        (4, 3), (4, 5), (4, 7),
        (6, 3), (6, 5),
        (8, 5), (8, 7),
    ])
    before = _snapshot(paper)
    found = max_edge_biclique(paper; samples=40, seed_size=1, both_sides=false, seed=1)
    @test found.left == [3, 4, 6]
    @test found.right == [3, 5]
    @test found.edges == 6
    @test is_biclique(paper, found.left, found.right)
    @test _snapshot(paper) == before

    full = _graph([(u, v) for u in 1:3 for v in 1:3])
    got = max_edge_biclique(full; samples=1, seed_size=1, seed=1)
    @test got.edges == 9
    @test got.left == [1, 2, 3]
    @test got.right == [1, 2, 3]
end

@testset "a private neighbor needs a larger seed or the other side" begin
    # {1,2} × {1,2} has 4 edges. Each left vertex also has a private right vertex,
    # so a seed of one left vertex never rebuilds the 2×2.
    g = _graph([
        (1, 1), (1, 2), (1, 3),
        (2, 1), (2, 2), (2, 4),
    ])
    one = max_edge_biclique(g; samples=20, seed_size=1, both_sides=false, seed=1)
    @test one.edges == 3
    pair = max_edge_biclique(g; samples=1, seed_size=2, both_sides=false, seed=1)
    @test pair.left == [1, 2]
    @test pair.right == [1, 2]
    @test pair.edges == 4
    other = max_edge_biclique(g; samples=30, seed_size=1, both_sides=true, seed=1)
    @test other.edges == 4
    @test is_biclique(g, other.left, other.right)
end

@testset "samples cover the optimum on tiny graphs" begin
    rng = Xoshiro(7)
    for _ in 1:12
        edges = Tuple{Int,Int}[]
        for u in 1:4, v in 1:4
            rand(rng) < 0.45 && push!(edges, (u, v))
        end
        g = _graph(edges)
        optimum = _brute_edges(g)
        best = 0
        # A discriminating seed can be any size up to the biclique. Cover each size.
        for p in 1:4
            found = max_edge_biclique(g; samples=40, seed_size=p, both_sides=true, seed=p)
            best = max(best, found.edges)
            found.edges == 0 && continue
            @test is_biclique(g, found.left, found.right)
        end
        @test best == optimum
    end
end

@testset "empty graph and bad arguments" begin
    g = BipartiteGraph{Nothing}()
    found = max_edge_biclique(g)
    @test found.left == Int[]
    @test found.right == Int[]
    @test found.edges == 0
    bare = _graph([(1, 1)])
    @test_throws ArgumentError max_edge_biclique(bare; beta=1.2)
    @test_throws ArgumentError max_edge_biclique(bare; samples=0)
    @test_throws ArgumentError max_edge_biclique(bare; seed_size=0)
    @test_throws ArgumentError max_edge_biclique(bare; max_samples=0)
end

@testset "raw and repaired core numbers" begin
    # Right vertices 1,2,3 are pairwise within 2. 10 and 11 are a second pair.
    # Distance is the absolute difference of the ids.
    edges = [(u, v) for u in 1:2 for v in (1, 2, 3, 10, 11)]
    g = _graph(edges)
    metric = (a, b) -> abs(a - b)
    split = score_biclique(g, [1, 2], [1, 2, 3, 10, 11], 2, 2.5, metric)
    @test !split.raw_ok
    @test split.raw_size == 0
    @test split.repaired_ok
    @test split.repaired_right == [1, 2, 3]
    @test split.repaired_left == [1, 2]
    @test split.repaired_size == 2
    @test is_biclique(g, split.repaired_left, split.repaired_right)
    @test _brute_core(g, 2, 2.5, metric) >= split.repaired_size

    whole = score_biclique(g, [1, 2], [1, 2, 3], 2, 2.5, metric)
    @test whole.raw_ok
    @test whole.raw_size == 2
    @test whole.repaired_right == [1, 2, 3]
    @test whole.repaired_size == 2
    @test _brute_core(g, 2, 2.5, metric) >= whole.raw_size

    # Distance exactly ε is outside the strict diameter test. θ = 2 then
    # rejects the singletons that remain; θ = 1 keeps one of them.
    touch = score_biclique(g, [1, 2], [1, 3], 2, 2, metric)
    @test !touch.raw_ok
    @test !touch.repaired_ok
    alone = score_biclique(g, [1, 2], [1, 3], 1, 2, metric)
    @test !alone.raw_ok
    @test alone.repaired_ok
    @test length(alone.repaired_right) == 1
    @test alone.repaired_size == 2

    lone = _graph([(1, 1)])
    @test !score_biclique(lone, [1], [1], 1, 0, metric).raw_ok
    one = score_biclique(lone, [1], [1], 1, 0.1, metric)
    @test one.raw_ok
    @test one.raw_size == 1

    # Five left vertices and one right vertex cannot meet θ = 2.
    wide = _graph([(u, 1) for u in 1:5])
    lopsided = score_biclique(wide, [1, 2, 3, 4, 5], [1], 2, 1, metric)
    @test !lopsided.raw_ok
    @test !lopsided.repaired_ok
    @test lopsided.raw_size == 0

    @test_throws ArgumentError score_biclique(lone, [1], [1], 0, 1, metric)
    @test_throws ArgumentError score_biclique(lone, [1], [1], 1, -0.1, metric)
    missing = score_biclique(lone, [1], [9], 1, 1, metric)
    @test !missing.raw_ok
    @test !missing.repaired_ok
end
