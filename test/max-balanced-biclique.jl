#=
Maximum balanced biclique (TSGR-MBBP).

  julia --project=. test/max-balanced-biclique.jl

The exact component search is checked against every balanced subset. The
full search is checked on the same graphs, including the case where tabu
search has to leave a size-1 dead end.
=#

include(joinpath(@__DIR__, "..", "src", "max-balanced-biclique.jl"))

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

function _brute_balanced(g)
    us = sort!(collect(keys(g.adjU)))
    vs = sort!(collect(keys(g.adjV)))
    best = 0
    nu, nv = length(us), length(vs)
    for mu in 0:((1 << nu) - 1)
        lu = [us[i] for i in 1:nu if mu & (1 << (i - 1)) != 0]
        k = length(lu)
        (k == 0 || k <= best || k > nv) && continue
        for mv in 0:((1 << nv) - 1)
            lv = Int[]
            for i in 1:nv
                mv & (1 << (i - 1)) == 0 && continue
                push!(lv, vs[i])
            end
            length(lv) == k || continue
            ok = true
            for u in lu
                nbrs = g.adjU[u]
                for v in lv
                    v in nbrs || (ok = false; break)
                end
                ok || break
            end
            ok && (best = k)
        end
    end
    return best
end

function _is_pair(g, left, right)
    isempty(left) && return isempty(right)
    for u in left
        haskey(g.adjU, u) || return false
        for v in right
            v in g.adjU[u] || return false
        end
    end
    return true
end

@testset "balanced predicate" begin
    g = _graph([(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)])
    @test is_balanced_biclique(g, [1, 2], [1, 2])
    @test !is_balanced_biclique(g, [1, 2, 3], [1, 2])
    @test !is_balanced_biclique(g, [1], [1, 2])
    @test !is_balanced_biclique(g, Int[], Int[])
    @test !is_balanced_biclique(g, [1, 1], [1, 2])
    @test_throws ArgumentError max_balanced_biclique(g; restarts=-1)
    @test_throws ArgumentError max_balanced_biclique(g; preset=:other)
    @test_throws ArgumentError max_balanced_biclique(g; depth=-1)
    @test_throws ArgumentError max_balanced_biclique(g; alpha=-0.1)
    @test_throws ArgumentError max_balanced_biclique(g; component_limit=-1)
    @test_throws ArgumentError max_balanced_biclique(g; exact_seconds=-1)
    @test_throws ArgumentError max_balanced_biclique(g; seconds=-1)
end

@testset "exact component search" begin
    rng = Xoshiro(4)
    for _ in 1:15
        edges = Tuple{Int,Int}[]
        for u in 1:4, v in 1:4
            rand(rng) < 0.5 && push!(edges, (u, v))
        end
        g = _graph(edges)
        optimum = _brute_balanced(g)
        us = sort!(collect(keys(g.adjU)))
        vs = sort!(collect(keys(g.adjV)))
        left, right, done = _mbb_exact(g.adjU, g.adjV, us, vs, 0, 5.0)
        @test done
        @test _is_pair(g, left, right)
        @test min(length(left), length(right)) == optimum
        if optimum > 0
            none, _, done_again = _mbb_exact(g.adjU, g.adjV, us, vs, optimum, 5.0)
            @test done_again
            @test isempty(none)
        end
    end
end

@testset "complete and lopsided graphs" begin
    full = _graph([(u, v) for u in 1:3 for v in 1:3])
    before = _snapshot(full)
    found = max_balanced_biclique(full; restarts=2, seed=1)
    @test found.left == [1, 2, 3]
    @test found.right == [1, 2, 3]
    @test found.balanced == 3
    @test found.proven
    @test is_balanced_biclique(full, found.left, found.right)
    @test _snapshot(full) == before

    # Four left vertices and two right vertices. The balanced size is 2.
    wide = _graph([(u, v) for u in 1:4 for v in 1:2])
    got = max_balanced_biclique(wide; restarts=2, seed=2)
    @test got.balanced == 2
    @test got.proven
    @test is_balanced_biclique(wide, got.left, got.right)
    @test length(got.left) == length(got.right)

    empty = BipartiteGraph{Nothing}()
    none = max_balanced_biclique(empty)
    @test none.balanced == 0
    @test none.proven
    @test none.left == Int[]
end

@testset "tabu search leaves the size-1 dead end" begin
    # Starting at left 3 and adding right 3 stops. The 3×2 biclique is the optimum.
    # Exact search is off, so only the restarts and the pushes can find it.
    g = _graph([
        (1, 1), (1, 2),
        (2, 1), (2, 2),
        (3, 1), (3, 2), (3, 3),
    ])
    optimum = _brute_balanced(g)
    @test optimum == 2
    found = max_balanced_biclique(g;
        restarts=25, depth=30, alpha=0.3, component_limit=0, seed=1)
    @test found.balanced == 2
    @test is_balanced_biclique(g, found.left, found.right)
end

@testset "reduction proves the optimum" begin
    rng = Xoshiro(11)
    for trial in 1:12
        edges = Tuple{Int,Int}[]
        for u in 1:4, v in 1:4
            rand(rng) < 0.45 && push!(edges, (u, v))
        end
        g = _graph(edges)
        optimum = _brute_balanced(g)
        found = max_balanced_biclique(g;
            restarts=8, depth=20, component_limit=20, exact_seconds=2, seed=trial)
        @test found.balanced == optimum
        optimum == 0 && continue
        @test is_balanced_biclique(g, found.left, found.right)
    end
end
