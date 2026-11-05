#=
Vertex δ-quasi-biclique heuristic.

  julia --project=. test/vertex-quasi.jl

The search is the greedy procedure from Liu, Li, and Wang, as outlined by
Ignatov, Ivanova, and Zamaletdinova (arXiv:2002.09880, Definition 7). It is
not an exact maximum. These checks pin the degree rule, the two-phase
search, and the choice of the larger orientation.
=#

include(joinpath(@__DIR__, "..", "src", "formulations", "vertex-quasi.jl"))

using Test

function _graph(edges)
    g = BipartiteGraph{Nothing}()
    for (u, v) in edges
        add_edge!(g, u, v, nothing)
    end
    return g
end

function _edge_count(g)
    n = 0
    for (_, nbrs) in g.adjU
        n += length(nbrs)
    end
    return n
end

"""No vertex outside the pair can be added without breaking the degree rule."""
function _locally_maximal(g, left, right, delta)
    L = Set(left)
    R = Set(right)
    for u in keys(g.adjU)
        u in L && continue
        is_vertex_quasi_biclique(g, vcat(left, u), right, delta) && return false
    end
    for v in keys(g.adjV)
        v in R && continue
        is_vertex_quasi_biclique(g, left, vcat(right, v), delta) && return false
    end
    return true
end

function _check(g, delta; tau=1, left, right)
    before = _edge_count(g)
    qb = vertex_quasi_biclique(g, delta; tau)
    @test qb.left == left
    @test qb.right == right
    @test is_vertex_quasi_biclique(g, qb.left, qb.right, delta)
    @test _locally_maximal(g, qb.left, qb.right, delta)
    @test _edge_count(g) == before
    return qb
end

@testset "degree rule" begin
    g = _graph([(1, 1), (1, 2), (2, 1)])
    @test is_vertex_quasi_biclique(g, [1, 2], [1, 2], 0.5)
    @test !is_vertex_quasi_biclique(g, [1, 2], [1, 2], 0)
    @test !is_vertex_quasi_biclique(g, Int[], [1], 0)
    @test !is_vertex_quasi_biclique(g, [1], Int[], 0)
    @test !is_vertex_quasi_biclique(g, [9], [1], 0)
    @test_throws ArgumentError is_vertex_quasi_biclique(g, [1], [1], 0.6)
    @test_throws ArgumentError vertex_quasi_biclique(g, -0.1)
    @test_throws ArgumentError vertex_quasi_biclique(g, 0; tau=0)
end

@testset "biclique and one missing edge" begin
    # The pendant edge is a smaller biclique. δ = 0 keeps the 2×2.
    bic = _graph([(1, 1), (1, 2), (2, 1), (2, 2), (3, 3)])
    _check(bic, 0; left=[1, 2], right=[1, 2])

    # One missing edge. δ = 0 drops a vertex; δ = 0.5 keeps everyone.
    hole = _graph([
        (1, 1), (1, 2), (1, 3),
        (2, 1), (2, 2), (2, 3),
        (3, 1), (3, 2),
    ])
    _check(hole, 0; left=[1, 2], right=[1, 2, 3])
    _check(hole, 0.5; left=[1, 2, 3], right=[1, 2, 3])

    single = _graph([(1, 10)])
    _check(single, 0; left=[1], right=[10])
end

@testset "orientations and tau" begin
    # The right-hand biclique is larger. Growing the left side first stops
    # on the smaller star; the swapped pass is the one that is returned.
    sided = _graph([
        (1, 1), (1, 2),
        (2, 10), (3, 10), (4, 10),
        (2, 11), (3, 11), (4, 11),
    ])
    _check(sided, 0; left=[2, 3, 4], right=[10, 11])

    # Forcing a second left vertex empties the opposite side of that pass.
    # The pair saved before that addition is the biclique on vertex 1.
    forced = _graph([(1, 10), (1, 11), (1, 12), (2, 20)])
    _check(forced, 0; tau=2, left=[1], right=[10, 11, 12])

    # Vertex 1 touches three private items and six shared ones. Vertices 2..9
    # touch only the shared items. Item 300 is a decoy star of ten private
    # users, so neither side's first vertex is the shared block. τ = 1 stops
    # on {1, 2} plus every item of 1. τ = 6 raises the degree rule enough to
    # drop the private items; expansion then takes the whole shared block.
    block_edges = Tuple{Int,Int}[]
    for v in 1:3
        push!(block_edges, (1, 100 + v))
    end
    for u in 1:9, v in 1:6
        push!(block_edges, (u, 200 + v))
    end
    for u in 10:19
        push!(block_edges, (u, 300))
    end
    block = _graph(block_edges)
    _check(block, 0.5; tau=1, left=[1, 2], right=vcat(101:103, 201:206))
    _check(block, 0.5; tau=6, left=collect(1:9), right=collect(201:206))
end

@testset "empty" begin
    bare = BipartiteGraph{Nothing}()
    add_u!(bare, 1)
    add_v!(bare, 2)
    qb = vertex_quasi_biclique(bare, 0.5)
    @test qb.left == Int[]
    @test qb.right == Int[]
end
