#=
Density and scaled density of an induced subgraph.

  julia --project=. test/metrics.jl
=#

include(joinpath(@__DIR__, "..", "src", "metrics", "density.jl"))
include(joinpath(@__DIR__, "..", "src", "metrics", "scaled-density.jl"))

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

@testset "density" begin
    g = _graph([(u, v) for u in 1:2 for v in 1:2])
    full = SubGraph(Set([1, 2]), Set([1, 2]))
    @test density(g, full) == 1.0
    @test scaled_density(g, full) == 2.0

    partial = _graph([(1, 1), (1, 2), (2, 1)])
    block = SubGraph(Set([1, 2]), Set([1, 2]))
    @test density(partial, block) == 0.75
    @test scaled_density(partial, block) == 1.5

    @test density(g, SubGraph(Set{Int}(), Set([1]))) == 0.0
    @test density(g, SubGraph(Set([1]), Set{Int}())) == 0.0
    @test scaled_density(g, SubGraph(Set{Int}(), Set([1]))) == 0.0
    @test scaled_density(g, SubGraph(Set([1]), Set{Int}())) == 0.0

    bigger = _graph([(1, 10), (1, 11), (2, 10)])
    before = _snapshot(bigger)
    one = SubGraph(Set([1]), Set([10]))
    @test density(bigger, one) == 1.0
    @test scaled_density(bigger, one) == 1.0
    @test density(bigger, SubGraph(Set([2]), Set([11]))) == 0.0
    @test scaled_density(bigger, SubGraph(Set([2]), Set([11]))) == 0.0
    both = SubGraph(Set([1, 2]), Set([10, 11]))
    @test density(bigger, both) == 0.75
    @test scaled_density(bigger, both) == 1.5
    @test _snapshot(bigger) == before

    fg = freeze(bigger)
    @test density(fg, one) == 1.0
    @test scaled_density(fg, both) == 1.5
end
