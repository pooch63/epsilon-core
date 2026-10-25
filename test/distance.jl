#=
Diffusion distance against the random-walk definition.

  julia --project=. test/distance.jl

The spectral map has to match Σ_k (P^t[i,k] − P^t[j,k])^2 / π_k, with P the
lazy walk on the whole graph or the two-step walk on the right side.
=#

include(joinpath(@__DIR__, "..", "src", "distance.jl"))

using LinearAlgebra
using Test

function _weighted(edges)
    g = BipartiteGraph{Float64}()
    for (u, v, w) in edges
        add_edge!(g, u, v, w)
    end
    return g
end

function _unweighted(edges)
    g = BipartiteGraph{Nothing}()
    for (u, v) in edges
        add_edge!(g, u, v, nothing)
    end
    return g
end

"""Direct Σ (P^t[i,k] − P^t[j,k])^2 / π_k, built without the spectral code."""
function _reference(edges, time::Int, walk::Symbol, r1::Int, r2::Int)
    r1 == r2 && return 0.0
    u_ids = sort!(unique(e[1] for e in edges))
    v_ids = sort!(unique(e[2] for e in edges))
    u_of = Dict(id => i for (i, id) in enumerate(u_ids))
    v_of = Dict(id => i for (i, id) in enumerate(v_ids))
    (haskey(v_of, r1) && haskey(v_of, r2)) || return Inf
    nU, nV = length(u_ids), length(v_ids)
    B = zeros(nU, nV)
    for (u, v, w) in edges
        B[u_of[u], v_of[v]] = w
    end
    du = vec(sum(B; dims=2))
    dv = vec(sum(B; dims=1))
    if walk === :step
        P = Diagonal(1 ./ dv) * (B' * (Diagonal(1 ./ du) * B))
        Pt = P^time
        π = dv ./ sum(dv)
        i, j = v_of[r1], v_of[r2]
        acc = 0.0
        for k in 1:nV
            acc += (Pt[i, k] - Pt[j, k])^2 / π[k]
        end
        return sqrt(acc)
    end
    n = nU + nV
    degree = [du; dv]
    A = zeros(n, n)
    A[1:nU, nU + 1:n] = B
    A[nU + 1:n, 1:nU] = B'
    P = (I(n) + A ./ degree) / 2
    Pt = P^time
    π = degree ./ sum(degree)
    i = nU + v_of[r1]
    j = nU + v_of[r2]
    acc = 0.0
    for k in 1:n
        acc += (Pt[i, k] - Pt[j, k])^2 / π[k]
    end
    return sqrt(acc)
end

function _agree(edges, walk::Symbol)
    g = _weighted(edges)
    fg = freeze(g)
    rights = sort!(unique(e[2] for e in edges))
    for time in 0:4
        for rank in (nothing, 1, 8)
            metric = diffusion_metric(g; time, rank, walk)
            for a in rights, b in rights
                got = metric(a, b)
                frozen = diffusion_distance(fg, a, b; time, rank, walk)
                @test frozen ≈ got atol=1e-8 rtol=1e-8
                rank === nothing || continue
                want = _reference(edges, time, walk, a, b)
                @test got ≈ want atol=1e-8 rtol=1e-8
            end
        end
        full = diffusion_metric(g; time, walk)
        short = diffusion_metric(g; time, rank=1, walk)
        for a in rights, b in rights
            @test short(a, b) <= full(a, b) + 1e-8
        end
    end
    return nothing
end

@testset "matches the walk" begin
    graphs = [
        [(1, 10, 1.0), (1, 11, 1.0), (2, 10, 1.0)],
        [(1, 10, 1.0), (1, 11, 1.0), (2, 10, 1.0), (2, 11, 1.0)],
        [(1, 10, 1.0), (1, 11, 1.0), (2, 11, 1.0), (2, 12, 1.0)],
        [(1, 10, 2.0), (1, 11, 0.5), (2, 11, 3.0), (3, 10, 1.0)],
        [(1, 1, 1.0), (1, 2, 2.0), (2, 1, 0.5)],
        [(1, 10, 1.0), (1, 11, 1.0), (2, 10, 1.0), (3, 12, 1.0), (4, 12, 2.0), (4, 13, 1.0)],
    ]
    for edges in graphs, walk in (:lazy, :step)
        _agree(edges, walk)
    end
end

@testset "closed forms" begin
    g = _unweighted([(1, 10), (1, 11), (2, 10), (2, 11)])
    # Lazy mass is still on the start vertex, so the two items differ by √2.
    @test diffusion_distance(g, 10, 11; time=1) ≈ sqrt(2) atol=1e-8
    # One two-step lands on the same distribution from either item.
    @test diffusion_distance(g, 10, 11; time=1, walk=:step) == 0.0
    @test diffusion_distance(g, 10, 11; time=0, walk=:step) > 0
    @test diffusion_distance(g, 10, 10; time=3) == 0.0
    @test diffusion_distance(freeze(g), 10, 11; time=1) ≈ diffusion_distance(g, 10, 11; time=1)

    weighted = _weighted([(1, 10, 2.0), (1, 11, 2.0), (2, 10, 0.5), (2, 11, 0.5)])
    @test diffusion_distance(weighted, 10, 11; time=3, walk=:step) == 0.0
    @test diffusion_distance(weighted, 10, 11; time=1, walk=:lazy) > 0
end

@testset "missing vertices" begin
    g = _unweighted([(1, 10)])
    add_v!(g, 99)
    @test diffusion_distance(g, 99, 99) == 0.0
    @test diffusion_distance(g, 10, 99) == Inf
    @test diffusion_distance(g, 10, 1000) == Inf
    @test diffusion_metric(g)(10, 10) == 0.0

    empty = BipartiteGraph{Nothing}()
    @test diffusion_distance(empty, 1, 1) == 0.0
    @test diffusion_distance(empty, 1, 2) == Inf
end

@testset "rejected diffusion arguments" begin
    g = _unweighted([(1, 10)])
    @test_throws ArgumentError diffusion_metric(g; time=-1)
    @test_throws ArgumentError diffusion_metric(g; rank=0)
    @test_throws ArgumentError diffusion_metric(g; walk=:both)
    bad = BipartiteGraph{Float64}()
    add_edge!(bad, 1, 10, 0.0)
    @test_throws ArgumentError diffusion_metric(bad)
    add_edge!(bad, 1, 11, -2.0)
    @test_throws ArgumentError diffusion_distance(bad, 10, 11)
end
