#=
The k-biplex heuristic against the definition, a direct peel, and exhaustive
search on tiny graphs.

  julia --project=. test/k-biplex.jl
=#

include(joinpath(@__DIR__, "..", "src", "k-biplex.jl"))

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

function _induced_edges(g, left, right)
    right_set = Set{Int}(right)
    edges = 0
    for u in left
        for v in g.adjU[u]
            v in right_set && (edges += 1)
        end
    end
    return edges
end

"""Same deletion rule as `_kbp_peel`, with degrees recomputed from scratch."""
function _slow_peel(g, k)
    left = Set(keys(g.adjU))
    right = Set(keys(g.adjV))
    function degree(id, is_left)
        nbrs = is_left ? g.adjU[id] : g.adjV[id]
        other = is_left ? right : left
        hits = 0
        for x in nbrs
            x in other && (hits += 1)
        end
        return hits
    end
    function biplex()
        for u in left
            length(right) - degree(u, true) > k && return false
        end
        for v in right
            length(left) - degree(v, false) > k && return false
        end
        return true
    end
    while !isempty(left) && !isempty(right) && !biplex()
        if length(left) > length(right)
            pick, best = nothing, typemax(Int)
            for u in left
                degree_u = degree(u, true)
                if degree_u < best || (degree_u == best && u < pick)
                    pick, best = u, degree_u
                end
            end
            delete!(left, pick)
        else
            pick, best = nothing, typemax(Int)
            for v in right
                degree_v = degree(v, false)
                if degree_v < best || (degree_v == best && v < pick)
                    pick, best = v, degree_v
                end
            end
            delete!(right, pick)
        end
    end
    return left, right
end

function _optimum(g, k, theta)
    left_ids = sort!(collect(keys(g.adjU)))
    right_ids = sort!(collect(keys(g.adjV)))
    n_left = length(left_ids)
    n_right = length(right_ids)
    (n_left > 12 || n_right > 12) && error("exhaustive check is only for graphs with ≤ 12 vertices on a side")
    best = 0
    for left_mask in 0:((1 << n_left) - 1)
        count_ones(left_mask) >= theta || continue
        left = Int[]
        for i in 1:n_left
            left_mask & (1 << (i - 1)) == 0 && continue
            push!(left, left_ids[i])
        end
        for right_mask in 0:((1 << n_right) - 1)
            count_ones(right_mask) >= theta || continue
            right = Int[]
            for j in 1:n_right
                right_mask & (1 << (j - 1)) == 0 && continue
                push!(right, right_ids[j])
            end
            is_k_biplex(g, left, right, k) || continue
            edges = _induced_edges(g, left, right)
            edges > best && (best = edges)
        end
    end
    return best
end

function _check_result(g, found, k, theta)
    @test found.edges == _induced_edges(g, found.left, found.right)
    @test issorted(found.left)
    @test issorted(found.right)
    if found.edges == 0
        @test isempty(found.left)
        @test isempty(found.right)
    else
        @test is_k_biplex(g, found.left, found.right, k)
        @test length(found.left) >= theta
        @test length(found.right) >= theta
    end
    return found
end

@testset "heap and peel" begin
    heap = NTuple{3,Int}[]
    items = NTuple{3,Int}[(d, id, 1) for id in 1:30 for d in (id * 5 % 11,)]
    for item in items
        _kbp_hpush!(heap, item)
    end
    previous = nothing
    for _ in 1:length(items)
        item = _kbp_hpop!(heap)
        previous === nothing || @test _kbp_before(previous, item)
        previous = item
    end
    @test isempty(heap)

    rng = MersenneTwister(1)
    for _ in 1:20
        edges = Tuple{Int,Int}[]
        for u in 1:6, v in 1:6
            rand(rng) < 0.45 && push!(edges, (u, v))
        end
        g = _graph(edges)
        for k in 0:2
            peeled_left, peeled_right, _ = _kbp_peel(g, k)
            slow_left, slow_right = _slow_peel(g, k)
            @test Set(peeled_left) == slow_left
            @test Set(peeled_right) == slow_right
        end
    end
end

@testset "definition" begin
    g = _graph([(1, 1), (1, 2), (2, 1), (2, 2), (3, 1)])
    @test is_k_biplex(g, [1, 2], [1, 2], 0)
    @test !is_k_biplex(g, [1, 2, 3], [1, 2], 0)
    @test is_k_biplex(g, [1, 2, 3], [1, 2], 1)
    @test !is_k_biplex(g, [1, 1], [1], 1)
    @test !is_k_biplex(g, [9], [1], 1)
    @test is_k_biplex(g, Int[], Int[], 0)
end

@testset "already a k-biplex" begin
    clique = _graph([(u, v) for u in 1:4 for v in 10:13])
    before = _snapshot(clique)
    found = k_biplex(clique, 1, 3)
    @test found.left == [1, 2, 3, 4]
    @test found.right == [10, 11, 12, 13]
    @test found.edges == 16
    _check_result(clique, found, 1, 3)
    @test _snapshot(clique) == before

    # One missing edge: each of its endpoints has a single non-neighbor.
    almost = _graph([(u, v) for u in 1:3 for v in 1:3 if (u, v) != (3, 3)])
    found = k_biplex(almost, 1, 3)
    @test found.left == [1, 2, 3]
    @test found.right == [1, 2, 3]
    @test found.edges == 8
    _check_result(almost, found, 1, 3)
end

@testset "deletion keeps the larger biclique" begin
    # Two disjoint bicliques. Minimum-degree deletion eats the K_{2,2} and
    # leaves the K_{3,3}.
    edges = [(u, v) for u in 1:2 for v in 1:2]
    append!(edges, [(u, v) for u in 3:5 for v in 3:5])
    g = _graph(edges)
    found = k_biplex(g, 0, 2)
    @test found.left == [3, 4, 5]
    @test found.right == [3, 4, 5]
    @test found.edges == 9
    _check_result(g, found, 0, 2)
    @test found.edges == _optimum(g, 0, 2)
end

@testset "expansion restores a biclique deletion cut into" begin
    # K_{4,2} plus a pendant right vertex on left vertex 4. The sides start
    # unbalanced, so deletion removes left vertex 1 before the pendant and
    # finishes on a 6-edge biclique. The 3-hop expansion puts the deleted
    # vertex back and drops the pendant, which is the 8-edge biclique.
    edges = [(u, v) for u in 1:4 for v in 1:2]
    push!(edges, (4, 3))
    g = _graph(edges)
    found = k_biplex(g, 0, 2)
    @test found.left == [1, 2, 3, 4]
    @test found.right == [1, 2]
    @test found.edges == 8
    _check_result(g, found, 0, 2)
    @test found.edges == _optimum(g, 0, 2)

    # The same graph has no 1-biplex with three vertices on both sides:
    # the pendant right vertex misses three left vertices.
    none = k_biplex(g, 1, 3)
    @test none.left == Int[]
    @test none.right == Int[]
    @test none.edges == 0
    @test _optimum(g, 1, 3) == 0
end

@testset "no feasible biplex" begin
    g = _graph([(1, 1), (1, 2), (2, 2)])
    found = k_biplex(g, 1, 3)
    @test found.left == Int[]
    @test found.right == Int[]
    @test found.edges == 0
    empty = BipartiteGraph{Nothing}()
    @test k_biplex(empty, 1, 1).edges == 0
end

@testset "within the exhaustive optimum" begin
    graphs = [
        _graph([(u, v) for u in 1:3 for v in 1:3]),
        _graph([(1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 3), (3, 2)]),
        _graph([(1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (3, 1), (4, 1), (4, 4)]),
        _graph([(u, v) for u in 1:4 for v in 1:4 if u != v]),
    ]
    rng = MersenneTwister(7)
    for _ in 1:12
        edges = Tuple{Int,Int}[]
        for u in 1:4, v in 1:4
            rand(rng) < 0.55 && push!(edges, (u, v))
        end
        push!(graphs, _graph(edges))
    end
    for g in graphs
        for k in 0:2
            theta = max(1, 2k - 1)
            found = k_biplex(g, k, theta)
            _check_result(g, found, k, theta)
            @test found.edges <= _optimum(g, k, theta)
            again = k_biplex(g, k, theta)
            @test again.left == found.left
            @test again.right == found.right
            @test again.edges == found.edges
        end
    end
end

@testset "rejected arguments" begin
    g = _graph([(1, 1)])
    @test_throws ArgumentError k_biplex(g, -1, 1)
    @test_throws ArgumentError k_biplex(g, 1, 0)
    @test !is_k_biplex(g, [1], [1], -1)
end
