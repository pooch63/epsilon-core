#=
Heuristic against the overlap coefficient, plus a few checks a custom
distance can stage more easily than overlap can. The diffusion test uses
the two-step walk, where identical neighborhoods are distance 0.

  julia --project=. test/heuristic.jl

Overlap distance is `1 - |N(u) ∩ N(v)| / min(|N(u)|, |N(v)|)`. It is 0 when
one neighborhood contains the other. On each tiny graph the search is also
compared with an exhaustive enumeration of right-hand subsets, which is the
largest core under the same distance.
=#

include(joinpath(@__DIR__, "..", "src", "heuristic.jl"))
include(joinpath(@__DIR__, "..", "src", "distance.jl"))

using Test

function _graph(edges)
    g = BipartiteGraph{Nothing}()
    for (u, v) in edges
        add_edge!(g, u, v, nothing)
    end
    return g
end

"""Size of the largest core, by trying every subset of the right side."""
function _optimum(g::BipartiteGraph, theta::Int, epsilon::Real, metric)
    rights = sort!(collect(keys(g.adjV)))
    n = length(rights)
    n > 16 && error("exhaustive check is only for graphs with ≤ 16 right vertices")
    best = 0
    for mask in 0:((1 << n) - 1)
        group = Int[]
        for i in 1:n
            mask & (1 << (i - 1)) == 0 && continue
            push!(group, rights[i])
        end
        ok = true
        for a in 1:length(group), b in (a + 1):length(group)
            if _pair_distance(metric, group[a], group[b]) > epsilon
                ok = false
                break
            end
        end
        ok || continue
        covered = 0
        wanted = Set(group)
        for (_, nbrs) in g.adjU
            hits = 0
            for v in nbrs
                v in wanted && (hits += 1)
            end
            hits >= theta && (covered += 1)
        end
        covered > best && (best = covered)
    end
    return best
end

function _check(g, theta, epsilon, metric; left, right, pruned)
    core = cohesive_core(g, theta, epsilon, metric)
    @test core.left == left
    @test core.right == right
    @test core.pruned == pruned
    @test is_cohesive_core(g, core.left, core.right, theta, epsilon, metric)
    @test length(core.left) == _optimum(g, theta, epsilon, metric)
    return core
end

@testset "peel" begin
    g = _graph([(1, 1), (2, 10), (2, 11), (3, 10), (3, 11)])
    stats = peel_core!(g, 2)
    @test stats.removed_left == 1
    @test stats.removed_right == 1
    @test stats.sweeps == 1
    @test !haskey(g.adjU, 1)
    @test !haskey(g.adjV, 1)
    @test sort!(collect(keys(g.adjU))) == [2, 3]
    @test sort!(collect(keys(g.adjV))) == [10, 11]

    # Removing an isolated right vertex does not drop any left degree, so a
    # second pass has nothing to do. A graph already above theta is untouched.
    clean = _graph([(1, 1), (1, 2), (2, 1), (2, 2)])
    again = peel_core!(clean, 2)
    @test again.removed_left == 0
    @test again.removed_right == 0
    @test again.sweeps == 0
end

@testset "overlap coefficient" begin
    # Private items sit inside the popular item's neighborhood, so overlap
    # distance is 0, but they share no buyers with each other (distance 1).
    # The later biclique is the larger core. The popular item is tried first
    # because it touches the most left vertices.
    mixed = _graph([
        (1, 1), (2, 1), (3, 1), (4, 1), (5, 1),
        (1, 11), (2, 12), (3, 13), (4, 14), (5, 15),
        (20, 30), (21, 30), (22, 30),
        (20, 31), (21, 31), (22, 31),
        (20, 32), (21, 32), (22, 32),
    ])
    mixed_d = overlap_metric(mixed)
    @test overlap_distance(mixed, 1, 11) == 0.0
    @test overlap_distance(mixed, 11, 12) == 1.0
    @test overlap_distance(mixed, 30, 31) == 0.0
    @test overlap_distance(mixed, 1, 30) == 1.0
    nU = length(mixed.adjU)
    _check(mixed, 2, 0.0, mixed_d;
        left=[20, 21, 22], right=[30, 31], pruned=7)
    @test length(mixed.adjU) == nU

    # Buyers below theta still count toward overlap. After they are peeled,
    # items 1 and 2 would look identical; on the original graph they are 0.5
    # apart, so ε = 0 admits no core and ε = 0.5 admits user 1.
    bridge = _graph([(1, 1), (1, 2), (2, 2), (3, 1)])
    @test overlap_distance(bridge, 1, 2) == 0.5
    bridge_d = overlap_metric(bridge)
    _check(bridge, 2, 0.0, bridge_d; left=Int[], right=Int[], pruned=2)
    _check(bridge, 2, 0.5, bridge_d; left=[1], right=[1, 2], pruned=1)

    # Each pair of items shares one of two buyers. Overlap distance is 0.5,
    # and the three users form a core that is not a biclique.
    triangle = _graph([(1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (3, 1)])
    @test overlap_distance(triangle, 1, 2) == 0.5
    tri_d = overlap_metric(triangle)
    _check(triangle, 2, 0.5, tri_d; left=[1, 2, 3], right=[1, 2, 3], pruned=2)
    @test !is_cohesive_core(triangle, [1, 2, 3], [1, 2, 3], 2, 0.0, tri_d)
    @test !is_cohesive_core(triangle, [1], [1, 2, 3], 2, 0.5, tri_d)

    # A user of degree 1 cannot qualify and must not appear in the core.
    with_stray = _graph([
        (1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (3, 1),
        (4, 4),
    ])
    _check(with_stray, 2, 0.5, overlap_metric(with_stray);
        left=[1, 2, 3], right=[1, 2, 3], pruned=2)
end

@testset "closeness weight" begin
    # θ = 3. After the seed and item 2, users 1 and 2 are one neighbor from
    # the threshold and users 3 and 4 are two away. Item 4 and item 3 each
    # help two users, but item 4 helps the closer pair. Item 3 is not similar
    # to item 4, so taking the closer item is what produces the core.
    # Item 5 keeps users 3 and 4 from being peeled and is similar to nobody.
    g = _graph([
        (1, 1), (2, 1), (3, 1), (4, 1),
        (1, 2), (2, 2),
        (3, 3), (4, 3),
        (1, 4), (2, 4),
        (3, 5), (4, 5),
    ])
    pairs = Set([(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)])
    metric(a::Int, b::Int) = a == b ? 0.0 : ((min(a, b), max(a, b)) in pairs ? 0.0 : 1.0)
    _check(g, 3, 0.0, metric; left=[1, 2], right=[1, 2, 4], pruned=4)
end

@testset "diffusion walk" begin
    # Same buyers on 10 and 11, so one step of the two-step walk cannot
    # tell them apart. Item 12 is not in that group. User 4 is below theta.
    g = _graph([
        (1, 10), (2, 10), (3, 10),
        (1, 11), (2, 11), (3, 11),
        (1, 12), (4, 12),
    ])
    step = diffusion_metric(g; time=1, walk=:step)
    @test step(10, 11) == 0.0
    @test step(10, 12) > 0
    _check(g, 2, 0.0, step; left=[1, 2, 3], right=[10, 11], pruned=2)
end

@testset "rejected arguments" begin
    g = _graph([(1, 1)])
    metric(_, _) = 0.0
    @test_throws ArgumentError cohesive_core(g, 0, 0.0, metric)
    @test_throws ArgumentError cohesive_core(g, 1, -0.1, metric)
    @test_throws ArgumentError peel_core!(g, 0)
end
