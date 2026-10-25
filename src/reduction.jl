#=
Degree reduction.

One pass deletes every left-side vertex whose degree is strictly below
`theta`. Right-side vertices stay, including degree 1 and any that become
isolated once their left neighbors are gone. Left vertices only touch the
right side, so deleting them does not change any other left degree.

  include(joinpath(@__DIR__, "reduction.jl"))
  g, _ = load_bipartite_graph(ensure_dataset!("konect-small/euroroad"; data_root=DATA_ROOT))
  degree_reduce!(g, 3)

  julia --project=. src/reduction.jl
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

using Printf

const REDUCTION_DATASETS = ["konect-small/euroroad", "konect-small/moreno"]
const REDUCTION_THETAS = [1, 2, 3, 4, 5, 8]

"""Structural copy of a metadata-free bipartite graph (adjacency only)."""
function copy_graph(g::BipartiteGraph{Nothing})
    h = BipartiteGraph{Nothing}()
    for (u, nbrs) in g.adjU
        h.adjU[u] = copy(nbrs)
    end
    for (v, nbrs) in g.adjV
        h.adjV[v] = copy(nbrs)
    end
    return h
end

function _low_degree_left(g::BipartiteGraph, theta::Int)
    drop_u = Int[]
    for (u, nbrs) in g.adjU
        length(nbrs) < theta && push!(drop_u, u)
    end
    return drop_u
end

"""
    degree_reduce!(g, theta) -> NamedTuple

Delete every left-side vertex of degree `< theta` in one pass. Right-side
vertices are left in place. Mutates `g`. Returns the number removed and the
vertices still present (`|U| + |V|`).
"""
function degree_reduce!(g::BipartiteGraph, theta::Int)
    theta < 1 && return (; removed=0, remaining=length(g.adjU) + length(g.adjV))
    removed = 0
    for u in _low_degree_left(g, theta)
        rem_u!(g, u) && (removed += 1)
    end
    return (; removed, remaining=length(g.adjU) + length(g.adjV))
end

_median(xs) = (s = sort(xs); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)

"""
Run `degree_reduce!` on a fresh copy. `ms` is the median peel time over `reps`
runs (graph copy excluded). `pct` is the share of left-side vertices removed.
"""
function reduction_trial(g::BipartiteGraph{Nothing}, theta::Int; reps::Int=5)
    nU0 = length(g.adjU)
    nV = length(g.adjV)
    samples = Vector{Float64}(undef, reps)
    removed = remaining_u = 0
    for r in 1:reps
        work = copy_graph(g)
        t0 = time_ns()
        stats = degree_reduce!(work, theta)
        samples[r] = (time_ns() - t0) / 1e6
        removed = stats.removed
        remaining_u = length(work.adjU)
    end
    pct = nU0 == 0 ? 0.0 : 100 * removed / nU0
    return (; nU0, nV, removed, remaining_u, pct, ms=_median(samples))
end

function _load_mutable(key::AbstractString)
    path = ensure_dataset!(key; data_root=DATA_ROOT)
    g, nE = load_bipartite_graph(path)
    return g, nE
end

"""Print left-side removal percentage and peel time for each dataset and theta."""
function report_degree_reduction(;
    datasets::Vector{String}=REDUCTION_DATASETS,
    thetas::Vector{Int}=REDUCTION_THETAS,
    reps::Int=5,
)
    println("Degree reduction: one pass, drop left-side vertices with deg < theta.")
    println("Right-side vertices stay, including degree 1. Removed % is of |U|.")
    println("Time is the median of $reps peels after one warmup, in milliseconds. Copy and load are excluded.")
    println()

    # Compile before the first measured cell.
    tiny = BipartiteGraph{Nothing}()
    add_edge!(tiny, 1, 1, nothing)
    degree_reduce!(tiny, 1)

    for key in datasets
        g, nE = _load_mutable(key)
        nU, nV = length(g.adjU), length(g.adjV)
        println("$key   |U|=$nU  |V|=$nV  |E|=$nE")
        @printf("  %5s  %8s  %10s  %8s  %8s  %10s\n",
            "theta", "removed", "removed %", "|U| left", "|V|", "ms")
        for theta in thetas
            trial = reduction_trial(g, theta; reps)
            @printf("  %5d  %8d  %9.2f%%  %8d  %8d  %10.3f\n",
                theta, trial.removed, trial.pct, trial.remaining_u, trial.nV, trial.ms)
        end
        println()
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    report_degree_reduction()
end
