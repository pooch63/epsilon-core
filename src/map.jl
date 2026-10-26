#=
Right-side similarity graph.

  include(joinpath(@__DIR__, "map.jl"))
  include(joinpath(@__DIR__, "sp.jl"))

  h = map_graph(g, 2, shortest_path_metric(g, 2))
  h = map_graph(fg, 0.5, (r1, r2) -> abs(x[r1] - x[r2]))

  julia --project=. src/map.jl konect-small/euroroad --theta=2 --epsilon=2

`h` is an undirected graph on exactly the right-side vertices of `g`. Two
distinct vertices are adjacent when `metric(r1, r2) ≤ epsilon`, and the edge
weight is that distance. `metric` is any callable of two right-side ids; it is
invoked once per unordered pair, with the smaller id first. `g` may be a
`BipartiteGraph` or a `FrozenBipartite`.

Run as a script, this loads the dataset you name, degree-reduces the left
side (`reduction.jl`), builds the shortest-path similarity graph, and prints
how many nodes and edges remain. The full pipeline — left reduction, one
mapping, similarity-degree pruning, left reduction again, then saving both
graphs — is `bin/process.jl`.
=#

isdefined(@__MODULE__, :__CORE_JL__) || include(joinpath(@__DIR__, "core.jl"))

function _map_require_epsilon(epsilon::Real)
    if !(epsilon >= 0 && isfinite(epsilon))
        throw(ArgumentError("epsilon must be finite and ≥ 0"))
    end
    return epsilon
end

_right_ids(g::BipartiteGraph) = sort!(collect(keys(g.adjV)))
_right_ids(g::FrozenBipartite) = sort(g.v_ids)

"""
    map_graph(g, epsilon, metric) -> UnipartiteGraph{Float64}

Unipartite graph on every right-side vertex of `g`. Distinct vertices `r1`
and `r2` are joined when `metric(r1, r2) ≤ epsilon`. The edge weight is that
distance. `metric` is called once per unordered pair, smaller id first.
"""
function map_graph(g::Union{BipartiteGraph,FrozenBipartite}, epsilon::Real, metric)
    _map_require_epsilon(epsilon)
    h = UnipartiteGraph{Float64}()
    rights = _right_ids(g)
    for r in rights
        add_vertex!(h, r)
    end
    n = length(rights)
    for i in 1:n
        r1 = rights[i]
        for j in (i + 1):n
            r2 = rights[j]
            d = metric(r1, r2)
            d <= epsilon && add_edge!(h, r1, r2, convert(Float64, d))
        end
    end
    return h
end

# Top-level so the included methods exist before `map_dataset` calls them.
isdefined(@__MODULE__, :degree_reduce!) || include(joinpath(@__DIR__, "reduction.jl"))
isdefined(@__MODULE__, :shortest_path_metric) || include(joinpath(@__DIR__, "sp.jl"))

"""
    similarity_reduce!(h, g, theta) -> NamedTuple

Delete every vertex of the similarity graph `h` whose degree is `< theta - 1`.
Those vertices are the right side of `g`. A sweep reads degrees as they stand
at the start, then deletes; the next sweep sees the edges that deletion
removed. Repeats until every remaining vertex has degree ≥ `theta - 1`. Each
deleted vertex is removed from `h` and from the right side of `g`. The
distance metric is not run again.
"""
function similarity_reduce!(h::UnipartiteGraph, g::BipartiteGraph, theta::Int)
    theta < 1 && return (; removed=0, sweeps=0, remaining=num_vertices(h))
    removed = 0
    sweeps = 0
    while true
        drop = Int[]
        for (v, nbrs) in h.adj
            length(nbrs) < theta - 1 && push!(drop, v)
        end
        isempty(drop) && break
        sweeps += 1
        for v in drop
            rem_vertex!(h, v) && (removed += 1)
            rem_v!(g, v)
        end
    end
    return (; removed, sweeps, remaining=num_vertices(h))
end

"""
    sync_map!(h, g) -> Int

Remove similarity-graph vertices that are no longer on the right side of `g`.
Their edges are deleted with them. Distances are not recomputed.
"""
function sync_map!(h::UnipartiteGraph, g::BipartiteGraph)
    removed = 0
    for v in collect(keys(h.adj))
        haskey(g.adjV, v) && continue
        rem_vertex!(h, v) && (removed += 1)
    end
    return removed
end

"""
    reduce_then_map(g, theta, epsilon) -> NamedTuple

Left-side degree reduction, then one shortest-path mapping of the right side,
then similarity-degree pruning, then the left-side reduction again. Later
pruning deletes vertices and the edges incident to them. The mapping is not
rebuilt.
"""
function reduce_then_map(g::BipartiteGraph, theta::Int, epsilon::Real)
    theta < 1 && throw(ArgumentError("theta must be ≥ 1"))
    _map_require_epsilon(epsilon)
    first = degree_reduce!(g, theta)
    h = map_graph(g, epsilon, shortest_path_metric(g, epsilon))
    peeled = similarity_reduce!(h, g, theta)
    second = degree_reduce!(g, theta)
    dropped = sync_map!(h, g)
    return (; h, first, peeled, second, dropped)
end

_bipartite_edges(g::BipartiteGraph) = sum(length(nbrs) for nbrs in values(g.adjU); init=0)

"""Directory `data/<dataset>/reduced/theta-N-epsilon-X/`."""
function reduced_output_dir(key::AbstractString, theta::Int, epsilon::Real; data_root::AbstractString=DATA_ROOT)
    eps = epsilon == floor(epsilon) ? string(Int(epsilon)) : string(epsilon)
    return joinpath(dataset_dir(key; data_root), "reduced", "theta-$theta-epsilon-$eps")
end

"""Write `u,v` edges. Isolated vertices are omitted; degree pruning has already dropped them when `theta ≥ 1`."""
function save_bipartite_graph(path::AbstractString, g::BipartiteGraph)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "u,v")
        for u in sort!(collect(keys(g.adjU)))
            for v in sort!(collect(g.adjU[u]))
                println(io, u, ',', v)
            end
        end
    end
    return path
end

"""Write each undirected edge once, as `u,v,w` with `u < v`."""
function save_unipartite_graph(path::AbstractString, h::UnipartiteGraph)
    mkpath(dirname(path))
    edges = Tuple{Int,Int,Float64}[]
    for (u, nbrs) in h.adj
        for (v, w) in nbrs
            u < v && push!(edges, (u, v, Float64(w)))
        end
    end
    sort!(edges)
    open(path, "w") do io
        println(io, "u,v,w")
        for (u, v, w) in edges
            println(io, u, ',', v, ',', w)
        end
    end
    return path
end

function _map_flag(arg::AbstractString)
    value = strip(split(arg, '=', limit=2)[2])
    isempty(value) && throw(ArgumentError("$arg is missing a value"))
    return value
end

function _print_map_usage()
    println(stderr, "Usage:")
    println(stderr, "  julia --project=. src/map.jl <provider>/<dataset> --theta=N --epsilon=X")
    println(stderr, "")
    println(stderr, "Options:")
    println(stderr, "  --theta=N         drop left-side vertices of degree < N before mapping")
    println(stderr, "  --epsilon=X       connect right-side vertices at shortest-path distance ≤ X")
    println(stderr, "")
    println(stderr, "Example:")
    println(stderr, "  julia --project=. src/map.jl konect-small/euroroad --theta=2 --epsilon=2")
end

"""
    map_dataset(key, theta, epsilon)

Load `key`, delete left-side vertices of degree `< theta`, then build the
right-side shortest-path graph at `epsilon`.
"""
function map_dataset(key::AbstractString, theta::Int, epsilon::Real)
    theta < 1 && throw(ArgumentError("theta must be ≥ 1"))
    path = ensure_dataset!(key; data_root=DATA_ROOT)
    g, _ = load_bipartite_graph(path)
    stats = degree_reduce!(g, theta)
    nE = sum(length(nbrs) for nbrs in values(g.adjU); init=0)
    h = map_graph(g, epsilon, shortest_path_metric(g, epsilon))
    return (;
        key, theta, epsilon,
        removed=stats.removed,
        bipartite_nodes=stats.remaining, bipartite_edges=nE,
        nodes=num_vertices(h), edges=num_edges(h),
    )
end

function _parse_map_args(args)
    dataset = nothing
    theta = nothing
    epsilon = nothing
    for arg in args
        if startswith(arg, "--theta=")
            theta = parse(Int, _map_flag(arg))
        elseif startswith(arg, "--epsilon=")
            epsilon = parse(Float64, _map_flag(arg))
        elseif arg in ("--theta", "--epsilon")
            throw(ArgumentError("$arg requires a value, as in --theta=2"))
        elseif startswith(arg, "--")
            throw(ArgumentError("Unknown flag: $arg"))
        elseif dataset === nothing
            dataset = arg
        else
            throw(ArgumentError("Unexpected argument: $arg"))
        end
    end
    if dataset === nothing || theta === nothing || epsilon === nothing
        throw(ArgumentError("Pass <provider>/<dataset> --theta=N --epsilon=X"))
    end
    theta < 1 && throw(ArgumentError("theta must be ≥ 1"))
    return (; dataset, theta, epsilon)
end

function _print_map_report(result)
    eps = result.epsilon == floor(result.epsilon) ? string(Int(result.epsilon)) : string(result.epsilon)
    println(result.key)
    println("  degree reduction  theta=$(result.theta)  removed $(result.removed) left vertices")
    println("  bipartite         $(result.bipartite_nodes) nodes, $(result.bipartite_edges) edges")
    println("  map               $(result.nodes) nodes, $(result.edges) edges  (shortest path ≤ $eps)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if any(a -> a in ("--help", "-h"), ARGS)
        _print_map_usage()
    else
        cfg = try
            _parse_map_args(ARGS)
        catch e
            e isa ArgumentError || rethrow()
            println(stderr, e.msg)
            println(stderr)
            _print_map_usage()
            exit(1)
        end
        _print_map_report(map_dataset(cfg.dataset, cfg.theta, cfg.epsilon))
    end
end
