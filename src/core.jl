#=
=================================================================================
epsilon-core entry point
=================================================================================
Include this from any script to get the shared GraphLogic toolkit:

  include(joinpath(@__DIR__, "..", "src", "core.jl"))
  fg = load_dataset("konect-small/euroroad")   # downloads on first use

GraphLogic provides BipartiteGraph / FrozenBipartite / SubGraph, the Subgraph
helpers, load_bipartite_graph, dataset discovery, and the download providers.
See ../graph-logic/README.md for the full API.
=#

const __CORE_JL__ = true

if Base.find_package("GraphLogic") === nothing
    error("GraphLogic is not installed. From the epsilon-core root run:\n" *
          "  julia --project=. -e 'using Pkg; Pkg.instantiate()'")
end

using GraphLogic

const ROOT = dirname(@__DIR__)
const DATA_ROOT = joinpath(ROOT, "data")

# CLI args with `--data-root=data/` filled in, so bin/ scripts work from any cwd.
with_data_root(args) =
    any(startswith("--data-root="), args) ? args : [args; "--data-root=$DATA_ROOT"]
