#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    findmcis, mcissize,
    findmces, mcessize


mutable struct MCSResult
    mapping::Dict{Int,Int}
    isvalidinput::Bool
    istimedout::Bool
    istargetreached::Bool
end


"""
    findmcis(G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> MCSResult

Compute maximum common induced subgraph between `G` and `H`.


## Keyword arguments:

- connectivity: `:Connected`, `:Disconnected`  (for mode=:MCS)
- constraint: `TopologicalConstraint`, `DiameterRestriction` (for mode=:MCS)
- timeout(Int): return suboptimal result so far if the execution time exceeded
the given value (in second).
- threshold(Int): return suboptimal result so far if the given mcs size
achieved.
- nodematcher(Function): node matcher function that takes two node indices as
arguments.
- edgematcher(Function): edge matcher function that takes two edge indices as
arguments.
- theta(Int):
    distance mismatch tolerance in topologically constrainted MCS
- diameter(Int):
    diameter size in MCS with diameter restriction
"""
function findmcis(G::UndirectedGraph, H::UndirectedGraph; kwargs...)
    res = MCSResult(Dict{Int,Int}(), false, false, false)
    (nodecount(G) == 0 || nodecount(H) == 0) && return res
    res.isvalidinput = true
    # Generate modular product
    if haskey(kwargs, :constraint)
        # TODO: constraints
        eflt = topologicalconstraint(G, H; kwargs...)
    end
    prod = modularproduct(G, H; kwargs...)
    # Clique detection
    state = FindCliqueState{ModularProduct}(prod; kwargs...)
    for nodes in state.channel
        if length(nodes) > length(res.mapping)
            res.mapping = Dict(
                nodeattr(prod, n).g => nodeattr(prod, n).h for n in nodes)
        end
    end
    res.istimedout = state.status == :timedout
    res.istargetreached = state.status == :targetreached
    return res
end


"""
    findmces(G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> MCSResult

Compute maximum common edge induced subgraph between `G` and `H`.
"""
function findmces(G::UndirectedGraph, H::UndirectedGraph; kwargs...)
    res = MCSResult(Dict{Int,Int}(), false, false, false)
    (edgecount(G) == 0 || edgecount(H) == 0) && return res
    res.isvalidinput = true
    lg = linegraph(G)
    lh = linegraph(H)
    nmatch = (g, h) -> true
    ematch = (g, h) -> true
    if haskey(kwargs, :nodematcher)
        ematch = lgedgematcher(lg, lh, kwargs[:nodematcher])
        if haskey(kwargs, :edgematcher)
            nmatch = lgnodematcher(
                lg, lh, kwargs[:nodematcher], kwargs[:edgematcher])
        end
    end
    # Generate modular product
    if haskey(kwargs, :constraint)
        # TODO: constraints
        eflt = topologicalconstraint(G, H; kwargs...)
    end
    eflt = modprodedgefilter(lg, lh, ematch)
    prod = modularproduct(lg, lh, nodematcher=nmatch, edgefilter=eflt)
    # Clique detection
    state = FindCliqueState{ModularProduct}(prod; kwargs...)
    for edges in state.channel
        length(edges) > length(res.mapping) || continue
        mp = Dict(nodeattr(prod, e).g => nodeattr(prod, e).h for e in edges)
        delta_y_correction!(mp, G, H)
        if length(mp) > length(res.mapping)
            res.mapping = mp
        end
    end
    res.istimedout = state.status == :timedout
    res.istargetreached = state.status == :targetreached
    return res
end


function modprodedgefilter(G, H, edgematcher)
    return function (g1, g2, h1, h2)
        if hasedge(G, g1, g2) != hasedge(H, h1, h2)
            return false
        elseif !hasedge(G, g1, g2) && !hasedge(H, h1, h2)
            return true
        else
            return edgematcher(findedgekey(G, g1, g2), findedgekey(H, h1, h2))
        end
    end
end


"""
    mcissize(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Int

Return the maximum common induced subgraph size (number of nodes).
"""
mcissize(G, H; kwargs...) = length(findmcis(G, H; kwargs...).mapping)


"""
    mcessize(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Int

Return the maximum common edge induced subgraph size (number of edges).
"""
mcessize(G, H; kwargs...) = length(findmces(G, H; kwargs...).mapping)
