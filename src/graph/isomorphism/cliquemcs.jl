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

- connected(Bool): if true, apply connected MCS constraint.
- topological(Bool): if true, apply topological constraint.
- diameter(Int): distance cutoff for topological constraint.
- tolerance(Int): distance mismatch tolerance for topological constraint.
- timeout(Int): abort calculation and return suboptimal results so far if the
execution time has reached the given value (default=60, in seconds).
- targetsize(Int): abort calculation and return suboptimal result so far if the
given mcs size achieved.
- nodematcher(Function): node matcher function that takes two node indices as
arguments.
- edgematcher(Function): edge matcher function that takes two edge indices as
arguments.
"""
function findmcis(G::UndirectedGraph, H::UndirectedGraph;
        nodematcher=(g,h)->true, edgematcher=(g,h)->true, kwargs...)
    res = MCSResult(Dict{Int,Int}(), false, false, false)
    (nodecount(G) == 0 || nodecount(H) == 0) && return res
    res.isvalidinput = true
    # Generate modular product
    if haskey(kwargs, :topological) && kwargs[:topological]
        eflt = tpconstraintfilter(G, H, edgematcher; kwargs...)
    else
        eflt = modprodedgefilter(G, H, edgematcher)
    end
    prod = modularproduct(
        G, H, nodematcher=nodematcher, edgefilter=eflt; kwargs...)
    # Clique detection
    state = FindCliqueState{ModularProduct}(prod; kwargs...)
    if haskey(kwargs, :connected) && kwargs[:connected]
        expand!(state, nodeset(prod), nodeset(prod), nodeset(prod))
    else
        expand!(state, nodeset(prod), nodeset(prod))
    end
    for nodes in state.cliques
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
function findmces(G::UndirectedGraph, H::UndirectedGraph;
        nodematcher=(g,h)->true, edgematcher=(g,h)->true, kwargs...)
    res = MCSResult(Dict{Int,Int}(), false, false, false)
    (edgecount(G) == 0 || edgecount(H) == 0) && return res
    res.isvalidinput = true
    lg = linegraph(G)
    lh = linegraph(H)
    ematch = lgedgematcher(lg, lh, nodematcher)
    nmatch = lgnodematcher(lg, lh, nodematcher, edgematcher)
    # Generate modular product
    if haskey(kwargs, :topological) && kwargs[:topological]
        eflt = tpconstraintfilter(lg, lh, ematch; kwargs...)
    else
        eflt = modprodedgefilter(lg, lh, ematch)
    end
    prod = modularproduct(lg, lh, nodematcher=nmatch, edgefilter=eflt)
    # Clique detection
    state = FindCliqueState{ModularProduct}(prod; kwargs...)
    if haskey(kwargs, :connected) && kwargs[:connected]
        expand!(state, nodeset(prod), nodeset(prod), nodeset(prod))
    else
        expand!(state, nodeset(prod), nodeset(prod))
    end
    for edges in state.cliques
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
        hasedge(G, g1, g2) == hasedge(H, h1, h2) || return false
        !hasedge(G, g1, g2) && !hasedge(H, h1, h2) && return true
        return edgematcher(findedgekey(G, g1, g2), findedgekey(H, h1, h2))
    end
end


function tpconstraintfilter(
        G, H, edgematcher; diameter=typemax(Int), tolerance=0, kwargs...)
    gdist = Dict(n => bfsdepth(adjacencies, G, n) for n in nodeset(G))
    hdist = Dict(n => bfsdepth(adjacencies, H, n) for n in nodeset(H))
    return function (g1, g2, h1, h2)
        gdist[g1][g2] > diameter && return false
        hdist[h1][h2] > diameter && return false
        abs(gdist[g1][g2] - hdist[h1][h2]) > tolerance && return false
        !hasedge(G, g1, g2) && !hasedge(H, h1, h2) && return true
        return edgematcher(findedgekey(G, g1, g2), findedgekey(H, h1, h2))
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
