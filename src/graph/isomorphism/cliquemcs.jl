#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    findmcis, mcissize,
    findmces, mcessize



"""
    findmcis(G::UndirectedGraph, H::UndirectedGraph; kwargs...
        ) -> Tuple{Dict{Int,Int},Symbol}

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
    mapping = Dict{Int,Int}()
    status = :invalidinput
    (nodecount(G) == 0 || nodecount(H) == 0) && return (mapping, status)
    # Generate modular product
    if haskey(kwargs, :topological) && kwargs[:topological]
        eflt = tpconstraintfilter(G, H, edgematcher; kwargs...)
    else
        eflt = modprodedgefilter(G, H, edgematcher)
    end
    prod = modularproduct(
        G, H, nodematcher=nodematcher, edgefilter=eflt; kwargs...)
    # Clique detection
    if haskey(kwargs, :connected) && kwargs[:connected]
        (cliques, status) = maximalconncliques(prod; kwargs...)
    else
        (cliques, status) = maximalcliques(prod; kwargs...)
    end
    maxclique = sortstablemax(collect(cliques), by=length, init=[])
    mapping = Dict(
        nodeattr(prod, n).g => nodeattr(prod, n).h for n in maxclique)
    return (mapping, status)
end



"""
    findmces(G::UndirectedGraph, H::UndirectedGraph; kwargs...
        ) -> Tuple{Dict{Int,Int},Symbol}

Compute maximum common edge induced subgraph between `G` and `H`.
"""
function findmces(G::UndirectedGraph, H::UndirectedGraph;
        nodematcher=(g,h)->true, edgematcher=(g,h)->true, kwargs...)
    mapping = Dict{Int,Int}()
    status = :invalidinput
    (edgecount(G) == 0 || edgecount(H) == 0) && return (mapping, status)
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
    if haskey(kwargs, :connected) && kwargs[:connected]
        (cliques, status) = maximalconncliques(prod; kwargs...)
    else
        (cliques, status) = maximalcliques(prod; kwargs...)
    end
    for edges in cliques
        length(edges) > length(mapping) || continue
        mp = Dict(nodeattr(prod, e).g => nodeattr(prod, e).h for e in edges)
        delta_y_correction!(mp, G, H)
        if length(mp) > length(mapping)
            mapping = mp
        end
    end
    return (mapping, status)
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
mcissize(G, H; kwargs...) = length(findmcis(G, H; kwargs...)[1])


"""
    mcessize(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Int

Return the maximum common edge induced subgraph size (number of edges).
"""
mcessize(G, H; kwargs...) = length(findmces(G, H; kwargs...)[1])
