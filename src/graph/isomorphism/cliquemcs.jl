#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    maxcommonsubgraph, MaxCommonSubgraphResult


struct MaxCommonSubgraphResult
    mapping::Dict{Int,Int}
    status::Symbol
end

Base.size(result::MaxCommonSubgraphResult) = length(result.mapping)


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
        hasedge(G, g1, g2) && hasedge(H, h1, h2) || return false
        return edgematcher(findedgekey(G, g1, g2), findedgekey(H, h1, h2))
    end
end


"""
    maxcommonsubgraph(G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> MaxCommonSubgraphResult

Compute maximum common edge induced subgraph between `G` and `H`.
"""
function maxcommonsubgraph(G::UndirectedGraph, H::UndirectedGraph, matchtype;
        nodematcher=(g,h)->true, edgematcher=(g,h)->true, kwargs...)
    mapping = Dict{Int,Int}()
    status = :invalidinput
    if matchtype == :nodeinduced
        (nodecount(G) == 0 || nodecount(H) == 0) && return MaxCommonSubgraphResult(mapping, status)
        G_ = G
        H_ = H
        nmatch = nodematcher
        ematch = edgematcher
    else
        (edgecount(G) == 0 || edgecount(H) == 0) && return MaxCommonSubgraphResult(mapping, status)
        G_ = linegraph(G)
        H_ = linegraph(H)
        nmatch = lgnodematcher(G_, H_, nodematcher, edgematcher)
        ematch = lgedgematcher(G_, H_, nodematcher)
    end
    # Generate modular product
    if haskey(kwargs, :topological) && kwargs[:topological]
        eflt = tpconstraintfilter(G_, H_, ematch; kwargs...)
    else
        eflt = modprodedgefilter(G_, H_, ematch)
    end
    prod = modularproduct(G_, H_, nodematcher=nmatch, edgefilter=eflt)
    # Clique detection
    if haskey(kwargs, :connected) && kwargs[:connected]
        cqstate = maximalconncliques(prod; kwargs...)
    else
        cqstate = maximalcliques(prod; kwargs...)
    end
    if matchtype == :nodeinduced
        maxclique = maximumclique(cqstate)
        mapping = Dict(nodeattr(prod, n).g => nodeattr(prod, n).h for n in maxclique)
    else
        @assert length(mapping) == 0
        for edges in cqstate.cliques
            length(edges) > length(mapping) || continue
            mp = Dict(nodeattr(prod, e).g => nodeattr(prod, e).h for e in edges)
            delta_y_correction!(mp, G, H)
            if length(mp) > length(mapping)
                mapping = mp
            end
        end
    end
    return MaxCommonSubgraphResult(mapping, cqstate.status)
end
