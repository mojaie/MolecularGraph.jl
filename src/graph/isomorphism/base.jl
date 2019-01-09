#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    is_isomorphic,
    is_subgraph,
    is_edge_subgraph,
    isomorphism,
    edgeisomorphism,
    isomorphismiter,
    edgeisomorphismiter,
    maximumcommonsubgraph,
    delta_y_correction!,
    delta_y_mismatch,
    lgnodematcher,
    lgedgematcher,
    connectivity

"""
mode: Isomorphism(VF2), Subgraph(VF2), MCS(Clique)
subgraphtype: NodeInduced, EdgeInduced  # not for Isomorphism
algorithm: VF2, Clique
connectivity: Connected, Disconnected  # MCS only
constraint: TopologicalConstraint, DiameterRestriction  # MCS only

# parameters

timeout  # get suboptimal result
threshold  # get suboptimal result
nodematcher  # atom
edgematcher  # bond, Isomorph and Subgraph
connectivity # MCS only
mandatory  # VF2 only
forbidden  # VF2 only
theta  # Topological only
diameter  # Diameter only
"""


function is_isomorphic(G, H; kwargs...)
    return isomorphism(
        G, H; mode=:Isomorphism, algorithm=:VF2, kwargs...) !== nothing
end


function is_subgraph(G, H; kwargs...)
    """ True if G is an induced subgraph of H"""
    return isomorphism(
        H, G; mode=:Subgraph, subgraphtype=:NodeInduced, algorithm=:VF2,
        kwargs...) !== nothing
end


function is_edge_subgraph(G, H; kwargs...)
    """ True if G is an induced subgraph of H"""
    return isomorphism(
        H, G; mode=:Subgraph, subgraphtype=:EdgeInduced, algorithm=:VF2,
        kwargs...) !== nothing
end


function isomorphism(G, H;
        algorithm=:VF2, subgraphtype=:NodeInduced, kwargs...)
    if algorithm == :VF2
        if subgraphtype == :NodeInduced
            isomorphismvf2(G, H, subgraphtype=:NodeInduced; kwargs...)
        elseif subgraphtype == :EdgeInduced
            edgeisomorphismvf2(G, H, subgraphtype=:EdgeInduced; kwargs...)
        end
    elseif algorithm == :Clique
        """
        # TODO: not implemented yet
        if subgraphtype == :NodeInduced
            isomorphismclique(G, H; kwargs...)
        elseif subgraphtype == :EdgeInduced
            edgeisomorphismclique(G, H; kwargs...)
        end
        """
    end
end

edgeisomorphism(
    G, H; kwargs...) = isomorphism(G, H, subgraphtype=:EdgeInduced; kwargs...)


function isomorphismiter(G, H;
        algorithm=:VF2, subgraphtype=:NodeInduced, kwargs...)
    if algorithm == :VF2
        if subgraphtype == :NodeInduced
            isomorphismitervf2(G, H, subgraphtype=:NodeInduced; kwargs...)
        elseif subgraphtype == :EdgeInduced
            edgeisomorphismitervf2(G, H, subgraphtype=:EdgeInduced; kwargs...)
        end
    elseif algorithm == :Clique
        # TODO: not implemented yet
        """
        if subgraphtype == :NodeInduced
            isomorphismiterclique(G, H; kwargs...)
        elseif subgraphtype == :EdgeInduced
            edgeisomorphismiterclique(G, H; kwargs...)
        end
        """
    end
end

edgeisomorphismiter(G, H; kwargs...) = isomorphismiter(
    G, H, subgraphtype=:EdgeInduced; kwargs...)


function maximumcommonsubgraph(
        G, H; subgraphtype=:NodeInduced, algorithm=:Clique, kwargs...)
    if algorithm == :VF2
        # TODO: not implemented yet
        """
        if subgraphtype == :NodeInduced
            nodemcsvf2(G, H; kwargs...)
        elseif subgraphtype == :EdgeInduced
            edgemcsvf2(G, H; kwargs...)
        end
        """
    elseif algorithm == :Clique
        if subgraphtype == :NodeInduced
            nodemcsclique(G, H; kwargs...)
        elseif subgraphtype == :EdgeInduced
            edgemcsclique(G, H; kwargs...)
        end
    end
end


function delta_y_correction!(mapping, G, H)
    gsub = edgesubgraph(G, keys(mapping))
    hsub = edgesubgraph(H, values(mapping))
    revmap = Dict(v => k for (k, v) in mapping)
    for e in delta_y_edges(mapping, gsub, hsub)
        delete!(mapping, e)
    end
    for e in delta_y_edges(revmap, hsub, gsub)
        delete!(mapping, revmap[e])
    end
end

function delta_y_edges(mapping, gsub, hsub)
    dys = Int[]
    for gn in triangles(gsub)
        g_edges = edgekeys(nodesubgraph(gsub, gn))
        h_edges = [mapping[e] for e in g_edges]
        if nodecount(edgesubgraph(hsub, h_edges)) != 3
            push!(dys, pop!(g_edges))
        end
    end
    return dys
end

function delta_y_mismatch(G, H, mapping)
    # Delta-Y check for edge induced subgraph isomorphism
    gsub = edgesubgraph(G, keys(mapping))
    hsub = edgesubgraph(H, values(mapping))
    revmap = Dict(v => k for (k, v) in mapping)
    return delta_y(gsub, hsub, mapping) || delta_y(hsub, gsub, revmap)
end


function delta_y(gsub, hsub, mapping)
    for gn in triangles(gsub)
        g_edges = edgekeys(nodesubgraph(gsub, gn))
        h_edges = [mapping[e] for e in g_edges]
        if nodecount(edgesubgraph(hsub, h_edges)) != 3
            return true
        end
    end
    return false
end


function lgnodematcher(G::LineGraph, H::LineGraph,
                       nodematcher::Function, edgematcher::Function)
    return function (g, h)
        if !edgematcher(g, h)
            return false
        end
        ge = getnode(G, g)
        he = getnode(H, h)
        m1 = nodematcher(ge.n1, he.n1) && nodematcher(ge.n2, he.n2)
        m2 = nodematcher(ge.n1, he.n2) && nodematcher(ge.n2, he.n1)
        return m1 || m2
    end
end


function lgedgematcher(G::LineGraph, H::LineGraph, nodematcher::Function)
    return (g, h) -> nodematcher(getedge(G, g).node, getedge(H, h).node)
end
