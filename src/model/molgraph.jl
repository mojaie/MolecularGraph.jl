#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    GraphMol, QueryMol,
    graphmol, querymol,
    SDFile, SMILES, SMARTS


struct GraphMol{A<:Atom,B<:Bond} <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{A}
    edgeattrs::Vector{B}
    cache::Dict{Symbol,Any}
    attributes::Dict{Symbol,Any}
end

"""
    graphmol() -> GraphMol

Generate empty `GraphMol`.
"""
graphmol(::Type{A}, ::Type{B}
    ) where {A<:Atom,B<:Bond} = GraphMol{A,B}([], [], [], [], Dict(), Dict())

"""
    graphmol(atoms::Vector{Atom}, bonds::Vector{Bond}) -> GraphMol

Generate `GraphMol` that has the given atom objects and edge objects.
"""
function graphmol(
        edges, atoms::Vector{A}, bonds::Vector{B}) where {A<:Atom,B<:Bond}
    nbrmap = [Dict{Int,Int}() for i in 1:length(atoms)]
    edges = collect(edges)
    for (i, (u, v)) in enumerate(edges)
        nbrmap[u][i] = v
        nbrmap[v][i] = u
    end
    return GraphMol(
        nbrmap, edges, atoms, bonds, Dict{Symbol,Any}(), Dict{Symbol,Any}())
end

"""
    graphmol(mol::GraphMol) -> GraphMol

Copy `GraphMol`.
"""
function graphmol(graph::GraphMol{A,B}) where {A<:Atom,B<:Bond}
    newg = graphmol(A, B)
    for nbr in graph.neighbormap
        push!(newg.neighbormap, copy(nbr))
    end
    append!(newg.edges, graph.edges)
    append!(newg.nodeattrs, graph.nodeattrs)
    append!(newg.edgeattrs, graph.edgeattrs)
    merge!(newg.cache, graph.cache)
    merge!(newg.attributes, graph.attributes)
    return newg
end

"""
    graphmol(mol::SubgraphView{GraphMol}) -> GraphMol

Generate a new `GraphMol` from a substructure view.

Graph property caches and attributes are not inherited.
"""
function graphmol(view::SubgraphView)
    newg = graphmol(nodeattrtype(view), edgeattrtype(view))
    nkeys = sort(collect(nodeset(view)))
    ekeys = sort(collect(edgeset(view)))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        push!(newg.nodeattrs, nodeattr(view, n))
        push!(newg.neighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        (oldu, oldv) = getedge(view, e)
        u = nmap[oldu]
        v = nmap[oldv]
        push!(newg.edges, (u, v))
        push!(newg.edgeattrs, edgeattr(view, e))
        newg.neighbormap[u][i] = v
        newg.neighbormap[v][i] = u
    end
    return newg
end



struct QueryMol{A<:QueryAtom,B<:QueryBond} <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{A}
    edgeattrs::Vector{B}
    cache::Dict{Symbol,Any}
    attributes::Dict{Symbol,Any}
    connectivity::Vector{Vector{Int}}
end

"""
    querymol() -> QueryMol

Generate empty `QueryMol`.
"""
querymol(
    ::Type{A}, ::Type{B}
) where {A<:QueryAtom,B<:QueryBond} = QueryMol{A,B}(
    [], [], [], [], Dict(), Dict(), [])

"""
    querymol(atoms::Vector{Atom}, bonds::Vector{Bond}) -> GraphMol

Generate `QueryMol` that has the given atom objects and edge objects.
"""
function querymol(edges, atoms::Vector{A}, bonds::Vector{B},
        connectivity::Vector{Vector{Int}}) where {A<:QueryAtom,B<:QueryBond}
    nbrmap = [Dict{Int,Int}() for i in 1:length(atoms)]
    edges = collect(edges)
    for (i, (u, v)) in enumerate(edges)
        nbrmap[u][i] = v
        nbrmap[v][i] = u
    end
    return QueryMol(
        nbrmap, edges, atoms, bonds,
        Dict{Symbol,Any}(), Dict{Symbol,Any}(), connectivity)
end


# Aliases

SDFile = GraphMol{SDFileAtom,SDFileBond}
SMILES = GraphMol{SmilesAtom,SmilesBond}
SMARTS = QueryMol{SmartsAtom,SmartsBond}
