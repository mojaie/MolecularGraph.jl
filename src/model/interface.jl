#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

import Graphs:
    AbstractGraph, edgetype, nv, ne, vertices, edges, is_directed,
    has_vertex, has_edge, inneighbors, outneighbors

export
    AbstractMolGraph, SimpleMolGraph, AbstractReaction,
    ordered_neighbors, undirectededge,
    vproptype, eproptype, props, vprops, eprops, get_prop, has_prop, edge_rank,
    edge_neighbors, ordered_edge_neighbors


abstract type AbstractMolGraph{T} <: AbstractGraph{T} end
abstract type SimpleMolGraph{T,V,E} <: AbstractMolGraph{T} end  # mol graph that have SimpleGraph
abstract type AbstractReaction{T<:AbstractMolGraph} end


# Graphs.jl common interface

edgetype(g::AbstractMolGraph) = edgetype(g.graph)
edgetype(::Type{<:AbstractMolGraph{T}}) where T = Edge{T}
Base.eltype(g::AbstractMolGraph) = eltype(g.graph)
# https://github.com/JuliaGraphs/Graphs.jl/pull/144
Base.eltype(::Type{<:AbstractMolGraph{T}}) where T = T
nv(g::AbstractMolGraph) = nv(g.graph)
ne(g::AbstractMolGraph) = ne(g.graph)
vertices(g::AbstractMolGraph) = vertices(g.graph)
edges(g::AbstractMolGraph) = edges(g.graph)
is_directed(::Type{<:AbstractMolGraph{T}}) where T = false
has_vertex(g::AbstractMolGraph, x::Integer) = has_vertex(g.graph, x)
has_edge(g::AbstractMolGraph, s::Integer, d::Integer) = has_edge(g.graph, s, d)
inneighbors(g::AbstractMolGraph, v::Integer) = inneighbors(g.graph, v)
outneighbors(g::AbstractMolGraph, v::Integer) = outneighbors(g.graph, v)
Base.zero(::Type{G}) where G <: AbstractMolGraph = G()


# if specific index order is required.
# neighbors in SimpleGraph guarantees output index in rexicographic order now,
# but it is not formulated in API
ordered_neighbors = neighbors

"""
    undirectededge(::Type{T}, src, dst) where T <: Integer -> Edge{T}
    undirectededge(g::SimpleGraph{T}, src, dst) where T -> Edge{T}
    undirectededge(mol::AbstractMolGraph{T}, src, dst) where T -> Edge{T}

A workaround for UndirectedEdge that are not yet implemented in SimpleGraph
"""
undirectededge(::Type{T}, src, dst) where T <: Integer = src < dst ? Edge{T}(src, dst) : Edge{T}(dst, src)
undirectededge(g::SimpleGraph{T}, src, dst) where T = undirectededge(T, src, dst)
undirectededge(mol::AbstractMolGraph{T}, src, dst) where T = undirectededge(T, src, dst)


# SimpleMolGraph interface (node/edge primary attributes)
Base.copy(mol::SimpleMolGraph) = deepcopy(mol)
vproptype(::Type{<:SimpleMolGraph{T,V,E}}) where {T,V,E} = V
vproptype(mol::T) where T<:SimpleMolGraph = vproptype(T)
eproptype(::Type{<:SimpleMolGraph{T,V,E}}) where {T,V,E} = E
eproptype(mol::T) where T<:SimpleMolGraph = eproptype(T)
vprops(mol::SimpleMolGraph) = mol.vprops
eprops(mol::SimpleMolGraph) = mol.eprops
props(mol::SimpleMolGraph) = mol.gprops
props(mol::SimpleMolGraph, v::Integer) = vprops(mol)[v]
props(mol::SimpleMolGraph, u::Integer, v::Integer) = props(mol, undirectededge(mol, u, v))
get_prop(mol::SimpleMolGraph, prop::Symbol) = props(mol)[prop]
get_prop(mol::SimpleMolGraph, v::Integer, prop::Symbol) = props(mol, v)[prop]
get_prop(mol::SimpleMolGraph, e::Edge, prop::Symbol) = props(mol, e)[prop]
get_prop(mol::SimpleMolGraph, u::Integer, v::Integer, prop::Symbol) = props(mol, u, v)[prop]
has_prop(mol::SimpleMolGraph, prop::Symbol) = haskey(props(mol), prop)
edge_rank(mol::SimpleMolGraph, u::Integer, v::Integer) = edge_rank(mol, undirectededge(mol, u, v))


# convenient functions

edge_neighbors(mol::AbstractMolGraph, u::Integer, v::Integer) = (
    filter(n -> n != v, neighbors(mol, u)),
    filter(n -> n != u, neighbors(mol, v))
)
edge_neighbors(mol::AbstractMolGraph, e::Edge) = edge_neighbors(mol, src(e), dst(e))

ordered_edge_neighbors = edge_neighbors