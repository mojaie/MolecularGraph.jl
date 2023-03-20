#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

import Graphs:
    AbstractGraph, edgetype, nv, ne, vertices, edges, is_directed,
    has_vertex, has_edge, inneighbors, outneighbors,
    add_vertex!, add_edge!, rem_vertex!, rem_vertices!, rem_edge!,
    induced_subgraph

export
    AbstractMolGraph, SimpleMolGraph,
    MolGraph, SDFMolGraph, SMILESMolGraph,
    MolGraphGen,
    AbstractReaction, Reaction,
    to_dict, to_json,
    ordered_neighbors, undirectededge, edge_rank,
    vproptype, eproptype,
    props, vprops, eprops,
    get_prop, has_prop,
    add_edges!, rem_edges!,
    descriptors, set_descriptor!, has_descriptor, get_descriptor,
    init_node_descriptor, init_edge_descriptor,
    set_prop!,
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


# SimpleGraph interface
Base.copy(mol::SimpleMolGraph) = deepcopy(mol)
add_edge!(mol::SimpleMolGraph, u::Integer, v::Integer, prop
    ) = add_u_edge!(mol, undirectededge(mol, u, v), prop)
add_edge!(mol::SimpleMolGraph, e::Edge, prop
    ) = add_edge!(mol, src(e), dst(e), prop)
rem_edge!(mol::SimpleMolGraph, u::Integer, v::Integer) = rem_u_edge!(mol, undirectededge(mol, u, v))
rem_edge!(mol::SimpleMolGraph, e::Edge) = rem_edge!(mol, src(e), dst(e))
induced_subgraph(mol::T, vlist::AbstractVector{U}
    ) where {T<:SimpleMolGraph, U<:Integer} = _induced_subgraph(mol, vlist)
induced_subgraph(mol::T, elist::AbstractVector{U}
    ) where {T<:SimpleMolGraph, U<:Edge} = _induced_subgraph(mol, elist)



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
undirectededge(g, e) = undirectededge(g, src(e), dst(e))
undirectededge(e::Edge{T}) where T = undirectededge(T, src(e), dst(e))

# SimpleMolGraph interface (node/edge primary attributes)
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
# get_prop(mol::SimpleMolGraph, prop::Symbol, default) = get(props(mol), prop, default)
get_prop(mol::SimpleMolGraph, v::Integer, prop::Symbol) = props(mol, v)[prop]
get_prop(mol::SimpleMolGraph, e::Edge, prop::Symbol) = props(mol, e)[prop]
get_prop(mol::SimpleMolGraph, u::Integer, v::Integer, prop::Symbol) = props(mol, u, v)[prop]
has_prop(mol::SimpleMolGraph, prop::Symbol) = haskey(props(mol), prop)
edge_rank(mol::SimpleMolGraph, u::Integer, v::Integer) = edge_rank(mol, undirectededge(mol, u, v))


add_edges!(mol::SimpleMolGraph, elist, plist) = add_u_edges!(mol,
    [undirectededge(mol, src(e), dst(e)) for e in elist], plist)

function set_prop!(mol::SimpleMolGraph{T,V,E}, v::T, value::V) where {T,V,E}
    mol.vprops[v] = value
    return value
end

function set_prop!(mol::SimpleMolGraph{T,V,E}, e::Edge{T}, value::E) where {T,V,E}
    mol.eprops[e] = value
    return value
end

function set_prop!(mol::SimpleMolGraph, prop::Symbol, value)
    mol.gprops[prop] = value
    return value
end

function Base.show(io::IO, ::MIME"text/plain", g::SimpleMolGraph{T,V,E}) where {T,V,E}
    print(io, "{$(nv(g)), $(ne(g))} simple molecular graph $(typeof(g))")
end


Base.:(==)(g::AbstractMolGraph, h::AbstractMolGraph
    ) = g.graph == h.graph && g.vprops == h.vprops && g.eprops == h.eprops && g.gprops == h.gprops


# convenient functions

edge_neighbors(mol::AbstractMolGraph, u::Integer, v::Integer) = (
    filter(n -> n != v, neighbors(mol, u)),
    filter(n -> n != u, neighbors(mol, v))
)
edge_neighbors(mol::AbstractMolGraph, e::Edge) = edge_neighbors(mol, src(e), dst(e))

ordered_edge_neighbors = edge_neighbors


"""
    to_dict(mol::MolGraph) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
to_dict

"""
    to_json(mol::MolGraph) -> String

Convert molecule object into JSON compatible dictionary.
"""
to_json


# Reaction

struct Reaction{T} <: AbstractReaction{T}
    reactants::Vector{T}
    products::Vector{T}
    rprops::Dict{Symbol,Any}
end

Reaction{T}() where T = Reaction{T}([], [], Dict())

Base.eltype(::Type{Reaction{T}}) where T = T
Base.eltype(rxn::Reaction{T}) where T = T
