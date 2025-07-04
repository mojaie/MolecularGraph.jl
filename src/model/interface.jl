#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    AbstractMolGraph{T<:Integer} <: Graphs.AbstractGraph{T}

The base class of molecular graphs
"""
abstract type AbstractMolGraph{T<:Integer} <: Graphs.AbstractGraph{T} end


# Graphs.jl common interface

Base.eltype(g::AbstractMolGraph) = eltype(g.graph)
# https://github.com/JuliaGraphs/Graphs.jl/pull/144
Base.eltype(::Type{<:AbstractMolGraph{T}}) where T = T
Base.zero(::Type{G}) where G <: AbstractMolGraph = G()
Graphs.edgetype(g::AbstractMolGraph) = edgetype(g.graph)
Graphs.edgetype(::Type{<:AbstractMolGraph{T}}) where T = Edge{T}
Graphs.nv(g::AbstractMolGraph) = nv(g.graph)
Graphs.ne(g::AbstractMolGraph) = ne(g.graph)
Graphs.vertices(g::AbstractMolGraph) = vertices(g.graph)
Graphs.edges(g::AbstractMolGraph) = edges(g.graph)
Graphs.is_directed(::Type{<:AbstractMolGraph}) = false
Graphs.has_vertex(g::AbstractMolGraph, x::Integer) = has_vertex(g.graph, x)
Graphs.has_edge(g::AbstractMolGraph, s::Integer, d::Integer) = has_edge(g.graph, s, d)
Graphs.inneighbors(g::AbstractMolGraph, v::Integer) = inneighbors(g.graph, v)
Graphs.outneighbors(g::AbstractMolGraph, v::Integer) = outneighbors(g.graph, v)


"""
    AbstractReaction

The base class of reactions
"""
abstract type AbstractReaction end


"""
    AbstractProperty

The base class of graph-level property of the molecular model.
"""
abstract type AbstractProperty end

function Base.:(==)(g::AbstractProperty, h::AbstractProperty)
    for sym in fieldnames(typeof(g))
        getproperty(g, sym) == getproperty(h, sym) || return false
    end
    return true
end


"""
    SimpleMolProperty{T<:Integer} <: AbstractProperty

The base class of graph-level property for SimpleMolGraph.
"""
abstract type SimpleMolProperty{T<:Integer} <: AbstractProperty end

Base.eltype(::Type{<:SimpleMolProperty{T}}) where T = T


"""
    AbstractState

The base class of molecular model states.
"""
abstract type AbstractState end


"""
    AbstractQueryNode

The base class of query node.
"""
abstract type AbstractQueryNode end

Base.getindex(elem::AbstractQueryNode, prop::Symbol) = getproperty(elem, prop)

Base.:(==)(x::T, y::T) where T <: AbstractQueryNode = all(
    [getfield(x, f) == getfield(y, f) for f in fieldnames(T)])

function Base.hash(elem::T, h::UInt) where T <: AbstractQueryNode
    for name in fieldnames(T)
        val = getfield(elem, name)
        h = hash(val, h)
    end
    return h
end


"""
    AbstractElement

The base class of vertex properties (atom) and edge properties (bond).
"""
abstract type AbstractElement end

Base.getindex(elem::AbstractElement, prop::Symbol) = getproperty(elem, prop)

Base.:(==)(x::T, y::T) where T <: AbstractElement = all(
    [getfield(x, f) == getfield(y, f) for f in fieldnames(T)])

function Base.hash(elem::T, h::UInt) where T <: AbstractElement
    for name in fieldnames(T)
        val = getfield(elem, name)
        h = hash(val, h)
    end
    return h
end


"""
    AbstractAtom <: AbstractElement

The base class of vertex properties (atom).
"""
abstract type AbstractAtom <: AbstractElement end


"""
    AbstractBond <: AbstractElement

The base class of edge properties (bond).
"""
abstract type AbstractBond <: AbstractElement end


"""
    QueryTree{T<:Integer,U<:AbstractQueryNode} <: AbstractElement

The base class of molecular query trees.
"""
abstract type QueryTree{T<:Integer,U<:AbstractQueryNode} <: AbstractElement end

Base.eltype(::Type{<:QueryTree{T,U}}) where {T,U} = T
Base.eltype(qtree::T) where T<:QueryTree = eltype(T)
vproptype(::Type{<:QueryTree{T,U}}) where {T,U} = U
vproptype(qtree::T) where T<:QueryTree = vproptype(T)


"""
    SimpleMolGraph{T<:Integer} <: AbstractMolGraph{T}

The base class of molecular graph models based on SimpleGraph
"""
abstract type SimpleMolGraph{T<:Integer} <: AbstractMolGraph{T} end


# SimpleMolGraph interface

Graphs.add_edge!(mol::SimpleMolGraph, u::Integer, v::Integer, prop
    ) = add_u_edge!(mol, u_edge(mol, u, v), prop)
Graphs.add_edge!(mol::SimpleMolGraph, e::Edge, prop
    ) = add_edge!(mol, src(e), dst(e), prop)
Graphs.rem_edge!(mol::SimpleMolGraph, u::Integer, v::Integer) = rem_u_edge!(mol, u_edge(mol, u, v))
Graphs.rem_edge!(mol::SimpleMolGraph, e::Edge) = rem_edge!(mol, src(e), dst(e))
Graphs.induced_subgraph(mol::T, vlist::AbstractVector{U}
    ) where {T<:SimpleMolGraph, U<:Integer} = _induced_subgraph(mol, vlist)
Graphs.induced_subgraph(mol::T, elist::AbstractVector{U}
    ) where {T<:SimpleMolGraph, U<:Edge} = _induced_subgraph(mol, elist)

function Base.show(io::IO, ::MIME"text/plain", g::SimpleMolGraph)
    print(io, "{$(nv(g)), $(ne(g))} simple molecular graph $(typeof(g))")
end


"""
    vproptype(::Type{SimpleMolGraph}) -> Type
    vproptype(mol::SimpleMolGraph) -> Type

Return the type of vertex properties
"""
vproptype(::Type{T}) where T<:SimpleMolGraph = vproptype(T)
vproptype(mol::T) where T<:SimpleMolGraph = vproptype(T)


"""
    eproptype(::Type{SimpleMolGraph}) -> Type
    eproptype(mol::SimpleMolGraph) -> Type

Return the type of edge properties
"""
eproptype(::Type{T}) where T<:SimpleMolGraph = eproptype(T)
eproptype(mol::T) where T<:SimpleMolGraph = eproptype(T)


"""
    u_edge(::Type{T}, src, dst) where T <: Integer -> Edge{T}
    u_edge(g::SimpleGraph{T}, src, dst) where T -> Edge{T}
    u_edge(mol::AbstractMolGraph{T}, src, dst) where T -> Edge{T}

A workaround for UndirectedEdge that are not yet implemented in SimpleGraph
"""
u_edge(::Type{T}, src::T, dst::T) where T<:Integer = src < dst ? Edge{T}(src, dst) : Edge{T}(dst, src)
u_edge(g::SimpleGraph{T}, src::T, dst::T) where T<:Integer = u_edge(T, src, dst)
u_edge(g::SimpleGraph, e::Edge) = u_edge(g, src(e), dst(e))
u_edge(e::Edge{T}) where T<:Integer = u_edge(T, src(e), dst(e))
u_edge(mol::AbstractMolGraph{T}, src::T, dst::T) where T<:Integer = u_edge(T, src, dst)
u_edge(mol::AbstractMolGraph, e::Edge) = u_edge(mol, src(e), dst(e))



"""
    edge_neighbors(g::SimpleGraph{T}, u::Integer, v::Integer) where T -> Tuple{Vector{T},Vector{T}}
    edge_neighbors(g::SimpleGraph{T}, e::Edge) where T -> Tuple{Vector{T},Vector{T}}

Return neighbors of the source and destination vertices of the edge, respectively.
"""
edge_neighbors(g::SimpleGraph, u::Integer, v::Integer) = (
    filter(n -> n != v, neighbors(g, u)),
    filter(n -> n != u, neighbors(g, v))
)
edge_neighbors(g::SimpleGraph, e::Edge) = edge_neighbors(g, src(e), dst(e))
edge_neighbors(mol::AbstractMolGraph, u, v) = edge_neighbors(mol.graph, u, v)
edge_neighbors(mol::AbstractMolGraph, e) = edge_neighbors(mol.graph, e)


# if specific index order is required.
# neighbors in SimpleGraph guarantees output index in lexicographic order now,
# but it is not formulated in API
ordered_neighbors = neighbors
ordered_edge_neighbors = edge_neighbors



# Reactions

@kwdef mutable struct ReactionProperty <: AbstractProperty
    # Reaction-level metadata properties (e.g. SDFile option fields)
    metadata::OrderedDict{String,String} = OrderedDict()
    # Parse errors
    logs::Dict{String,String} = Dict()
end


"""
    Reaction{T<:AbstractMolGraph}

Reaction type.
"""
struct Reaction{T<:AbstractMolGraph} <: AbstractReaction
    reactants::Vector{T}
    products::Vector{T}
    rprops::ReactionProperty
end

Reaction{T}() where T = Reaction{T}([], [], ReactionProperty())

Base.eltype(rxn::Reaction) = eltype(rxn.reactants)
Base.eltype(::Type{Reaction{T}}) where T = T