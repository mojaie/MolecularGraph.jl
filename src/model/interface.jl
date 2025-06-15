#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    AbstractMolGraph{T} <: Graphs.AbstractGraph{T}

The base class of molecular graphs
"""
abstract type AbstractMolGraph{T<:Integer} <: Graphs.AbstractGraph{T} end

"""
    SimpleMolGraph{T,V,E} <: AbstractMolGraph{T}

The base class of molecular graph models based on SimpleGraph
"""
abstract type SimpleMolGraph{T<:Integer,V,E} <: AbstractMolGraph{T} end

"""
    AbstractReaction{T<:AbstractMolGraph}

The base class of reactions
"""
abstract type AbstractReaction{T<:AbstractMolGraph} end


"""
    QueryTree{T,U}

The base class of reactions
"""
abstract type QueryTree{T<:Integer,U} end

const QueryMolGraph = SimpleMolGraph{<:Integer,<:QueryTree,<:QueryTree}


# Graphs.jl common interface

Graphs.edgetype(g::AbstractMolGraph) = edgetype(g.graph)
Graphs.edgetype(::Type{<:AbstractMolGraph{T}}) where T = Edge{T}
Base.eltype(g::AbstractMolGraph) = eltype(g.graph)
# https://github.com/JuliaGraphs/Graphs.jl/pull/144
Base.eltype(::Type{<:AbstractMolGraph{T}}) where T = T
Graphs.nv(g::AbstractMolGraph) = nv(g.graph)
Graphs.ne(g::AbstractMolGraph) = ne(g.graph)
Graphs.vertices(g::AbstractMolGraph) = vertices(g.graph)
Graphs.edges(g::AbstractMolGraph) = edges(g.graph)
Graphs.is_directed(::Type{<:AbstractMolGraph{T}}) where T = false
Graphs.has_vertex(g::AbstractMolGraph, x::Integer) = has_vertex(g.graph, x)
Graphs.has_edge(g::AbstractMolGraph, s::Integer, d::Integer) = has_edge(g.graph, s, d)
Graphs.inneighbors(g::AbstractMolGraph, v::Integer) = inneighbors(g.graph, v)
Graphs.outneighbors(g::AbstractMolGraph, v::Integer) = outneighbors(g.graph, v)
Base.zero(::Type{G}) where G <: AbstractMolGraph = G()


# SimpleMolGraph interface

Base.copy(mol::SimpleMolGraph) = deepcopy(mol)
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


"""
    edge_neighbors(g::SimpleGraph{T}, u::Integer, v::Integer) where T -> Tuple{Vector{T},Vector{T}}
    edge_neighbors(g::SimpleGraph{T}, e::Edge) where T -> Tuple{Vector{T},Vector{T}}

Return neighbors of the source and destination vertices of the edge, respectively
"""
edge_neighbors(g::SimpleGraph, u::Integer, v::Integer) = (
    filter(n -> n != v, neighbors(g, u)),
    filter(n -> n != u, neighbors(g, v))
)
edge_neighbors(g::SimpleGraph, e::Edge) = edge_neighbors(g, src(e), dst(e))

# if specific index order is required.
# neighbors in SimpleGraph guarantees output index in lexicographic order now,
# but it is not formulated in API
ordered_neighbors = neighbors
ordered_edge_neighbors = edge_neighbors


# SimpleMolGraph interface (node/edge primary attributes)

u_edge(mol::AbstractMolGraph{T}, src::T, dst::T) where T<:Integer = u_edge(T, src, dst)
u_edge(mol::AbstractMolGraph, e::Edge) = u_edge(mol, src(e), dst(e))
get_prop(mol::SimpleMolGraph, prop::Symbol) = getproperty(mol.gprops, prop)
get_prop(mol::SimpleMolGraph, v::Integer, prop::Symbol) = props(mol, v)[prop]
get_prop(mol::SimpleMolGraph, e::Edge, prop::Symbol) = props(mol, e)[prop]
get_prop(mol::SimpleMolGraph, u::Integer, v::Integer, prop::Symbol) = props(mol, u, v)[prop]
edge_rank(mol::SimpleMolGraph, e::Edge) = mol.state.edge_rank[e]
edge_rank(mol::SimpleMolGraph, u::Integer, v::Integer) = edge_rank(mol, u_edge(mol, u, v))
edge_neighbors(mol::AbstractMolGraph, u, v) = edge_neighbors(mol.graph, u, v)
edge_neighbors(mol::AbstractMolGraph, e) = edge_neighbors(mol.graph, e)

function set_prop!(mol::SimpleMolGraph{T,V,E}, v::T, value::V) where {T,V,E}
    mol.vprops[v] = value
    notify_updates!(mol)
end

function set_prop!(mol::SimpleMolGraph{T,V,E}, e::Edge{T}, value::E) where {T,V,E}
    mol.eprops[e] = value
    notify_updates!(mol)
end


function Base.show(io::IO, ::MIME"text/plain", g::SimpleMolGraph)
    print(io, "{$(nv(g)), $(ne(g))} simple molecular graph $(typeof(g))")
end

Base.:(==)(g::SimpleMolGraph, h::SimpleMolGraph
    ) = g.graph == h.graph && g.vprops == h.vprops && g.eprops == h.eprops && g.gprops == h.gprops


@kwdef mutable struct ReactionProperty
    # Reaction-level metadata properties (e.g. SDFile option fields)
    metadata::OrderedDict{String,String} = OrderedDict()
    # Parse errors
    logs::Dict{String,String} = Dict()
end

function Base.:(==)(g::ReactionProperty, h::ReactionProperty)
    for sym in fieldnames(typeof(g))
        getproperty(g, sym) == getproperty(h, sym) || return false
    end
    return true
end

"""
    Reaction{T}

Reaction type.
"""
struct Reaction{T} <: AbstractReaction{T}
    reactants::Vector{T}
    products::Vector{T}
    rprops::ReactionProperty
end

Reaction{T}() where T = Reaction{T}([], [], ReactionProperty())

Base.eltype(::Type{Reaction{T}}) where T = T
Base.eltype(rxn::Reaction{T}) where T = T
