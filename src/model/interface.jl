#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# Registered gprop types and converters for (de)serialization
const ELEMENT_TYPE_REGISTRY = Dict{String,Type}(
    "Int" => Int,
    "Any" => Any
)
const PROPERTY_TYPE_REGISTRY = Dict{String,Function}()


abstract type AbstractMolGraph{T} <: AbstractGraph{T} end
abstract type SimpleMolGraph{T,V,E} <: AbstractMolGraph{T} end  # mol graph that have SimpleGraph
abstract type AbstractReaction{T<:AbstractMolGraph} end


# Graph-level metadata properties (e.g. SDFile option fields)

struct Metadata <: AbstractDict{String,Any}
    mapping::OrderedDict{String,Any}
end

Metadata() = Metadata(OrderedDict{String,Any}())
Metadata(data::Vector) = Metadata(OrderedDict{String,Any}(String(row[1]) => row[2] for row in data))
Metadata(data::Dict) = Metadata([[String(k), v] for (k, v) in data])

Base.iterate(meta::Metadata, i...) = iterate(meta.mapping, i...)
Base.length(meta::Metadata) = length(meta.mapping)
Base.get(meta::Metadata, k, v) = get(meta.mapping, k, v)
Base.setindex!(meta::Metadata, v, k) = setindex!(meta.mapping, v, k)
Base.delete!(meta::Metadata, k) = delete!(meta.mapping, k)
to_dict(::Val{:default}, key::Symbol, meta::Metadata) = Dict{String,Any}(
    "key" => string(key),
    "type" => "Metadata",
    "data" => [[i, val] for (i, val) in meta]
)
PROPERTY_TYPE_REGISTRY["Metadata"] = (T, data) -> Metadata(data)


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


# SimpleGraph interface

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
u_edge(::Type{T}, src, dst) where T <: Integer = src < dst ? Edge{T}(src, dst) : Edge{T}(dst, src)
u_edge(g::SimpleGraph{T}, src, dst) where T = u_edge(T, src, dst)
u_edge(g, e) = u_edge(g, src(e), dst(e))
u_edge(e::Edge{T}) where T = u_edge(T, src(e), dst(e))

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

u_edge(mol::AbstractMolGraph{T}, src, dst) where T = u_edge(T, src, dst)
vproptype(::Type{<:SimpleMolGraph{T,V,E}}) where {T,V,E} = V
vproptype(mol::T) where T<:SimpleMolGraph = vproptype(T)
eproptype(::Type{<:SimpleMolGraph{T,V,E}}) where {T,V,E} = E
eproptype(mol::T) where T<:SimpleMolGraph = eproptype(T)
props(mol::SimpleMolGraph, v::Integer) = mol.vprops[v]
props(mol::SimpleMolGraph, e::Edge) = mol.eprops[e]
props(mol::SimpleMolGraph, u::Integer, v::Integer) = props(mol, u_edge(mol, u, v))
get_prop(mol::SimpleMolGraph, prop::Symbol) = mol.gprops[prop]
get_prop(mol::SimpleMolGraph, prop::String) = get_prop(mol, :metadata)[prop]
get_prop(mol::SimpleMolGraph, v::Integer, prop::Symbol) = props(mol, v)[prop]
get_prop(mol::SimpleMolGraph, e::Edge, prop::Symbol) = props(mol, e)[prop]
get_prop(mol::SimpleMolGraph, u::Integer, v::Integer, prop::Symbol) = props(mol, u, v)[prop]
has_prop(mol::SimpleMolGraph, prop::Symbol) = haskey(mol.gprops, prop)
has_prop(mol::SimpleMolGraph, prop::String) = haskey(get_prop(mol, :metadata), prop)
edge_rank(mol::SimpleMolGraph, e::Edge) = mol.edge_rank[e]
edge_rank(mol::SimpleMolGraph, u::Integer, v::Integer) = edge_rank(mol, u_edge(mol, u, v))
edge_neighbors(mol::AbstractMolGraph, u, v) = edge_neighbors(mol.graph, u, v)
edge_neighbors(mol::AbstractMolGraph, e) = edge_neighbors(mol.graph, e)

function set_prop!(mol::SimpleMolGraph{T,V,E}, v::T, value::V) where {T,V,E}
    mol.vprops[v] = value
    set_state!(mol, :has_updates, true)
end

function set_prop!(mol::SimpleMolGraph{T,V,E}, e::Edge{T}, value::E) where {T,V,E}
    mol.eprops[e] = value
    set_state!(mol, :has_updates, true)
end

function set_prop!(mol::SimpleMolGraph, prop::String, value)
    # Metadata update would not affect graph state
    mol.gprops[:metadata][prop] = value
end

function set_prop!(mol::SimpleMolGraph, prop::Symbol, value)
    # Note: this should not be called manually
    mol.gprops[prop] = value
    set_state!(mol, :has_updates, true)
end

Base.getindex(mol::SimpleMolGraph, k::String) = get_prop(mol, k)
Base.setindex!(mol::SimpleMolGraph, v, k::String) = set_prop!(mol, k, v)

function Base.show(io::IO, ::MIME"text/plain", g::SimpleMolGraph{T,V,E}) where {T,V,E}
    print(io, "{$(nv(g)), $(ne(g))} simple molecular graph $(typeof(g))")
end

Base.:(==)(g::SimpleMolGraph, h::SimpleMolGraph
    ) = g.graph == h.graph && g.vprops == h.vprops && g.eprops == h.eprops && g.gprops == h.gprops


function remap_gprops(mol::SimpleMolGraph, vmap)  # vmap[old] -> new
    newgp = Dict{Symbol,Any}()
    for (k, v) in mol.gprops
        newgp[k] = applicable(remap, v, vmap) ? remap(v, vmap) : v
    end
    return newgp
end

function remap_gprops!(mol::SimpleMolGraph, vmap)  # vmap[old] -> new
    for (k, v) in mol.gprops
        mol.gprops[k] = applicable(remap, v, vmap) ? remap(v, vmap) : v
    end
end


"""
    Reaction{T}

Reaction type.
"""
struct Reaction{T} <: AbstractReaction{T}
    reactants::Vector{T}
    products::Vector{T}
    rprops::Dict{Symbol,Any}
end

Reaction{T}() where T = Reaction{T}([], [], Dict())

Base.eltype(::Type{Reaction{T}}) where T = T
Base.eltype(rxn::Reaction{T}) where T = T
