#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# ReactiveMolGraph interfaces


# Update mechanisms


function update_edge_rank!(mol::ReactiveMolGraph)
    empty!(mol.state.edge_rank)
    for (i, e) in enumerate(edges(mol.graph))
        mol.state.edge_rank[e] = i
    end
    return
end


function dispatch_update!(mol::ReactiveMolGraph)
    mol.state.has_updates || return
    update_edge_rank!(mol)  # necessary
    mol.state.has_updates = false  # call before update to avoid infinite roop
    mol.state.on_update(mol)
    mol.state.has_new_edges = false
    return
end


function notify_updates!(mol::ReactiveMolGraph)
    # TODO: flag to inactivate auto-update
    mol.state.has_updates = true
    return
end

function notify_new_edges!(mol::ReactiveMolGraph)
    # TODO: flag to inactivate auto-update
    mol.state.has_new_edges = true
    return
end


function initialize!(mol::ReactiveMolGraph)
    mol.state.initialized || mol.state.on_init(mol)
    mol.state.initialized = true
    dispatch_update!(mol)  # checking on-update callback errors
    return
end


function default_on_init!(mol::ReactiveMolGraph)
    return  # No initialization by default
end


function default_on_update!(mol::ReactiveMolGraph)
    return  # No auto-update callbacks by default
end


Base.:(==)(g::ReactiveMolGraph, h::ReactiveMolGraph
    ) = g.graph == h.graph && g.vprops == h.vprops && g.eprops == h.eprops && g.gprops == h.gprops

vproptype(::Type{<:ReactiveMolGraph{T,V,E}}) where {T,V,E} = V
vproptype(::Type{T}) where T<:SimpleMolGraph = vproptype(T)
vproptype(mol::T) where T<:SimpleMolGraph = vproptype(T)
eproptype(::Type{<:ReactiveMolGraph{T,V,E}}) where {T,V,E} = E
eproptype(::Type{T}) where T<:SimpleMolGraph = eproptype(T)
eproptype(mol::T) where T<:SimpleMolGraph = eproptype(T)

props(mol::ReactiveMolGraph, v::Integer) = mol.vprops[v]
props(mol::ReactiveMolGraph, e::Edge) = mol.eprops[e]
props(mol::ReactiveMolGraph, u::Integer, v::Integer) = props(mol, u_edge(mol, u, v))

get_prop(mol::ReactiveMolGraph, v::Integer, prop::Symbol) = props(mol, v)[prop]
get_prop(mol::ReactiveMolGraph, e::Edge, prop::Symbol) = props(mol, e)[prop]
get_prop(mol::ReactiveMolGraph, u::Integer, v::Integer, prop::Symbol) = props(mol, u, v)[prop]



"""
    edge_rank(mol::SimpleMolGraph, e::Edge) -> Integer
    edge_rank(mol::SimpleMolGraph, u::Integer, v::Integer) -> Integer

A workaround for Edge property
"""
function edge_rank(mol::ReactiveMolGraph, e::Edge)
    dispatch_update!(mol)
    return mol.state.edge_rank[e]
end

edge_rank(mol::ReactiveMolGraph, u::Integer, v::Integer) = edge_rank(mol, u_edge(mol, u, v))


# Edit graph topology and properties

function set_prop!(mol::ReactiveMolGraph{T,V,E}, v::T, value::V) where {T,V,E}
    mol.vprops[v] = value
    notify_updates!(mol)
end

function set_prop!(mol::ReactiveMolGraph{T,V,E}, e::Edge{T}, value::E) where {T,V,E}
    mol.eprops[e] = value
    notify_updates!(mol)
end


function add_u_edge!(mol::ReactiveMolGraph{T,V,E}, e::Edge{T}, prop::E) where {T,V,E}
    # Can be directly called if src < dst is guaranteed.
    add_edge!(mol.graph, e) || return false
    mol.eprops[e] = prop
    notify_updates!(mol)
    notify_new_edges!(mol)
    return true
end


function Graphs.add_vertex!(mol::ReactiveMolGraph{T,V,E}, prop::V) where {T,V,E}
    add_vertex!(mol.graph) || return false
    mol.vprops[nv(mol.graph)] = prop
    notify_updates!(mol)
    return true
end


function rem_u_edge!(mol::ReactiveMolGraph, e::Edge)
    # Can be directly called if src < dst is guaranteed.
    rem_edge!(mol.graph, e) || return false
    delete!(mol.eprops, e)
    notify_updates!(mol)
    return true
end


function Graphs.rem_vertex!(mol::ReactiveMolGraph, v::Integer)
    nv_ = nv(mol)
    edges_ = collect(edges(mol))
    rem_vertex!(mol.graph, v) || return false
    # last index node is re-indexed to the removed node
    nv_ == v || begin mol.vprops[v] = mol.vprops[nv_] end
    delete!(mol.vprops, nv_)
    for e in edges_
        (src(e) == v || dst(e) == v) && delete!(mol.eprops, e)
	    ((src(e) == v && dst(e) == nv_) || (src(e) == nv_ && dst(e) == v)) && continue
        nv_ == v && continue
        src(e) == nv_ && begin mol.eprops[u_edge(mol, v, dst(e))] = mol.eprops[e]; delete!(mol.eprops, e) end
        dst(e) == nv_ && begin mol.eprops[u_edge(mol, src(e), v)] = mol.eprops[e]; delete!(mol.eprops, e) end
    end
    # remap TODO: refactoring
    if nv_ != v
        mapper = Dict(i => i for i in vertices(mol))
        delete!(mapper, v)
        mapper[nv_ ] = v
        remap_gprops!(mol, mapper)
    end
    notify_updates!(mol)
    return true
end


function Graphs.rem_vertices!(mol::ReactiveMolGraph{T,V,E}, vs::Vector{T}) where {T,V,E}
    # TODO: if many vertices should be removed, induced_subgraph may be more efficient.
    vmap = rem_vertices!(mol.graph, vs)
    for (i, v) in enumerate(vmap)
        i == v && continue
        mol.vprops[i] = mol.vprops[v]
    end
    for v in (nv(mol)+1):(nv(mol)+length(vs))
        delete!(mol.vprops, v)
    end
    for e in edges(mol)
        _e = u_edge(T, vmap[src(e)], vmap[dst(e)])
        e == _e && continue
        mol.eprops[e] = mol.eprops[_e]
    end
    new_edges = collect(edges(mol))
    for e in keys(mol.eprops)
        e in new_edges || delete!(mol.eprops, e)
    end
    remap_gprops!(mol, Dict(v => i for (i, v) in enumerate(vmap)))
    notify_updates!(mol)
    return vmap
end


function _induced_subgraph(mol::T, vlist_or_elist) where {T<:ReactiveMolGraph}
    # In many cases, induced_subgraph(mol.graph, vlist_or_elist) is sufficient
    subg, vmap = induced_subgraph(mol.graph, vlist_or_elist)
    vps = Dict(v => mol.vprops[vmap[v]] for v in vertices(subg))
    eps = Dict(e => mol.eprops[u_edge(mol, vmap[src(e)], vmap[dst(e)])]
        for e in edges(subg))
    newgp = remap_gprops(mol, Dict(v => i for (i, v) in enumerate(vmap)))
    notify_updates!(mol)
    return T(subg, vps, eps, newgp, mol.state), vmap
end



mutable struct MolState{T,F1,F2} <: AbstractState
    initialized::Bool
    has_updates::Bool
    has_new_edges::Bool  # as a SSSR recalculation flag
    on_init::F1
    on_update::F2
    edge_rank::Dict{Edge{T},Int}
end


function MolState{T}(
        ; initialized=false, has_updates=true, has_new_edges=true,
        on_init=default_on_init!, on_update=default_on_update!,
        edge_rank=Dict{Edge{T},Int}()) where T
    return MolState{T,typeof(on_init),typeof(on_update)}(
        initialized, has_updates, has_new_edges, on_init, on_update, edge_rank)
end


function reactive_molgraph(
        g::SimpleGraph{T}, vprops::Dict{T,V}, eprops::Dict{Edge{T},E};
        gprops=MolProperty{T}(), kwargs...) where {T,V,E}
    if nv(g) > length(vprops)
        error("Mismatch in the number of nodes and node properties")
    elseif ne(g) != length(eprops)
        error("Mismatch in the number of edges and edge properties")
    end
    # expand fadjlist for vprops of isolated nodes
    for _ in nv(g):(length(vprops) - 1)
        push!(g.fadjlist, T[])
    end
    config = MolState{T}(;kwargs...)
    return (g, vprops, eprops, gprops, config)
end

# from edge and property list (sdftomol, smilestomol interface)

function reactive_molgraph(
        edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}
        ; kwargs...) where {T,V,E}
    g = SimpleGraph(edge_list)
    vps = Dict{T,V}(i => v for (i, v) in enumerate(vprop_list))
    # eprop_list in edge_list order
    eps = Dict{Edge{T},E}(e => eprop_list[i] for (i, e) in enumerate(edge_list))
    return reactive_molgraph(g, vps, eps; kwargs...)
end


"""
    MolGraph{T,V,E} <: ReactiveMolGraph{T,V,E}

Basic molecular graph type.
"""
struct MolGraph{T<:Integer,V<:AbstractAtom,E<:AbstractBond} <: ReactiveMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::MolProperty{T}
    state::MolState{T}
end

function MolGraph{T,V,E}(args...; kwargs...) where {T,V,E}
    mol = MolGraph{T,V,E}(reactive_molgraph(args...; kwargs...)...)
    mol.state.initialized || mol.state.on_init(mol)
    mol.state.initialized = true
    initialize!(mol)
    return mol
end

MolGraph(
    g::SimpleGraph{T}, vprops::Dict{T,V}, eprops::Dict{Edge{T},E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(g, vprops, eprops; kwargs...)
MolGraph(
    edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(edge_list, vprop_list, eprop_list; kwargs...)

MolGraph{T,V,E}() where {T,V,E} = MolGraph(SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}())
MolGraph() = MolGraph{Int,SDFAtom,SDFBond}()


# MolGraph type aliases

const SDFMolGraph = MolGraph{Int,SDFAtom,SDFBond}
const SMILESMolGraph = MolGraph{Int,SMILESAtom,SMILESBond}
const CommonChemMolGraph = MolGraph{Int,CommonChemAtom,CommonChemBond}
