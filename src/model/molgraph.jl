#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@kwdef mutable struct MolGraphState{T}
    initialized::Bool = false
    has_updates::Bool = true
    on_init::Function = default_on_init!
    on_update::Function = default_on_update!
    edge_rank::Dict{Edge{T},Int} = Dict{Edge{T},Int}()
end

function MolGraphState{T}(on_init, on_update) where T
    state = MolGraphState{T}()
    state.on_init = on_init
    state.on_update = on_update
    return state
end

function default_on_init!(mol)
    # No initialization by default
end

function default_on_update!(mol)
    # These two emthods are necessary
end


"""
    MolGraph{T,V,E} <: ReactiveMolGraph{T,V,E}

Basic molecular graph type.
"""
struct MolGraph{T,V,E} <: ReactiveMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::MolGraphProperty{T}
    state::MolGraphState{T}
end

function MolGraph{T,V,E}(
        g::SimpleGraph, vprops::Dict, eprops::Dict;
        gprops=MolGraphProperty{T}(),
        on_init=default_on_init!,
        on_update=default_on_update!, kwargs...) where {T,V,E}
    (nv(g) > length(vprops)
        && error("Mismatch in the number of nodes and node properties"))
    (ne(g) == length(eprops)
        || error("Mismatch in the number of edges and edge properties"))
    # expand fadjlist for vprops of isolated nodes
    for _ in nv(g):(length(vprops) - 1)
        push!(g.fadjlist, T[])
    end
    config = MolGraphState{T}(on_init, on_update)
    mol = MolGraph{T,V,E}(g, vprops, eprops, gprops, config)
    mol.state.initialized || mol.state.on_init(mol)
    mol.state.initialized = true
    dispatch_update!(mol)
    return mol
end

MolGraph(
    g::SimpleGraph{T}, vprops::Dict{T,V}, eprops::Dict{Edge{T},E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(g, vprops, eprops; kwargs...)


# from empty set

MolGraph{T,V,E}() where {T,V,E} = MolGraph(SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}())
MolGraph() = MolGraph{Int,Any,Any}()


# from edge and property list (sdftomol, smilestomol interface)

function MolGraph{T,V,E}(
        edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}; kwargs...
        ) where {T,V,E}
    g = SimpleGraph(edge_list)
    vps = Dict(i => v for (i, v) in enumerate(vprop_list))
    eps = Dict(e => eprop_list[i] for (i, e) in enumerate(edge_list))  # eprop_list in edge_list order
    return MolGraph{T,V,E}(g, vps, eps; kwargs...)
end

MolGraph(
    edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(edge_list, vprop_list, eprop_list; kwargs...)


# MolGraph type aliases

const SDFMolGraph = MolGraph{Int,SDFAtom,SDFBond}
const SMILESMolGraph = MolGraph{Int,SMILESAtom,SMILESBond}
const CommonChemMolGraph = MolGraph{Int,CommonChemAtom,CommonChemBond}


# Update mechanisms

function dispatch_update!(mol::ReactiveMolGraph)
    mol.state.has_updates || return
    update_edge_rank!(mol)  # necessary
    mol.state.has_updates = false  # call before update to avoid infinite roop
    mol.state.on_update(mol)
end

function notify_updates!(mol::ReactiveMolGraph)
    # TODO: flag to inactivate auto-update
    mol.state.has_updates = true
end

function update_edge_rank!(mol::ReactiveMolGraph)
    empty!(mol.state.edge_rank)
    for (i, e) in enumerate(edges(mol.graph))
        mol.state.edge_rank[e] = i
    end
end


# Edit graph

function add_u_edge!(mol::ReactiveMolGraph{T,V,E}, e::Edge{T}, prop::E) where {T,V,E}
    # Can be directly called if src < dst is guaranteed.
    add_edge!(mol.graph, e) || return false
    mol.eprops[e] = prop
    notify_updates!(mol)
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
