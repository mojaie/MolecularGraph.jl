#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function reconstruct!(::Val{T}, gprop, data) where T
    setproperty!(gprop, T, reconstruct(Val(T), gprop, data))
end

# vmap[old] -> new
function remap!(::Val{T}, gprop, vmap) where T
    setproperty!(gprop, T, remap(Val(T), gprop, vmap))
end


@kwdef mutable struct Descriptors{T}
    # cached relatively expensive descriptors
    sssr::Vector{Vector{T}} = []
    lone_pair::Vector{Int} = []
    apparent_valence::Vector{Int} = []
    valence::Vector{Int} = []
    is_ring_aromatic::Vector{Bool} = []
    # standardized atom charges and bond orders
    atom_charge::Vector{Int} = []
    bond_order::Vector{Int} = []
end

function Descriptors{T}(data::Dict{String,Any}) where T
    desc = Descriptors{T}()
    for (k, v) in data
        setproperty!(desc, k, v)
    end
    return desc
end



@kwdef mutable struct MolGraphProperty{T}
    stereocenter::Dict{T,Tuple{T,T,T,Bool}} = Dict()
    stereobond::Dict{Edge{T},Tuple{T,T,Bool}} = Dict()
    pyrrole_like::Vector{T} = []
    smarts_lexical_succ::Vector{Vector{T}} = []
    smarts_connectivity::Vector{Vector{T}} = []
    descriptors::Descriptors{T} = Descriptors{T}()
    coords2d::Vector{Vector{Point2d}} = []
    coords3d::Vector{Vector{Point3d}} = []
    # Graph-level metadata properties (e.g. SDFile option fields)
    metadata::OrderedDict{String,String} = OrderedDict()
    # Parse errors
    logs::Dict{String,String} = Dict()
end

function MolGraphProperty{T}(data::Dict{String,Any}) where T
    gprop = MolGraphProperty{T}()
    for sym in fieldnames(typeof(gprop))
        reconstruct!(Val(sym), gprop, data[string(sym)])
    end
    return gprop
end

function to_dict(::Val{:default}, gprop::MolGraphProperty)
    data = Dict{String,Any}()
    for sym in fieldnames(typeof(gprop))
        data[string(sym)] = to_dict(Val(sym), getproperty(gprop, sym))
    end
    """
    data["coords2d"] = [[collect(cd) for cd in cds] for cds in gprop.coords2d]
    data["coords3d"] = [[collect(cd) for cd in cds] for cds in gprop.coords3d]
    data["metadata"] = [[String(k), v] for (k, v) in gprop.metadata]
    data["error_sdfilereader"] = gprop.error_sdfilereader
    """
    return data
end


@kwdef mutable struct MolGraphState{T}
    initialized::Bool = false
    has_updates::Bool = true
    updater::Function = default_updater!
    on_init::Function = default_on_init!
    edge_rank::Dict{Edge{T},Int} = Dict{Edge{T},Int}()
end

function default_updater!(mol)
    update_edge_rank!(mol)
    reset_updates!(mol)
end

function default_on_init!(mol)
    # No initialization by default, just set flag
    set_state!(mol, :initialized, true)
end

function remap_gprops(mol::SimpleMolGraph{T,V,E}, vmap) where {T,V,E}
    gprop = MolGraphProperty{T}()
    for k in fieldnames(typeof(mol.gprops))
        setproperty!(gprop, k, remap(Val(k), mol.gprops, vmap))
    end
    return gprop
end

function remap_gprops!(mol::SimpleMolGraph, vmap)
    for k in fieldnames(typeof(mol.gprops))
        remap!(Val(k), mol.gprops, vmap)
    end
end


function to_dict(::Val{T}, ::Val{:descriptors}, gprop::MolGraphProperty) where T
    data = Dict{String,Any}()
    for sym in fieldnames(typeof(gprop.descriptors))
        data[string(sym)] = getproperty(gprop.descriptors, sym)
    end
    return data
end

function reconstruct(::Val{:descriptors}, gprop::MolGraphProperty{T}, data) where T
    desc = Descriptor{T}()
    for sym in fieldnames(typeof(gprop.descriptors))
        setproperty!(desc, sym, data[string(sym)])
    end
    return desc
end
# recalculate
remap(::Val{:descriptors}, gprop::MolGraphProperty, vmap) = gprop.descriptors


to_dict(
    ::Val{:default}, ::Val{:metadata}, gprop::MolGraphProperty) = gprop.metadata
reconstruct(::Val{:metadata}, gprop::MolGraphProperty, data) = data
remap(::Val{:metadata}, gprop::MolGraphProperty, vmap) = gprop.metadata

to_dict(
    ::Val{:default}, ::Val{:logs}, gprop::MolGraphProperty) = gprop.logs
reconstruct(::Val{:logs}, gprop::MolGraphProperty, data) = data
remap(::Val{:logs}, gprop::MolGraphProperty, vmap) = gprop.logs


"""
    MolGraph{T,V,E} <: SimpleMolGraph{T,V,E}

Basic molecular graph type.
"""
struct MolGraph{T,V,E} <: SimpleMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::MolGraphProperty{T}
    state::MolGraphState{T}
end

function MolGraph{T,V,E}(
        g::SimpleGraph, vprop::Dict, eprop::Dict;
        gprop=MolGraphProperty{T}(), config=MolGraphState{T}(), kwargs...
        ) where {T,V,E}
    (nv(g) > length(vprop)
        && error("Mismatch in the number of nodes and node properties"))
    (ne(g) == length(eprop)
        || error("Mismatch in the number of edges and edge properties"))
    # expand fadjlist for vprops of isolated nodes
    for _ in nv(g):(length(vprop) - 1)
        push!(g.fadjlist, T[])
    end
    mol = MolGraph{T,V,E}(g, vprop, eprop, gprop, config)
    mol.state.initialized || mol.state.on_init(mol)
    mol.state.has_updates && mol.state.updater(mol)
    return mol
end

MolGraph(
    g::SimpleGraph{T}, vprop::Dict{T,V}, eprop::Dict{Edge{T},E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(g, vprop, eprop; kwargs...)


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

dispatch!(mol::MolGraph, event::Symbol) = getproperty(mol.state, event)(mol)

function dispatch_update!(mol::MolGraph)
    mol.has_updates || return
    mol.state.updater(mol)
end

get_state(mol::MolGraph, sym::Symbol) = getproperty(mol.state, sym)
set_state!(mol::MolGraph, sym::Symbol, value) = setproperty!(mol.state, sym, value)

get_cache(mol::MolGraph, sym::Symbol) = getproperty(mol.state.caches, sym)
set_cache!(mol::MolGraph, sym::Symbol, value) = setproperty!(mol.state.caches, sym, value)

reset_updates!(mol) = setproperty!(mol.state, :has_updates, false)
has_updates(mol) = mol.state.has_updates
function notify_updates!(mol)
    # TODO: flag to inactivate update
    setproperty!(mol.state, :has_updates, true)
end

function update_edge_rank!(mol::MolGraph)
    empty!(mol.state.edge_rank)
    for (i, e) in enumerate(edges(mol.graph))
        mol.state.edge_rank[e] = i
    end
end


# Edit graph

function add_u_edge!(mol::MolGraph{T,V,E}, e::Edge, prop::E) where {T,V,E}
    # Can be directly called if src < dst is guaranteed.
    add_edge!(mol.graph, e) || return false
    mol.eprops[e] = prop
    notify_updates!(mol)
    return true
end


function Graphs.add_vertex!(mol::MolGraph{T,V,E}, prop::V) where {T,V,E}
    add_vertex!(mol.graph) || return false
    mol.vprops[nv(mol.graph)] = prop
    notify_updates!(mol)
    return true
end


function rem_u_edge!(mol::MolGraph, e::Edge)
    # Can be directly called if src < dst is guaranteed.
    rem_edge!(mol.graph, e) || return false
    delete!(mol.eprops, e)
    notify_updates!(mol)
    return true
end


function Graphs.rem_vertex!(mol::MolGraph, v::Integer)
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


function Graphs.rem_vertices!(mol::MolGraph{T,V,E}, vs::Vector{T}) where {T,V,E}
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


function _induced_subgraph(mol::T, vlist_or_elist) where {T<:MolGraph}
    # In many cases, induced_subgraph(mol.graph, vlist_or_elist) is sufficient
    subg, vmap = induced_subgraph(mol.graph, vlist_or_elist)
    vps = Dict(v => mol.vprops[vmap[v]] for v in vertices(subg))
    eps = Dict(e => mol.eprops[u_edge(mol, vmap[src(e)], vmap[dst(e)])]
        for e in edges(subg))
    newgp = remap_gprops(mol, Dict(v => i for (i, v) in enumerate(vmap)))
    notify_updates!(mol)
    return T(subg, vps, eps, gprop=newgp, state=mol.state), vmap
end
