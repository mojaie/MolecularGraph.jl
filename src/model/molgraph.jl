#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    MolGraph{T,V,E} <: SimpleMolGraph{T,V,E}

Basic molecular graph type.
"""
struct MolGraph{T,V,E} <: SimpleMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::Dict{Symbol,Any}
    state::Dict{Symbol,Any}
    edge_rank::Dict{Edge{T},Int}
end

function MolGraph{T,V,E}(
        g::SimpleGraph, vprop_map::Dict, eprop_map::Dict;
        gprop_map::Dict=Dict(), config_map::Dict=Dict(), kwargs...
        ) where {T,V,E}
    (nv(g) > length(vprop_map)
        && error("Mismatch in the number of nodes and node properties"))
    (ne(g) == length(eprop_map)
        || error("Mismatch in the number of edges and edge properties"))
    # expand fadjlist for vprops of isolated nodes
    for _ in nv(g):(length(vprop_map) - 1)
        push!(g.fadjlist, T[])
    end
    default_config = Dict(
        :initialized => false,
        :has_updates => true,
        :updater => mol -> (update_edge_rank!(mol); clear_caches!(mol); set_state!(mol, :has_updates, false)),
        :on_init => mol -> set_state!(mol, :initialized, true),
        :caches => Dict{Symbol,Any}()
    )
    merge!(default_config, config_map)
    mol = MolGraph{T,V,E}(g, vprop_map, eprop_map, gprop_map, default_config, Dict{Edge{T},T}())
    default_config[:initialized] || dispatch!(mol, :on_init)
    default_config[:has_updates] && dispatch!(mol, :updater)
    return mol
end

MolGraph(
    g::SimpleGraph{T}, vprop_map::Dict{T,V}, eprop_map::Dict{Edge{T},E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(g, vprop_map, eprop_map; kwargs...)


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

dispatch!(mol, event) = mol.state[event](mol)

get_state(mol::MolGraph, sym::Symbol) = mol.state[sym]
has_state(mol::MolGraph, sym::Symbol) = haskey(mol.state, sym)
set_state!(mol::MolGraph, sym::Symbol, value) = begin mol.state[sym] = value end

get_cache(mol::MolGraph, sym::Symbol) = mol.state[:caches][sym]
has_cache(mol::MolGraph, sym::Symbol) = haskey(mol.state[:caches], sym)
set_cache!(mol::MolGraph, sym::Symbol, value) = begin mol.state[:caches][sym] = value end
clear_caches!(mol::MolGraph) = empty!(mol.state[:caches])

function update_edge_rank!(mol::MolGraph)
    empty!(mol.edge_rank)
    for (i, e) in enumerate(edges(mol.graph))
        mol.edge_rank[e] = i
    end
end


# Edit graph

function add_u_edge!(mol::MolGraph{T,V,E}, e::Edge, prop::E) where {T,V,E}
    # Can be directly called if src < dst is guaranteed.
    add_edge!(mol.graph, e) || return false
    mol.eprops[e] = prop
    mol.state[:has_updates] = true
    return true
end


function add_vertex!(mol::MolGraph{T,V,E}, prop::V) where {T,V,E}
    add_vertex!(mol.graph) || return false
    mol.vprops[nv(mol.graph)] = prop
    mol.state[:has_updates] = true
    return true
end


function rem_u_edge!(mol::MolGraph, e::Edge)
    # Can be directly called if src < dst is guaranteed.
    rem_edge!(mol.graph, e) || return false
    delete!(mol.eprops, e)
    mol.state[:has_updates] = true
    return true
end


function rem_vertex!(mol::MolGraph, v::Integer)
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
    mol.state[:has_updates] = true
    return true
end


function rem_vertices!(mol::MolGraph{T,V,E}, vs::Vector{T}) where {T,V,E}
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
    mol.state[:has_updates] = true
    return vmap
end


function _induced_subgraph(mol::T, vlist_or_elist) where {T<:MolGraph}
    # In many cases, induced_subgraph(mol.graph, vlist_or_elist) is sufficient
    subg, vmap = induced_subgraph(mol.graph, vlist_or_elist)
    vps = Dict(v => mol.vprops[vmap[v]] for v in vertices(subg))
    eps = Dict(e => mol.eprops[u_edge(mol, vmap[src(e)], vmap[dst(e)])]
        for e in edges(subg))
    newgp = remap_gprops(mol, Dict(v => i for (i, v) in enumerate(vmap)))
    mol.state[:has_updates] = true
    return T(subg, vps, eps, gprop_map=newgp, config_map=mol.state), vmap
end
