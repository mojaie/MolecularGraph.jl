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

function MolGraph{T,V,E}(g::SimpleGraph,
        vprop_map::Dict, eprop_map::Dict, gprop_map::Dict, config_map::Dict) where {T,V,E}
    (nv(g) > length(vprop_map)
        && error("Mismatch in the number of nodes and node properties"))
    (ne(g) == length(eprop_map)
        || error("Mismatch in the number of edges and edge properties"))
    # expand fadjlist for vprops of isolated nodes
    for _ in nv(g):(length(vprop_map) - 1)
        push!(g.fadjlist, T[])
    end
    # edge_rank mapping
    er = Dict{Edge{T},T}()
    for (i, e) in enumerate(edges(g))
        er[e] = i
    end
    default_config = Dict(
        :has_updates => false,
        :on_update => mol -> (update_edge_rank!(mol); mol.state[:has_updates] = false)
    )
    merge!(default_config, config_map)
    mol = MolGraph{T,V,E}(g, vprop_map, eprop_map, gprop_map, default_config, er)
    dispatch!(mol, :on_update)
    return mol
end

MolGraph(g::SimpleGraph{T}, 
    vprop_map::Dict{T,V}, eprop_map::Dict{Edge{T},E}, gprop_map::Dict=Dict(), config_map::Dict=Dict()
) where {T,V,E} = MolGraph{T,V,E}(g, vprop_map, eprop_map, gprop_map, config_map)


# from empty set

MolGraph{T,V,E}() where {T,V,E} = MolGraph(SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}(), Dict())
MolGraph() = MolGraph{Int,Any,Any}()


# from edge and property list (sdftomol, smilestomol interface)

function MolGraph{T,V,E}(edge_list::Vector{Edge{T}},
        vprop_list::Vector{V}, eprop_list::Vector{E}, gprop_map::Dict=Dict(), config_map::Dict=Dict()
        ) where {T,V,E}
    g = SimpleGraph(edge_list)
    vps = Dict(i => v for (i, v) in enumerate(vprop_list))
    eps = Dict(e => eprop_list[i] for (i, e) in enumerate(edge_list))  # eprop_list in edge_list order
    return MolGraph{T,V,E}(g, vps, eps, gprop_map, config_map)
end

MolGraph(edge_list::Vector{Edge{T}},
    vprop_list::Vector{V}, eprop_list::Vector{E}, gprop_map::Dict=Dict(), config_map::Dict=Dict()
) where {T,V,E} = MolGraph{T,V,E}(edge_list, vprop_list, eprop_list, gprop_map, config_map)


# from dict (deserialize)

function MolGraph(data::Dict)
    eltype = eval(Meta.parse(data["eltype"]))
    vproptype = eval(Meta.parse(data["vproptype"]))
    eproptype = eval(Meta.parse(data["eproptype"]))
    g = SimpleGraph([Edge(e...) for e in data["graph"]])
    vps = Dict(i => vproptype(vp) for (i, vp) in enumerate(data["vprops"]))
    eps = Dict(e => eproptype(ep) for (e, ep) in zip(edges(g), data["eprops"]))
    gps = Dict(Symbol(k) => v for (k, v) in data["gprops"])
    return MolGraph(g, vps, eps, gps)
end

MolGraph(json::String) = MolGraph(JSON.parse(json))


# MolGraph type aliases

const SDFMolGraph = MolGraph{Int,SDFAtom,SDFBond}
const SMILESMolGraph = MolGraph{Int,SMILESAtom,SMILESBond}


to_dict(mol::MolGraph) = Dict(
    "eltype" => string(eltype(mol)),
    "vproptype" => string(vproptype(mol)),
    "eproptype" => string(eproptype(mol)),
    "graph" => [[src(e), dst(e)] for e in edges(mol)],
    "vprops" => [to_dict(props(mol, i)) for i in vertices(mol)],
    "eprops" => [to_dict(props(mol, e)) for e in edges(mol)],
    "gprops" => Dict(string(k) => v for (k, v) in mol.gprops)
)

dispatch!(mol, event) = mol.state[event](mol)
# descriptors(mol::MolGraph) = mol.descriptors
get_state(mol::MolGraph, sym::Symbol) = mol.state[sym]
has_state(mol::MolGraph, sym::Symbol) = haskey(mol.state, sym)
init_node_descriptor(::Type{T}, mol::MolGraph) where T = Vector{T}(undef, nv(mol))
init_node_descriptor(::Type{T}, mol::MolGraph) where T >: Nothing = Vector{T}(nothing, nv(mol))
init_node_descriptor(::Type{T}, mol::MolGraph) where T <: Number = fill(zero(T), nv(mol))
init_node_descriptor(::Type{Vector{T}}, mol::MolGraph) where T = [T[] for _ in vertices(mol)]
init_edge_descriptor(::Type{T}, mol::MolGraph) where T = Vector{T}(undef, ne(mol))
init_edge_descriptor(::Type{T}, mol::MolGraph) where T >: Nothing = Vector{T}(nothing, ne(mol))
init_edge_descriptor(::Type{T}, mol::MolGraph) where T <: Number = fill(zero(T), ne(mol))
init_edge_descriptor(::Type{Vector{T}}, mol::MolGraph) where T = [T[] for _ in edges(mol)]

function set_state!(mol::MolGraph, sym::Symbol, value)
    mol.state[sym] = value
    return value
end


function update_edge_rank!(mol::MolGraph)
    empty!(mol.edge_rank)
    for (i, e) in enumerate(edges(mol.graph))
        mol.edge_rank[e] = i
    end
end


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
    mol.vprops[v] = mol.vprops[nv_]
    delete!(mol.vprops, nv_)
    for e in edges_
        (src(e) == v || dst(e) == v) && delete!(mol.eprops, e)
        src(e) == nv_ && (mol.eprops[undirectededge(mol, v, dst(e))] = mol.eprops[e]; delete!(mol.eprops, e))
        dst(e) == nv_ && (mol.eprops[undirectededge(mol, src(e), v)] = mol.eprops[e]; delete!(mol.eprops, e))
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
        delete!(mol.vprops, v)
    end
    for e in edges(mol)
        _e = undirectededge(T, vmap[src(e)], vmap[dst(e)])
        e == _e && continue
        mol.eprops[e] = mol.eprops[_e]
        delete!(mol.eprops, _e)
    end
    mol.state[:has_updates] = true
    return vmap
end


function _induced_subgraph(mol::T, vlist_or_elist) where {T<:MolGraph}
    # In many cases, induced_subgraph(mol.graph, vlist_or_elist) is sufficient
    subg, vmap = induced_subgraph(mol.graph, vlist_or_elist)
    vps = Dict(v => mol.vprops[vmap[v]] for v in vertices(subg))
    eps = Dict(e => mol.eprops[undirectededge(mol, vmap[src(e)], vmap[dst(e)])]
        for e in edges(subg))
    return T(subg, vps, eps, mol.gprops, mol.state), vmap
end
