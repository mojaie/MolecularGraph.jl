#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    MolGraphGen{T,V,E} <: SimpleMolGraph{T,V,E}

MolGraph for manual moleculargraph generation.

This does not support descriptor array but enables faster addition and removal
of vertices and edges.
"""
struct MolGraphGen{T,V,E} <: SimpleMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::Dict{Symbol,Any}

    function MolGraphGen{T,V,E}(g::SimpleGraph,
            vprop_map::Dict, eprop_map::Dict, gprop_map::Dict) where {T,V,E}
        (nv(g) > length(vprop_map)
            && throw(ErrorException("Mismatch in the number of nodes and node properties")))
        (ne(g) == length(eprop_map)
            || throw(ErrorException("Mismatch in the number of edges and edge properties")))
        # expand fadjlist for vprops of isolated nodes
        for _ in nv(g):(length(vprop_map) - 1)
            push!(g.fadjlist, T[])
        end
        new(g, vprop_map, eprop_map, gprop_map)
    end
end

MolGraphGen{T,V,E}() where {T,V,E} = MolGraphGen{T,V,E}(
    SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}(), Dict())
MolGraphGen() = MolGraphGen{Int,Any,Any}()

function MolGraphGen{T,V,E}(g::SimpleGraph,
        vprop_list::Vector, eprop_list::Vector, gprop_map::Dict=Dict()
        ) where {T,V,E}
    vps = Dict(v => vprop_list[v] for v in vertices(g))
    eps = Dict(e => eprop_list[i] for (i, e) in enumerate(edges(g)))
    return MolGraphGen{T,V,E}(g, vps, eps, gprop_map)
end

MolGraphGen(g::SimpleGraph{T}, 
    vprop_list::Vector{V}, eprop_list::Vector{E}, gprop_map::Dict=Dict()
) where {T,V,E} = MolGraphGen{T,V,E}(g, vprop_list, eprop_list, gprop_map)

function MolGraphGen{T,V,E}(edge_list::Vector{Edge{T}},
        vprop_list::Vector{V}, eprop_list::Vector{E}, gprop_map::Dict=Dict()
        ) where {T,V,E}
    # sdftomol, smilestomol interface
    g = SimpleGraph(edge_list)
    vps = Dict(v => vprop_list[v] for v in vertices(g))
    eps = Dict(e => eprop_list[i] for (i, e) in enumerate(edge_list))
    return MolGraphGen{T,V,E}(g, vps, eps, gprop_map)
end

MolGraphGen(edge_list::Vector{Edge{T}},
    vprop_list::Vector{V}, eprop_list::Vector{E}, gprop_map::Dict=Dict()
) where {T,V,E} = MolGraphGen{T,V,E}(edge_list, vprop_list, eprop_list, gprop_map)


function MolGraphGen(mol::MolGraph)
    vps = Dict(v => mol.vprops[v] for v in vertices(mol))
    eps = Dict(e => mol.eprops[i] for (i, e) in enumerate(edges(mol)))
    return MolGraphGen(mol.graph, vps, eps, mol.gprops)
end

function MolGraph(mol::MolGraphGen)
    vps = [mol.vprops[v] for v in vertices(mol)]
    eps = [mol.eprops[e] for e in edges(mol)]
    return MolGraph(mol.graph, vps, eps, mol.gprops)
end


function MolGraphGen(data::Dict)
    eltype = eval(Meta.parse(data["eltype"]))
    vproptype = eval(Meta.parse(data["vproptype"]))
    eproptype = eval(Meta.parse(data["eproptype"]))
    g = SimpleGraph([Edge(e...) for e in data["graph"]])
    vps = [vproptype(vp) for vp in data["vprops"]]
    eps = [eproptype(ep) for ep in data["eprops"]]
    gps = Dict(Symbol(k) => v for (k, v) in data["gprops"])
    return MolGraphGen(g, vps, eps, gps)
end

MolGraphGen(json::String) = MolGraphGen(JSON.parse(json))


to_dict(mol::MolGraphGen) = to_dict(MolGraph(mol))
to_json(mol::MolGraphGen) = JSON.json(to_dict(mol))

props(mol::MolGraphGen, e::Edge) = mol.eprops[e]
edge_rank(mol::MolGraphGen, e::Edge) = e
has_descriptor(mol::MolGraphGen, desc::Symbol) = false


function add_u_edge!(mol::MolGraphGen{T,V,E}, e::Edge, prop::E) where {T,V,E}
    # Can be directly called if src < dst is guaranteed.
    add_edge!(mol.graph, e) || throw(ErrorException("failed to add edge $(e)"))
    mol.eprops[e] = prop
    # TODO: stereochemistry, empty descriptors
    return true
end


function add_u_edges!(mol::MolGraphGen, elist, plist)
    for (e, p) in zip(elist, plist)
        add_edge!(mol, e, p)
    end
    return true
end


function add_vertex!(mol::MolGraphGen{T,V,E}, prop::V) where {T,V,E}
    add_vertex!(mol.graph) || throw(ErrorException("failed to add vertex"))
    mol.vprops[nv(mol.graph)] = prop
    # TODO: stereochemistry, empty descriptors
    return true
end


function rem_u_edge!(mol::MolGraphGen, e::Edge)
    # Can be directly called if src < dst is guaranteed.
    rem_edge!(mol.graph, e) || throw(ErrorException("failed to remove edge $(e)"))
    delete!(mol.eprops, e)
    # TODO: stereochemistry, empty descriptors
    return true
end


function rem_edges!(mol::MolGraphGen, edges)
    for e in edges
        rem_edge!(mol, e)
    end
    return true
end


function rem_vertex!(mol::MolGraphGen, v::Integer)
    incs = [undirectededge(mol, v, nbr) for nbr in neighbors(mol, v)]
    rem_edges!(mol, incs)
    rem_vertex!(mol.graph, v) || throw(ErrorException("failed to remove vertex $(v)"))
    mol.vprops[v] = mol.vprops[nv(mol)]
    delete!(mol.vprops, nv(mol))
    # TODO: stereochemistry, empty descriptors
    return true
end


function rem_vertices!(mol::MolGraphGen{T,V,E}, vs::Vector{T}) where {T,V,E}
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
    # TODO: stereochemistry
    return vmap
end


function _induced_subgraph(mol::T, vlist_or_elist) where {T<:MolGraphGen}
    # In many cases, induced_subgraph(mol.graph, vlist_or_elist) is sufficient
    subg, vmap = induced_subgraph(mol.graph, vlist_or_elist)
    vps = Dict(v => mol.vprops[vmap[v]] for v in vertices(subg))
    eps = Dict(e => mol.eprops[undirectededge(mol, vmap[src(e)], vmap[dst(e)])]
        for e in edges(subg))
    # TODO: stereochemistry
    return T(subg, vps, eps, mol.gprops), vmap
end
