#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    MolGraph{T,V,E} <: SimpleMolGraph{T,V,E}

Basic molecular graph type.

Note that graph manipulation (add/remove vertices/edges) functions are not
efficient in `MolGraph`. For manual graph generation tasks, use `MolGraphGen` instead.
"""
struct MolGraph{T,V,E} <: SimpleMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Vector{V}
    eprops::Vector{E}
    gprops::Dict{Symbol,Any}
    descriptors::Dict{Symbol,Any}
    edge_rank::Dict{Edge{T},Int}

    function MolGraph{T,V,E}(g::SimpleGraph,
            vprop_list::Vector, eprop_list::Vector, gprop_map::Dict) where {T,V,E}
        (nv(g) > length(vprop_list)
            && throw(ErrorException("Mismatch in the number of nodes and node properties")))
        (ne(g) == length(eprop_list)
            || throw(ErrorException("Mismatch in the number of edges and edge properties")))
        # expand fadjlist for vprops of isolated nodes
        for _ in nv(g):(length(vprop_list) - 1)
            push!(g.fadjlist, T[])
        end
        # edge_rank mapping
        er = Dict{Edge{T},T}()
        for (i, e) in enumerate(edges(g))
            er[e] = i
        end
        new(g, vprop_list, eprop_list, gprop_map, Dict(), er)
    end
end

MolGraph(g::SimpleGraph{T}, 
    vprop_list::Vector{V}, eprop_list::Vector{E}, gprop_map::Dict=Dict()
) where {T,V,E} = MolGraph{T,V,E}(g, vprop_list, eprop_list, gprop_map)

MolGraph{T,V,E}() where {T,V,E} = MolGraph(SimpleGraph{T}(), V[], E[], Dict())
MolGraph() = MolGraph{Int,Any,Any}()

function MolGraph{T,V,E}(edge_list::Vector{Edge{T}},
        vprop_list::Vector{V}, eprop_list::Vector{E}, gprop_map::Dict=Dict()
        ) where {T,V,E}
    # reorder edge properties
    mapping = Dict(e => eprop_list[i] for (i, e) in enumerate(edge_list))
    g = SimpleGraph(edge_list)
    new_eprops = [mapping[e] for e in edges(g)]
    return MolGraph{T,V,E}(g, vprop_list, new_eprops, gprop_map)
end

MolGraph(edge_list::Vector{Edge{T}},
    vprop_list::Vector{V}, eprop_list::Vector{E}, gprop_map::Dict=Dict()
) where {T,V,E} = MolGraph{T,V,E}(edge_list, vprop_list, eprop_list, gprop_map)


function MolGraph(data::Dict)
    eltype = eval(Meta.parse(data["eltype"]))
    vproptype = eval(Meta.parse(data["vproptype"]))
    eproptype = eval(Meta.parse(data["eproptype"]))
    g = SimpleGraph([Edge(e...) for e in data["graph"]])
    vps = [vproptype(vp) for vp in data["vprops"]]
    eps = [eproptype(ep) for ep in data["eprops"]]
    gps = Dict(Symbol(k) => v for (k, v) in data["gprops"])
    return MolGraph(g, vps, eps, gps)
end

MolGraph(json::String) = MolGraph(JSON.parse(json))


# Aliases
SDFMolGraph = MolGraph{Int,SDFAtom,SDFBond}
SMILESMolGraph = MolGraph{Int,SMILESAtom,SMILESBond}


"""
    to_dict(mol::MolGraph) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
to_dict(mol::MolGraph) = Dict(
    "eltype" => string(eltype(mol)),
    "vproptype" => string(vproptype(mol)),
    "eproptype" => string(eproptype(mol)),
    "graph" => [[src(e), dst(e)] for e in edges(mol)],
    "vprops" => [to_dict(vp) for vp in mol.vprops],
    "eprops" => [to_dict(ep) for ep in mol.eprops],
    "gprops" => Dict(string(k) => v for (k, v) in mol.gprops)
)

to_json(mol::MolGraph) = JSON.json(to_dict(mol))

props(mol::MolGraph, e::Edge) = eprops(mol)[mol.edge_rank[e]]
edge_rank(mol::MolGraph, e::Edge) = mol.edge_rank[e]
descriptors(mol::MolGraph) = mol.descriptors
get_descriptor(mol::MolGraph, desc::Symbol) = descriptors(mol)[desc]
has_descriptor(mol::MolGraph, desc::Symbol) = haskey(descriptors(mol), desc)
init_node_descriptor(::Type{T}, mol::MolGraph) where T = Vector{T}(undef, nv(mol))
init_node_descriptor(::Type{T}, mol::MolGraph) where T >: Nothing = Vector{T}(nothing, nv(mol))
init_node_descriptor(::Type{T}, mol::MolGraph) where T <: Number = fill(zero(T), nv(mol))
init_node_descriptor(::Type{Vector{T}}, mol::MolGraph) where T = [T[] for _ in vertices(mol)]
init_edge_descriptor(::Type{T}, mol::MolGraph) where T = Vector{T}(undef, ne(mol))
init_edge_descriptor(::Type{T}, mol::MolGraph) where T >: Nothing = Vector{T}(nothing, ne(mol))
init_edge_descriptor(::Type{T}, mol::MolGraph) where T <: Number = fill(zero(T), ne(mol))
init_edge_descriptor(::Type{Vector{T}}, mol::MolGraph) where T = [T[] for _ in edges(mol)]

function set_descriptor!(mol::MolGraph, desc::Symbol, value)
    # delete!(mol.descriptors, desc)  # remove old descriptor
    mol.descriptors[desc] = value
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
    # Note: Not efficient. For frequent graph manipulation, use MolGraphGen instead.
    add_edge!(mol.graph, e) || throw(ErrorException("failed to add edge $(e)"))
    update_edge_rank!(mol)
    insert!(mol.eprops, mol.edge_rank[e], prop)
    # TODO: stereochemistry, empty descriptors
    return true
end


function add_u_edges!(mol::MolGraph, elist, plist)
    for e in elist
        add_edge!(mol.graph, e) || throw(ErrorException("failed to add edge $(e)"))
    end
    update_edge_rank!(mol)
    for (e, p) in zip(elist, plist)
        insert!(mol.eprops, mol.edge_rank[e], p)
    end
    return true
end


function add_vertex!(mol::MolGraph{T,V,E}, prop::V) where {T,V,E}
    add_vertex!(mol.graph) || throw(ErrorException("failed to add vertex"))
    push!(mol.vprops, prop)
    # TODO: stereochemistry, empty descriptors
    return true
end


function rem_u_edge!(mol::MolGraph, e::Edge)
    # Can be directly called if src < dst is guaranteed.
    # Not efficient. For frequent graph manipulation, use MolGraphGen instead.
    rem_edge!(mol.graph, e) || throw(ErrorException("failed to remove edge $(e)"))
    deleteat!(mol.eprops, mol.edge_rank[e])
    update_edge_rank!(mol)
    # TODO: stereochemistry, empty descriptors
    return true
end


function rem_edges!(mol::MolGraph, edges)
    for e in edges
        rem_edge!(mol.graph, e) || throw(ErrorException("failed to remove edge $(e)"))
    end
    deleteat!(mol.eprops, [mol.edge_rank[e] for e in edges])
    update_edge_rank!(mol)
    # TODO: stereochemistry, empty descriptors
    return true
end


function rem_vertex!(mol::MolGraph, v::Integer)
    # Not efficient. for frequent graph topology manipulation, use MolGraphGen instead.
    incs = [undirectededge(mol, v, nbr) for nbr in neighbors(mol, v)]
    rem_edges!(mol, incs)
    rem_vertex!(mol.graph, v) || throw(ErrorException("failed to remove vertex $(v)"))
    mol.vprops[v] = mol.vprops[length(mol.vprops)]
    deleteat!(mol.vprops, length(mol.vprops))
    update_edge_rank!(mol)
    # TODO: stereochemistry, empty descriptors
    return true
end


function rem_vertices!(mol::MolGraph{T,V,E}, vs::Vector{T}) where {T,V,E}
    # Not efficient. for frequent graph topology manipulation, use MolGraphGen instead.
    # TODO: if many vertices should be removed, induced_subgraph may be more efficient.
    vmap = rem_vertices!(mol.graph, vs)
    vps = mol.vprops[vmap]
    emap = [mol.edge_rank[undirectededge(mol, vmap[src(e)], vmap[dst(e)])]
        for e in edges(mol.graph)]
    eps = mol.eprops[emap]
    # TODO: compare empty!->append! and keep_order=true->deleteat!
    empty!(mol.vprops)
    empty!(mol.eprops)
    empty!(mol.descriptors)
    # TODO: stereochemistry
    append!(mol.vprops, vps)
    append!(mol.eprops, eps)
    update_edge_rank!(mol)
    return vmap
end


function _induced_subgraph(mol::T, vlist_or_elist) where {T<:MolGraph}
    # Note that this will return substructure without calculated descriptors.
    # In many cases, induced_subgraph(mol.graph, vlist_or_elist) is sufficient
    subg, vmap = induced_subgraph(mol.graph, vlist_or_elist)
    vps = mol.vprops[vmap]
    emap = [mol.edge_rank[undirectededge(mol, vmap[src(e)], vmap[dst(e)])]
        for e in edges(subg)]
    eps = mol.eprops[emap]
    # TODO: stereochemistry
    return T(subg, vps, eps, mol.gprops), vmap
end
