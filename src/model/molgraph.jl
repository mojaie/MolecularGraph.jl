#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

import Graphs:
    add_vertex!, add_edge!, rem_vertex!, rem_vertices!, rem_edge!

export
    MolGraph, SDFMolGraph, SMILESMolGraph, SMARTSMolGraph,
    EditableMolGraph, Reaction,
    to_dict, to_json, setcache!,
    descriptors, set_descriptor!, has_descriptor, get_descriptor


struct MolGraph{T,V,E} <: SimpleMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Vector{V}
    eprops::Vector{E}
    gprops::Dict{Symbol,Any}
    descriptors::Dict{Symbol,Any}
    edge_rank::Dict{Edge{T},T}

    function MolGraph{T,V,E}(
            g=SimpleGraph{T}(), vprops=V[], eprops=E[], gprops=Dict()) where {T,V,E}
        nv(g) > length(vprops) && throw(ErrorException("Mismatch in the number of nodes and node properties"))
        ne(g) == length(eprops) || throw(ErrorException("Mismatch in the number of edges and edge properties"))
        # expand fadjlist for vprops of isolated nodes
        for _ in nv(g):(length(vprops) - 1)
            push!(g.fadjlist, Vector{T}())
        end
        # edge_rank mapping
        er = Dict{Edge{T},T}()
        for (i, e) in enumerate(edges(g))
            er[e] = i
        end
        new(g, vprops, eprops, gprops, Dict(), er)
    end
end


struct EditableMolGraph{T,V,E} <: SimpleMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::Dict{Symbol,Any}

    function EditableMolGraph{T,V,E}(
            g=SimpleGraph{T}(), vprops=Dict{T,V}(), eprops=Dict{Edge{T},E}(), gprops=Dict()) where {T,V,E}
        nv(g) > length(vprops) && throw(ErrorException("Mismatch in the number of nodes and node properties"))
        ne(g) == length(eprops) || throw(ErrorException("Mismatch in the number of edges and edge properties"))
        # expand fadjlist for vprops of isolated nodes
        for _ in nv(g):(length(vprops) - 1)
            push!(g.fadjlist, Vector{T}())
        end
        new(g, vprops, eprops, gprops)
    end
end


SDFMolGraph = MolGraph{Int,SDFAtom,SDFBond}
SMILESMolGraph = MolGraph{Int,SMILESAtom,SMILESBond}
SMARTSMolGraph = MolGraph{Int,QueryTruthTable,QueryTruthTable}

MolGraph(g=SimpleGraph{Int}(), vprops=Dict{Symbol,Any}[], eprops=Dict{Symbol,Any}[], gprops=Dict()
    ) = MolGraph{eltype(g),eltype(vprops),eltype(eprops)}(g, vprops, eprops, gprops)

function MolGraph{T,V,E}(edge_list::AbstractVector{Edge{T}}, vprops=V[], eprops=E[], gprops=Dict()) where {T,V,E}
    # reorder edge properties
    mapping = Dict(e => eprops[i] for (i, e) in enumerate(edge_list))
    g = SimpleGraph(edge_list)
    new_eprops = E[]
    for e in edges(g)
        push!(new_eprops, mapping[e])
    end
    return MolGraph{T,V,E}(g, vprops, new_eprops, gprops)
end

MolGraph(edge_list::Vector, vprops=[], eprops=[], gprops=Dict()
    ) = MolGraph{eltype(eltype(edge_list)),eltype(vprops),eltype(eprops)}(edge_list, vprops, eprops, gprops)


function MolGraph(emol::EditableMolGraph)
    vps = []
    eps = []
    for v in vertices(mol)
        push!(vps, emol.vprops[v])
    end
    for e in edges(mol)
        push!(eps, emol.eprops[e])
    end
    return MolGraph(emol.graph, vps, eps, emol.gprops)
end


function MolGraph(data::Dict)
    vproptype = eval(Meta.parse(data["vproptype"]))
    eproptype = eval(Meta.parse(data["eproptype"]))
    g = SimpleGraph(Edge.(data["graph"]))
    vps = [vproptype(vp) for vp in data["vprops"]]
    eps = [eproptype(vp) for ep in data["eprops"]]
    gps = Dict(Symbol(k) => v for (k, v) in data["gprops"])
    return MolGraph(g, vps, eps, gps)
end

MolGraph(json::String) = MolGraph(JSON.parse(json))


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
    delete!(mol.descriptors, desc)  # remove old descriptor
    mol.descriptors[desc] = value
    return value
end

"""
    to_dict(mol::MolGraph) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
to_dict(mol::MolGraph) = Dict(
    "graph" => to_dict(mol.graph),
    "vproptype" => string(eltype(mol.vprops)),
    "eproptype" => string(eltype(mol.eprops)),
    "vprops" => [to_dict(vp) for vp in mol.vprops],
    "eprops" => [to_dict(ep) for vp in mol.eprops],
    "gprops" => Dict(string(k) => v for (k, v) in mol.gprops)
)

to_json(mol::MolGraph) = JSON.json(to_dict(mol))


EditableMolGraph(
    g=SimpleGraph{Int}(), vprops=Dict{Int,Dict{Symbol,Any}}(), eprops=Dict{Edge{Int},Dict{Symbol,Any}}(), gprops=Dict()
) = EditableMolGraph{eltype(g),valtype(vprops),valtype(eprops)}(g, vprops, eprops, gprops)

EditableMolGraph{T,V,E}(
    edge_list::AbstractVector{Edge{T}}, vprops=Dict{T,V}(),
    eprops=Dict{Edge{T},E}(), gprops=Dict()
) where {T,V,E} = EditableMolGraph{T,V,E}(SimpleGraph(edge_list), vprops, eprops, gprops)

EditableMolGraph(
    edge_list::AbstractVector{T}, vprops=Dict(), eprops=Dict(), gprops=Dict()
) where T = MolGraph{eltype(T),valtype(vprops),valtype(eprops)}(edge_list, vprops, eprops, gprops)

function EditableMolGraph{T,V,E}(
        edge_list::AbstractVector{Edge{T}}, vprops::AbstractVector{V},
        eprops::AbstractVector{E}, gprops=Dict()) where {T,V,E}
    # sdftomol, smilestomol interface
    # TODO: need edge sanitization? -> undirectededge(T, e)
    g = SimpleGraph(edge_list)
    vps = Dict(v => vprops[v] for v in vertices(g))
    eps = Dict(e => eprops[i] for (i, e) in enumerate(edge_list))  # reorder edge properties
    return EditableMolGraph{T,V,E}(g, vps, eps, gprops)
end


function EditableMolGraph(mol::MolGraph{T,V,E}) where {T,V,E}
    vps = Dict{T,V}()
    eps = Dict{Edge{T},E}()
    for v in vertices(mol)
        vps[v] = mol.vprops[v]
    end
    for (i, e) in enumerate(edges(mol))
        eps[e] = mol.eprops[i]
    end
    return EditableMolGraph(mol.graph, vps, eps, mol.gprops)
end


props(mol::EditableMolGraph, e::Edge) = mol.eprops[e]
edge_rank(mol::EditableMolGraph, e::Edge) = e
has_descriptor(mol::EditableMolGraph, desc::Symbol) = false
init_node_descriptor(::Type{D}, mol::EditableMolGraph{T,V,E}, length::Int
    ) where {D,T,V,E} = Dict{T,D}()
init_node_descriptor(::Type{D}, mol::EditableMolGraph{T,V,E}
    ) where {D>:Nothing,T,V,E} = Dict{T,D}(i => nothing for i in vertices(mol))
init_node_descriptor(::Type{D}, mol::EditableMolGraph{T,V,E}
    ) where {D<:Number,T,V,E} = Dict{T,D}(i => zero(D) for i in vertices(mol))
init_node_descriptor(::Type{Vector{D}}, mol::EditableMolGraph{T,V,E}
    ) where {D,T,V,E} = Dict{T,Vector{D}}(i => D[] for i in vertices(mol))
init_edge_descriptor(::Type{D}, mol::EditableMolGraph{T,V,E}
    ) where {D,T,V,E} = Dict{Edge{T},D}()
init_edge_descriptor(::Type{D}, mol::EditableMolGraph{T,V,E}
    ) where {D>:Nothing,T,V,E} = Dict{Edge{T},D}(e => nothing for e in edges(mol))
init_edge_descriptor(::Type{D}, mol::EditableMolGraph{T,V,E}
    ) where {D<:Number,T,V,E} = Dict{Edge{T},D}(e => zero(D) for e in edges(mol))
init_edge_descriptor(::Type{Vector{D}}, mol::EditableMolGraph{T,V,E}
    ) where {D,T,V,E} = Dict{T,Vector{D}}(e => D[] for e in edges(mol))

function add_edge!(mol::EditableMolGraph{T,V,E}, e::Edge, prop::E) where {T,V,E}
    add_edge!(mol.graph, e) || return false
    mol.eprops[e] = prop
    return true
end
add_edge!(mol::EditableMolGraph{T,V,E}, u::Integer, v::Integer, prop::E
    ) where {T,V,E} = add_edge!(mol, undirectededge(mol, u, v), prop)

function add_vertex!(mol::EditableMolGraph{T,V,E}, prop::V) where {T,V,E}
    add_vertex!(mol.graph) || return false
    mol.vprops[nv(mol.graph)] = prop
    return true
end

function rem_edge!(mol::EditableMolGraph, e::Edge)
    rem_edge!(mol.graph, e) || return false
    delete!(mol.eprops, e)
    return true
end
rem_edge!(mol::EditableMolGraph, u::Integer, v::Integer) = rem_edge!(mol, undirectededge(mol, u, v))

function rem_vertex!(mol::EditableMolGraph, v::Integer)
    for nbr in neighbors(mol, v)
        rem_edge!(mol, undirectededge(mol, v, nbr)) || return false
    end
    rem_vertex!(mol.graph, v) || return false
    mol.vprops[v] = mol.vprops[length(mol.vprops)]
    delete!(mol.vprops, length(mol.vprops))
    return true
end

function rem_vertices!(
        mol::EditableMolGraph{T,V,E}, vs::AbstractVector{<:Integer}) where {T,V,E}
    issubset(Set(vs), Set(vertices(mol))) || return false
    subg, vmap = induced_subgraph(mol.graph, vs)
    new_eprops = Dict{Edge{T},E}()
    for e in edges(subg)
        new_eprops[e] = mol.eprops[undirectededge(mol, vmap[src(e)], vmap[dst(e)])]
    end
    new_vprops = Dict{T,V}()
    for v in vertices(subg)
        new_vprops[v] = mol.vprops[vmap[v]]
    end
    mol.graph.ne = subg.ne
    mol.graph.fadjlist = subg.fadjlist
    empty!(mol.vprops)
    merge!(mol.vprops, new_vprops)
    empty!(mol.eprops)
    merge!(new_eprops, new_eprops)
    return true
end


struct Reaction{T} <: AbstractReaction{T}
    reactants::Vector{T}
    products::Vector{T}
    rprops::Dict{Symbol,Any}
end

Reaction{T}() where T = Reaction{T}([], [], Dict())

Base.eltype(::Type{Reaction{T}}) where T = T
Base.eltype(rxn::Reaction{T}) where T = T


""" TO BE REMOVED
function remapnodes(graph::GraphMol, mapping::Dict{Int,Int})
    newg = graphmol(graph)
    for i in 1:nv(graph)
        newg.nodeattrs[mapping[i]] = graph.nodeattrs[i]
        newg.neighbormap[mapping[i]] = Dict(
            k => mapping[v] for (k, v) in graph.neighbormap[i])
    end
    for i in 1:ne(graph)
        u, v = graph.edges[i]
        newg.edges[i] = (mapping[u], mapping[v])
        newg.edgeattrs[i] = graph.edgeattrs[i]
    end
    merge!(newg.attributes, graph.attributes)
    return newg
end


# Note: Unmarshal patch for Set
# TODO: get rid of sets. convert to vectors
function Unmarshal.unmarshal(::Type{Set{E}},
    parsedJson::Union{Vector, AbstractArray},
    verbose::Bool=false, verboseLvl::Int=0) where E
    if verbose
        prettyPrint(verboseLvl, "Set{\$E}")
        verboseLvl += 1
    end
    Set(Unmarshal.unmarshal(
        E, field, verbose, verboseLvl) for field in parsedJson)
end
"""
