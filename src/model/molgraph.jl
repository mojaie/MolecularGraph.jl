#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolGraph, SDFMolGraph, SMILESMolGraph, EditableMolGraph, Reaction,
    to_dict, to_json, setcache!,
    vproptype, eproptype,
    props, get_prop, eprop_iter, eprop_dict


struct MolGraph{T,A,B} <: OrderedMolGraph{T}
    graph::SimpleGraph{T}
    vprops::Vector{A}
    eprops::Vector{B}
    gprops::Dict{Symbol,Any}

    function MolGraph{T,A,B}(
            g::SimpleGraph{T}=SimpleGraph{T}(), vprops=A[], eprops=B[], gprops=Dict()) where {T,A,B}
        nv(g) > length(vprops) && throw(ErrorException("Mismatch in the number of nodes and node properties"))
        ne(g) == length(eprops) || throw(ErrorException("Mismatch in the number of edges and edge properties"))
        # expand fadjlist for vprops of isolated nodes
        for _ in nv(g):(length(vprops) - 1)
            push!(g.fadjlist, Vector{T}())
        end
        new(g, vprops, eprops, gprops)
    end
end

struct EditableMolGraph{T,A,B} <: AbstractMolGraph{T}
    graph::SimpleGraph{T}
    vprops::Dict{T,A}
    eprops::Dict{Edge{T},B}
    gprops::Dict{Symbol,Any}
end

SDFMolGraph = MolGraph{Int,SDFAtom,SDFBond}
SMILESMolGraph = MolGraph{Int,SMILESAtom,SMILESBond}

MolGraph(g=SimpleGraph{Int}(), vprops=Dict{Symbol,Any}[], eprops=Dict{Symbol,Any}[], gprops=Dict()
    ) = MolGraph{eltype(g),eltype(vprops),eltype(eprops)}(g, vprops, eprops, gprops)

function MolGraph{T,A,B}(edge_list::Vector{Edge{T}}, vprops=A[], eprops=B[], gprops=Dict()) where {T,A,B}
    # reorder edge properties
    mapping = Dict(e => eprops[i] for (i, e) in enumerate(edge_list))
    g = SimpleGraph(edge_list)
    new_eprops = B[]
    for e in edges(g)
        push!(new_eprops, mapping[e])
    end
    return MolGraph{T,A,B}(g, vprops, new_eprops, gprops)
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

"""
    MolGraph(mol::SubstructView{MolGraph}) -> GraphMol

Generate a new `GraphMol` from a substructure view.
"""

"""
function MolGraph(view::SubstructView{MolGraph})
    mol = zero(view.mol)
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(vertices(view))
        nmap[n] = i
        push!(mol.vprops, get_prop(view, n))
    end
    for e in edges(view)
        s = nmap[e.src]
        d = nmap[e.dst]
        add_edge!(mol.graph, s, d)
        push!(mol.eprops, get_prop(view, e))
    end
    return mol
end
"""

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

Base.eltype(::Type{MolGraph{T,A,B}}) where {T,A,B} = T  # not implemented in Graph.jl interface ?
edgetype(::Type{MolGraph{T,A,B}}) where {T,A,B} = Edge{T}

Base.zero(::Type{MolGraph{T,A,B}}) where {T,A,B} = MolGraph{T,A,B}()
Graphs.is_directed(::Type{MolGraph{T,A,B}}) where {T,A,B} = false

vproptype(::Type{MolGraph{T,A,B}}) where {T,A,B} = A
vproptype(mol::MolGraph) = vproptype(typeof(mol))
eproptype(::Type{MolGraph{T,A,B}}) where {T,A,B} = B
eproptype(mol::MolGraph) = eproptype(typeof(mol))

props(mol::MolGraph) = mol.gprops
props(mol::MolGraph, v::Integer) = mol.vprops[v]
props(mol::MolGraph, u::Integer, v::Integer) = mol.eprops[edge_rank(mol.graph, u, v)]
props(mol::MolGraph, e::Edge) = props(mol, src(e), dst(e))
get_prop(mol::MolGraph, prop::Symbol) = mol.gprops[prop]
get_prop(mol::MolGraph, v::Integer, prop::Symbol) = props(mol, v)[prop]
get_prop(mol::MolGraph, u::Integer, v::Integer, prop::Symbol) = props(mol, u, v)[prop]
get_prop(mol::MolGraph, e::Edge, prop::Symbol) = props(mol, e)[prop]
eprop_iter(mol::MolGraph) = zip(edges(mol), mol.eprops)
eprop_dict(mol::MolGraph) = Dict(eprop_iter(mol))

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


EditableMolGraph() = EditableMolGraph(
    SimpleGraph(), Dict{Int,Dict{Symbol,Any}}(), Dict{Edge{Int},Dict{Symbol,Any}}(), Dict())
EditableMolGraph{A,B}() where {A,B
    } = EditableMolGraph{A,B}(SimpleGraph(), Dict(), Dict(), Dict())

EditableMolGraph(g::SimpleGraph{Int}, vprops, eprops, gprops
    ) = EditableMolGraph{valtype(vprops),valtype(eprops)}(g, vprops, eprops, gprops)
EditableMolGraph(g::SimpleGraph{Int}, vprops, eprops
    ) = EditableMolGraph(g::SimpleGraph{Int}, vprops, eprops, Dict())
EditableMolGraph(edge_list::Vector{Edge}, vprops, eprops, gprops
    ) = EditableMolGraph(SimpleGraph(edge_list), vprops, eprops, gprops)
EditableMolGraph(edge_list::Vector{Edge}, vprops, eprops
    ) = EditableMolGraph(SimpleGraph(edge_list), vprops, eprops, Dict())

function EditableMolGraph(mol::MolGraph)
    vps = Dict{Int,valtype(mol.vprops)}()
    eps = Dict{Edge{Int},valtype(mol.eprops)}()
    for v in vertices(mol)
        vps[v] = mol.vprops[v]
    end
    for e in edges(mol)
        eps[e] = mol.vprops[e]
    end
    return EditableMolGraph(mol.graph, vps, eps, mol.gprops)
end

Base.zero(::Type{EditableMolGraph{T,A,B}}) where {T,A,B} = EditableMolGraph{T,A,B}()
Graphs.is_directed(::Type{EditableMolGraph{T,A,B}}) where {T,A,B} = false

vproptype(mol::EditableMolGraph) = valtype(mol.vprops)
eproptype(mol::EditableMolGraph) = valtype(mol.eprops)

props(mol::EditableMolGraph) = mol.gprops
props(mol::EditableMolGraph, v::Integer) = mol.vprops[v]
props(mol::EditableMolGraph, e::Edge) = mol.eprops[e]
props(mol::EditableMolGraph, u::Integer, v::Integer) = mol.eprops[undirectededge(u, v)]
get_prop(mol::EditableMolGraph, prop::Symbol) = mol.gprops[prop]
get_prop(mol::EditableMolGraph, v::Integer, prop::Symbol) = props(mol, v)[prop]
get_prop(mol::EditableMolGraph, e::Edge, prop::Symbol) = props(mol, e)[prop]
get_prop(mol::EditableMolGraph, u::Integer, v::Integer, prop::Symbol) = props(mol, u, v)[prop]
eprop_iter(mol::EditableMolGraph) = [mol.eprops[e] for e in edges(mol)]
eprop_dict(mol::EditableMolGraph) = mol.eprops


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


"""
    setcache!(graph::, key; kwargs...)

Set calculated property caches.
"""
function setcache!(mol::MolGraph, key; kwargs...)
    mol.gprops[:cache][key] = getfield(MolecularGraph, key)(mol; kwargs...)
end


function Base.getindex(graph::MolGraph, sym::Symbol)
    if haskey(graph.gprops[:cache], sym)
        return graph.gprops[:cache][sym]
    else
        return eval(Expr(:call, sym, graph))
    end
end

"""
function Base.getindex(view::SubstructView, sym::Symbol)
    if haskey(view.mol.gprops[:cache], sym)
        return view.mol.gprops[:cache][sym]
    else
        return eval(Expr(:call, sym, view.mol))
    end
end
"""
