#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    GraphMol, GraphReaction,
    graphmol, remapnodes, todict, tojson,
    SDFile, SMILES,
    getatom, getbond, hasbond,
    setatom!, setbond!, addatom!, addbond!,
    atomcount, bondcount

using JSON
import Unmarshal


struct GraphMol{A<:Atom,B<:Bond} <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{A}
    edgeattrs::Vector{B}
    cache::Dict{Symbol,Any}
    attributes::Dict{Symbol,Any}
end

struct GraphReaction
    reactants::Vector{GraphMol}
    products::Vector{GraphMol}
end

"""
    graphmol() -> GraphMol

Generate empty `GraphMol`.
"""
graphmol(::Type{A}, ::Type{B}
    ) where {A<:Atom,B<:Bond} = GraphMol{A,B}([], [], [], [], Dict(), Dict())

"""
    graphmol(atoms::Vector{Atom}, bonds::Vector{Bond}) -> GraphMol

Generate `GraphMol` that has the given atom objects and edge objects.
"""
function graphmol(
        edges, atoms::Vector{A}, bonds::Vector{B}) where {A<:Atom,B<:Bond}
    nbrmap = [Dict{Int,Int}() for i in 1:length(atoms)]
    edges = collect(edges)
    for (i, (u, v)) in enumerate(edges)
        nbrmap[u][i] = v
        nbrmap[v][i] = u
    end
    return GraphMol(
        nbrmap, edges, atoms, bonds, Dict{Symbol,Any}(), Dict{Symbol,Any}())
end

"""
    graphmol(mol::GraphMol; clearcache=true) -> GraphMol

Copy `GraphMol`. Cached properties will be removed.
"""
function graphmol(graph::GraphMol{A,B}; clearcache=true
        ) where {A<:Atom,B<:Bond}
    newg = graphmol(A, B)
    for nbr in graph.neighbormap
        push!(newg.neighbormap, copy(nbr))
    end
    append!(newg.edges, graph.edges)
    append!(newg.nodeattrs, graph.nodeattrs)
    append!(newg.edgeattrs, graph.edgeattrs)
    clearcache || merge!(newg.cache, graph.cache)
    merge!(newg.attributes, graph.attributes)
    return newg
end


Graph.clone(mol::GraphMol) = graphmol(mol, clearcache=false)


"""
    graphmol(mol::SubgraphView{GraphMol}) -> GraphMol

Generate a new `GraphMol` from a substructure view.

Graph property caches and attributes are not inherited.
"""
function graphmol(view::SubgraphView)
    newg = graphmol(nodeattrtype(view), edgeattrtype(view))
    nkeys = sort(collect(nodeset(view)))
    ekeys = sort(collect(edgeset(view)))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        push!(newg.nodeattrs, nodeattr(view, n))
        push!(newg.neighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        (oldu, oldv) = getedge(view, e)
        u = nmap[oldu]
        v = nmap[oldv]
        push!(newg.edges, (u, v))
        push!(newg.edgeattrs, edgeattr(view, e))
        newg.neighbormap[u][i] = v
        newg.neighbormap[v][i] = u
    end
    return newg
end


function remapnodes(graph::GraphMol, mapping::Dict{Int,Int})
    newg = graphmol(graph)
    for i in 1:nodecount(graph)
        newg.nodeattrs[mapping[i]] = graph.nodeattrs[i]
        newg.neighbormap[mapping[i]] = Dict(
            k => mapping[v] for (k, v) in graph.neighbormap[i])
    end
    for i in 1:edgecount(graph)
        u, v = graph.edges[i]
        newg.edges[i] = (mapping[u], mapping[v])
        newg.edgeattrs[i] = graph.edgeattrs[i]
    end
    merge!(newg.attributes, graph.attributes)
    return newg
end


"""
    todict(graph::GraphMol) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
todict(graph::GraphMol) = Dict(
    "edges" => graph.edges,
    "nodetype" => string(nodeattrtype(graph)),
    "edgetype" => string(edgeattrtype(graph)),
    "nodeattrs" => [todict(atom) for atom in graph.nodeattrs],
    "edgeattrs" => [todict(bond) for bond in graph.edgeattrs],
    "cache" => Dict(
        string(k) => (string(typeof(v)), v) for (k, v) in graph.cache),
    "attributes" => Dict(string(k) => v for (k, v) in graph.attributes)
)

tojson(graph::GraphMol) = JSON.json(todict(graph))



# Note: Unmarshal patch for Set
function Unmarshal.unmarshal(::Type{Set{E}},
    parsedJson::Union{Vector, AbstractArray},
    verbose::Bool=false, verboseLvl::Int=0) where E
    if verbose
        prettyPrint(verboseLvl, "Set{$E}")
        verboseLvl += 1
    end
    Set(Unmarshal.unmarshal(
        E, field, verbose, verboseLvl) for field in parsedJson)
end



function graphmol(data::Dict{String,Any})
    atomtype = eval(Meta.parse(data["nodetype"]))
    bondtype = eval(Meta.parse(data["edgetype"]))
    mol = graphmol(atomtype, bondtype)
    for a in data["nodeattrs"]
        push!(mol.nodeattrs, atomtype(a))
        push!(mol.neighbormap, Dict())
    end
    for (i, (u, v)) in enumerate(data["edges"])
        mol.neighbormap[u][i] = v
        mol.neighbormap[v][i] = u
        push!(mol.edgeattrs, bondtype(data["edgeattrs"][i]))
        push!(mol.edges, (u, v))
    end
    for (k, (vt, v)) in data["cache"]
        mol.cache[Symbol(k)] = Unmarshal.unmarshal(eval(Meta.parse(vt)), v)
    end
    for (k, v) in data["attributes"]
        mol.attributes[Symbol(k)] = v
    end
    return mol
end

graphmol(json::String) = graphmol(JSON.parse(json))


"""
    setcache!(graph::, key; kwargs...)

Set calculated property caches.
"""
function Graph.setcache!(mol::GraphMol, key; kwargs...)
    mol.cache[key] = getfield(MolecularGraph, key)(mol; kwargs...)
end


# Aliases

SDFile = GraphMol{SDFileAtom,SDFileBond}
SMILES = GraphMol{SmilesAtom,SmilesBond}


function Base.getindex(graph::GraphMol, sym::Symbol)
    if haskey(graph.cache, sym)
        return graph.cache[sym]
    else
        return eval(Expr(:call, sym, graph))
    end
end

function Base.getindex(view::SubgraphView, sym::Symbol)
    if haskey(view.graph.cache, sym)
        return view.graph.cache[sym]
    else
        return eval(Expr(:call, sym, view.graph))
    end
end
