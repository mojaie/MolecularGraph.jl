#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    GraphMol, QueryMol,
    graphmol, remapnodes, todict, tojson,
    querymol, 
    SDFile, SMILES, SMARTS,
    getatom, getbond, hasbond,
    setatom!, setbond!, addatom!, addbond!,
    atomcount, bondcount

import JSON

struct GraphMol{A<:Atom,B<:Bond} <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{A}
    edgeattrs::Vector{B}
    cache::Dict{Symbol,Any}
    attributes::Dict{Symbol,Any}
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
        "cache" => Dict(string(k) => v for (k, v) in graph.cache),
        "attributes" => Dict(string(k) => v for (k, v) in graph.attributes)
)

tojson(graph::GraphMol) = JSON.json(todict(graph))


function graphmol(data::Dict{String,Any})
    atomtype = eval(Meta.parse(data["nodetype"]))
    bondtype = eval(Meta.parse(data["edgetype"]))
    mol = graphmol(atomtype, bondtype)
    for a in data["nodeattrs"]
        push!(mol.nodeattrs, atomtype(a))
        push!(mol.neighbormap, Dict())
    end
    for (i, (u, v)) in enumerate(data["edges"])
        mol.neighbormap[u][v] = i
        mol.neighbormap[v][u] = i
        push!(mol.edgeattrs, bondtype(data["edgeattrs"][i]))
        push!(mol.edges, (u, v))
    end
    for (k, v) in data["cache"]
        mol.cache[Symbol(k)] = v
    end
    for (k, v) in data["attributes"]
        mol.attributes[Symbol(k)] = v
    end
    return mol
end

graphmol(json::String) = graphmol(JSON.parse(json))



struct QueryMol{A<:QueryAtom,B<:QueryBond} <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{A}
    edgeattrs::Vector{B}
    cache::Dict{Symbol,Any}
    attributes::Dict{Symbol,Any}
    connectivity::Vector{Vector{Int}}
end

"""
    querymol() -> QueryMol

Generate empty `QueryMol`.
"""
querymol(
    ::Type{A}, ::Type{B}
) where {A<:QueryAtom,B<:QueryBond} = QueryMol{A,B}(
    [], [], [], [], Dict(), Dict(), [])

"""
    querymol(atoms::Vector{Atom}, bonds::Vector{Bond}) -> GraphMol

Generate `QueryMol` that has the given atom objects and edge objects.
"""
function querymol(edges, atoms::Vector{A}, bonds::Vector{B},
        connectivity::Vector{Vector{Int}}) where {A<:QueryAtom,B<:QueryBond}
    nbrmap = [Dict{Int,Int}() for i in 1:length(atoms)]
    edges = collect(edges)
    for (i, (u, v)) in enumerate(edges)
        nbrmap[u][i] = v
        nbrmap[v][i] = u
    end
    return QueryMol(
        nbrmap, edges, atoms, bonds,
        Dict{Symbol,Any}(), Dict{Symbol,Any}(), connectivity)
end


# Aliases

SDFile = GraphMol{SDFileAtom,SDFileBond}
SMILES = GraphMol{SmilesAtom,SmilesBond}
SMARTS = QueryMol{SmartsAtom,SmartsBond}


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



Molecule = Union{GraphMol,QueryMol}
IAtom = Union{Atom,QueryAtom}
IBond = Union{Bond,QueryBond}

# Aliases

getatom(mol::Molecule, i::Int) = nodeattr(mol, i)
getbond(mol::Molecule, i::Int) = edgeattr(mol, i)
getbond(mol::Molecule, u::Int, v::Int) = edgeattr(mol, u, v)
hasbond(mol::Molecule, u::Int, v::Int) = hasedge(mol, u, v)

setatom!(mol::Molecule, i::Int, attr::IAtom) = setnodeattr!(mol, i, attr)
setbond!(mol::Molecule, i::Int, attr::IBond) = setedgeattr!(mol, i, attr)
addatom!(mol::Molecule, attr::IAtom) = addnode!(mol, attr)
addbond!(mol::Molecule, u::Int, v::Int, attr::IBond) = addedge!(mol, u, v, attr)

atomcount(mol::Molecule) = nodecount(mol)
bondcount(mol::Molecule) = edgecount(mol)
