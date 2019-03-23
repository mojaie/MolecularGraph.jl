#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MapMol, VectorMol, QueryMol,
    mapmol, vectormol, querymol,
    SDFile, SMILES, SMARTS


struct MapMol{A<:Atom,B<:Bond} <: MolGraph
    graph::MapGraph{A,B}
    attribute::Dict{Symbol,String}

    function MapMol{A,B}(graph::Graph) where {A<:Atom,B<:Bond}
        new(graph, Dict())
    end
end

MapMol{A,B}() where {A<:Atom,B<:Bond} = MapMol{A,B}(mapgraph(A,B))

"""
    mapmol(::Type{A}, ::Type{B}) where {A<:Atom,B<:Bond} -> MapMol{N,E}

Generate empty `MapMol` that has atoms and bonds with the given types.
"""
mapmol(::Type{A}, ::Type{B}) where {A<:Atom,B<:Bond} = MapMol{A,B}()

"""
    mapmol(atoms::Vector{A}, bonds::Vector{B}) -> MapMol{A,B}

Generate `MapMol` that has the given atom objects and edge objects.
"""
function mapmol(atoms::Vector{A}, bonds::Vector{B}) where {A<:Atom,B<:Bond}
    mol = MapMol{A,B}()
    for (i, atom) in enumerate(atoms)
        mol.graph.nodes[i] = atom
        mol.graph.neighbormap[i] = Dict()
    end
    for (i, bond) in enumerate(bonds)
        mol.graph.edges[i] = bond
        mol.graph.neighbormap[bond.u][bond.v] = i
        mol.graph.neighbormap[bond.v][bond.u] = i
    end
    return mol
end


"""
    mapmol(mol::MolGraph{A,B}; clone=false) -> MapMol{A,B}

Convert the given molecule into a new `MapMol`. See [`mapgraph`](@ref)
for the details.
"""
function mapmol(mol::MolGraph)
    A = nodetype(mol)
    B = edgetype(mol)
    newmol = MapMol{A,B}(mapgraph(mol))
    merge!(newmol.attribute, mol.attribute)
    return newmol
end



struct QueryMol{A<:QueryAtom,B<:QueryBond} <: MolGraph
    graph::MapGraph{A,B}
    connectivity::Array{Array{Int}}

    function QueryMol{A,B}() where {A<:QueryAtom,B<:QueryBond}
        new(mapgraph(A,B), [])
    end
end

"""
    querymol(::Type{A}, ::Type{B}
        ) where {A<:QueryAtom,B<:QueryBond} -> QueryMol{N,E}

Generate empty `QueryMol` that has atoms and bonds with the given types.
"""
querymol(::Type{A}, ::Type{B}
    ) where {A<:QueryAtom,B<:QueryBond} = QueryMol{A,B}()


struct VectorMol{A<:Atom,B<:Bond} <: MolGraph
    graph::VectorGraph{A,B}
    vector::Dict{Symbol,Vector}
    annotation::Dict{Symbol,Annotation}
    coords::Dict{Symbol,Coordinates}
    attribute::Dict{Symbol,String}

    function VectorMol{A,B}(graph::Graph) where {A<:Atom,B<:Bond}
        new(graph, Dict(), Dict(), Dict(), Dict())
    end
end

VectorMol{A,B}() where {A<:Atom,B<:Bond} = VectorMol{A,B}(vectorgraph(A,B))


"""
    vectormol(::Type{A}, ::Type{B}) where {A<:Atom,B<:Bond} -> VectorMol{A,B}

Generate empty `VectorMol` that has atoms and bonds with the given types.
"""
vectormol(::Type{A}, ::Type{B}) where {A<:Atom,B<:Bond} = VectorMol{A,B}()

"""
    vectormol(atoms::Vector{A}, bonds::Vector{B}) -> VectorMol{A,B}

Generate `VectorMol` that has the given atom objects and edge objects.
"""
function vectormol(nodes::Vector{A}, edges::Vector{B}) where {A<:Atom,B<:Bond}
    mol = VectorMol{A,B}()
    for (i, node) in enumerate(nodes)
        push!(mol.graph.nodes, node)
        push!(mol.graph.neighbormap, Dict())
    end
    for (i, edge) in enumerate(edges)
        push!(mol.graph.edges, edge)
        mol.graph.neighbormap[edge.u][edge.v] = i
        mol.graph.neighbormap[edge.v][edge.u] = i
    end
    return mol
end

"""
    vectormol(mol::MolGraph; clone=false) -> VectorMol

Convert the given molecule into a new `VectorMol`. See [`vectorgraph`](@ref)
for the details.
"""
function vectormol(mol::MolGraph)
    A = nodetype(mol)
    B = edgetype(mol)
    newmol = VectorMol{A,B}(vectorgraph(mol))
    merge!(newmol.attribute, mol.attribute)
    return newmol
end


# Aliases
# TODO: use traits


SDFile = MapMol{SDFileAtom,SDFileBond}
SMILES = MapMol{SmilesAtom,SmilesBond}
SMARTS = QueryMol{SmartsAtom,SmartsBond}


mapmol(::Type{T}) where T <: SDFile = mapmol(SDFileAtom,SDFileBond)
mapmol(::Type{T}) where T <: SMILES = mapmol(SmilesAtom,SmilesBond)
querymol(::Type{T}) where T <: SMARTS = querymol(SmartsAtom,SmartsBond)
