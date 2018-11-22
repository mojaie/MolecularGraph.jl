#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    GMapMol,
    GQueryMol,
    GVectorMol,
    SDFile,
    SMILES,
    SMARTS,
    vectormol,
    nullmol,
    getatom,
    getbond,
    neighbors,
    neighborcount,
    degree,
    atomcount,
    bondcount,
    updateatom!,
    updatebond!,
    required_annotation

import ..Graph: neighbors, neighborcount


struct GMapMol{A<:Atom,B<:Bond} <: MapMol
    graph::MutableUDGraph{A,B}
    annotation::Dict{Symbol, Annotation}
    attribute::Dict

    function GMapMol{A,B}() where {A<:Atom,B<:Bond}
        new(MutableUDGraph{A,B}(), Dict(), Dict())
    end
end

function GMapMol{A,B}(nodes::Vector{A}, edges::Vector{B}
        ) where {A<:Atom,B<:Bond}
    mol = GMapMol{A,B}()
    for (i, a) in enumerate(nodes)
        updateatom!(mol, a, i)
    end
    for (i, b) in enumerate(edges)
        updatebond!(mol, b, i)
    end
    mol
end


struct GQueryMol{A<:QueryAtom,B<:QueryBond} <: QueryMol
    graph::MutableUDGraph{A,B}
    connectivity::Array{Array{Int}}
    attribute::Dict

    function GQueryMol{A,B}() where {A<:QueryAtom,B<:QueryBond}
        new(MutableUDGraph{A,B}(), [], Dict())
    end
end


struct GVectorMol{A<:Atom,B<:Bond} <: VectorMol
    graph::UDGraph{A,B}
    v::Dict{Symbol, Array}
    annotation::Dict{Symbol, Annotation}
    attribute::Dict
end


function GVectorMol{A,B}(nodes::Vector{A}, edges::Vector{B}
        ) where {A<:Atom,B<:Bond}
    # do not use `fill`
    adj = [Dict() for i in 1:length(nodes)]
    for (i, e) in enumerate(edges)
        adj[e.u][e.v] = i
        adj[e.v][e.u] = i
    end
    GVectorMol{A,B}(
        UDGraph{A,B}(nodes, edges, adj), Dict(), Dict(), Dict()
    )
end


# Aliases
SDFile = GMapMol{SDFileAtom,SDFileBond}
SMILES = GMapMol{SmilesAtom,SmilesBond}
SMARTS = GQueryMol{SmartsAtom,SmartsBond}


nullmol(::Type{T}) where T <: SDFile = SDFile()
nullmol(::Type{T}) where T <: SMILES = SMILES()
nullmol(::Type{T}) where T <: SMARTS = SMARTS()


function vectormol(mol::GMapMol{A,B}) where {A<:Atom,B<:Bond}
    GVectorMol{A,B}(
        UDGraph{A,B}(mol.graph), Dict(), mol.annotation, mol.attribute
    )
end


getatom(mol::AbstractMol, idx) = getnode(mol.graph, idx)

getbond(mol::AbstractMol, u, v) = getedge(mol.graph, u, v)
getbond(mol::AbstractMol, idx) = getedge(mol.graph, idx)

neighbors(mol::AbstractMol, idx) = neighbors(mol.graph, idx)

atomcount(mol::AbstractMol) = nodecount(mol.graph)
bondcount(mol::AbstractMol) = edgecount(mol.graph)

neighborcount(mol::AbstractMol, idx) = length(neighbors(mol.graph, idx))
degree = neighborcount


function updateatom!(mol::AbstractMapMol, atom, idx)
    updatenode!(mol.graph, atom, idx)
end


function updatebond!(mol::AbstractMapMol, bond, idx)
    updateedge!(mol.graph, bond, idx)
end

updatebond!(m::AbstractMapMol, bond, u, v) = updatedge!(m.graph, bond, u, v)


function unlinkatom!(mol::AbstractMapMol, idx)
    unlinknode!(mol.graph, idx)
end


function unlinkbond!(mol::AbstractMapMol, idx)
    unlinkedge!(mol.graph, idx)
end

function unlinkbond!(mol::AbstractMapMol, u, v)
    unlinkedge!(mol.graph, u, v)
end


function required_annotation(mol::AbstractMol, annot)
    if !(annot in keys(mol.annotation))
        throw(ErrorException("$(annot) is not available"))
    end
end
