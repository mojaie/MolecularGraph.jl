#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    GeneralMapMol, GeneralVectorMol,
    ConnectedQueryMol, DisconnectedQueryMol,
    SDFile, SMILES, ConnectedSMARTS, SMARTS,
    vectormol, nullmol


struct GeneralMapMol{A<:Atom,B<:Bond} <: MapMolGraph
    graph::MapUDGraph{A,B}
    annotation::Dict{Symbol, Annotation}
    attribute::Dict

    function GeneralMapMol{A,B}() where {A<:Atom,B<:Bond}
        new(MapUDGraph{A,B}(), Dict(), Dict())
    end
end

function GeneralMapMol{A,B}(nodes::Vector{A}, edges::Vector{B}
        ) where {A<:Atom,B<:Bond}
    mol = GeneralMapMol{A,B}()
    for (i, a) in enumerate(nodes)
        updateatom!(mol, a, i)
    end
    for (i, b) in enumerate(edges)
        updatebond!(mol, b, i)
    end
    return mol
end


struct ConnectedQueryMol{A<:QueryAtom,B<:QueryBond} <: QueryMolGraph
    graph::MapUDGraph{A,B}
    attribute::Dict

    function ConnectedQueryMol{A,B}() where {A<:QueryAtom,B<:QueryBond}
        new(MapUDGraph{A,B}(), Dict())
    end
end


struct DisconnectedQueryMol{A<:QueryAtom,B<:QueryBond} <: QueryMolGraph
    graph::MapUDGraph{A,B}
    connectivity::Array{Array{Int}}
    attribute::Dict

    function DisconnectedQueryMol{A,B}() where {A<:QueryAtom,B<:QueryBond}
        new(MapUDGraph{A,B}(), [], Dict())
    end
end


struct GeneralVectorMol{A<:Atom,B<:Bond} <: VectorMolGraph
    graph::VectorUDGraph{A,B}
    vector::Dict{Symbol, Array}
    annotation::Dict{Symbol, Annotation}
    attribute::Dict
end


function GeneralVectorMol{A,B}(nodes::Vector{A}, edges::Vector{B}
        ) where {A<:Atom,B<:Bond}
    # do not use `fill`
    adj = [Dict() for i in 1:length(nodes)]
    for (i, e) in enumerate(edges)
        adj[e.u][e.v] = i
        adj[e.v][e.u] = i
    end
    return GeneralVectorMol{A,B}(
        VectorUDGraph{A,B}(nodes, edges, adj), Dict(), Dict(), Dict()
    )
end


function vectormol(mol::GeneralMapMol{A,B}) where {A<:Atom,B<:Bond}
    return GeneralVectorMol{A,B}(
        VectorUDGraph{A,B}(mol.graph), Dict(), mol.annotation, mol.attribute
    )
end


# Aliases

# TODO: use traits
SDFile = GeneralMapMol{SDFileAtom,SDFileBond}
SMILES = GeneralMapMol{SmilesAtom,SmilesBond}
ConnectedSMARTS = ConnectedQueryMol{SmartsAtom,SmartsBond}
SMARTS = DisconnectedQueryMol{SmartsAtom,SmartsBond}


nullmol(::Type{T}) where T <: SDFile = SDFile()
nullmol(::Type{T}) where T <: SMILES = SMILES()
nullmol(::Type{T}) where T <: SMARTS = SMARTS()
