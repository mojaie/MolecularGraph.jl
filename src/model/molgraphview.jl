#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolGraphView,
    atomsubstr,
    bondsubstr,
    required_annotation


struct MolGraphView{G<:SubgraphView,M<:AbstractMol} <: AbstractMol
    graph::G
    molecule::M
end


function atomsubstr(mol::AbstractMol, atoms)
    return MolGraphView(nodesubgraph(mol.graph, atoms), mol)
end


function bondsubstr(mol::AbstractMol, bonds)
    return MolGraphView(edgesubgraph(mol.graph, bonds), mol)
end


required_annotation(
    view::MolGraphView, annot) = required_annotation(view.molecule, annot)
