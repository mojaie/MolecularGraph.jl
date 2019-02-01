#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolGraphView,
    atomsubstr,
    bondsubstr


struct MolGraphView{G<:SubgraphView,M<:MolGraph} <: MolGraph
    graph::G
    molecule::M
end


function atomsubstr(mol::MolGraph, atoms)
    return MolGraphView(nodesubgraph(mol.graph, atoms), mol)
end


function bondsubstr(mol::MolGraph, bonds)
    return MolGraphView(edgesubgraph(mol.graph, bonds), mol)
end
