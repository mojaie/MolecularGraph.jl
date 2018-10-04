#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using GraphMol.GraphModel


mutable struct MolecularGraph
    graph::UndirectedGraph
    descriptors::set
    function MolecularGraph(graph::UndirectedGraph)
        initialize!(new(), graph)
    end
end


MolecularGraph() = MolecularGraph(UndirectedGraph())


function initialize!(mol::MolecularGraph,
                     graph::UndirectedGraph)
    mol.graph = graph
    mol
end


function getatom(mol::MolecularGraph, idx::Int32)
    mol.graph.nodes[idx]["atom"]
end


function updateAtom!(mol::MolecularGraph, idx::Int32, atom::Atom)
    mol.graph.nodes[idx]["atom"] = atom
end


function required_descriptor(mol::MolecularGraph, desc::AbstractString)
    if desc ∉ mol.descriptors
        throw(ErrorException("$() is not assigned"))
    end
end
