module MoleculeGraph

using UndirectedGraph
using Atom

struct MoleculeGraph
    graph::Graph
    function MoleculeGraph(graph::Graph)
        new(graph)
    end
end


function atom(mol::MoleculeGraph, idx::Int32)
    mol.graph.nodes[idx]["atom"]
end


function updateAtom!(mol::MoleculeGraph, idx::Int32, atom::Atom)
    mol.graph.nodes[idx]["atom"] = atom
end
