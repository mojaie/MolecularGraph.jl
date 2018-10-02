

using UndirectedGraph
using Atom

struct MolecularGraph
    graph::Graph
    function MolecularGraph(graph::Graph)
        new(graph)
    end
end


function atom(mol::MolecularGraph, idx::Int32)
    mol.graph.nodes[idx]["atom"]
end


function updateAtom!(mol::MolecularGraph, idx::Int32, atom::Atom)
    mol.graph.nodes[idx]["atom"] = atom
end
