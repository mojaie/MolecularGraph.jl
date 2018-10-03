
using GraphMol.GraphModel


mutable struct MolecularGraph
    graph::UndirectedGraph
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


function atom(mol::MolecularGraph, idx::Int32)
    mol.graph.nodes[idx]["atom"]
end


function updateAtom!(mol::MolecularGraph, idx::Int32, atom::Atom)
    mol.graph.nodes[idx]["atom"] = atom
end
