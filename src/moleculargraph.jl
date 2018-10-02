
using GraphMol.GraphModel
using GraphMol.MolecularModel

struct MolecularGraph
    graph::UndirectedGraph
    function MolecularGraph(graph::UndirectedGraph)
        new(graph)
    end
end


function atom(mol::MolecularGraph, idx::Int32)
    mol.graph.nodes[idx]["atom"]
end


function updateAtom!(mol::MolecularGraph, idx::Int32, atom::Atom)
    mol.graph.nodes[idx]["atom"] = atom
end
