
import YAML

const ptab = YAML.load(open("./params/periodictable.yaml"))


mutable struct Atom
    symbol::Char
    number::Int
    std_weight::Real
    name::String
    color::Tuple
    function Atom(symbol::Char)
        Atom.symbol = ptab[symbol]
        Atom.std_weight = ptab[std_weight]
        Atom.name = ptab[name]
        Atom.color = ptab[color]
        new(graph)
    end
end


function atom(mol::MoleculeGraph, idx::Int32)
    mol.graph.nodes[idx]["atom"]
end


function updateAtom!(mol::MoleculeGraph, idx::Int32, atom::Atom)
    mol.graph.nodes[idx]["atom"] = atom
end
