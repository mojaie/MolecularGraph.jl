
import YAML

ppath = joinpath(dirname(@__FILE__),"..", "params", "periodictable.yaml")
const ptab = YAML.load(open(ppath))


mutable struct Atom
    symbol::String
    number::Int
    std_weight::Real
    name::String
    color::Tuple
    function Atom(symbol::String)
        Atom.symbol = ptab[symbol]
        Atom.std_weight = ptab[std_weight]
        Atom.name = ptab[name]
        Atom.color = ptab[color]
        new(graph)
    end
end
