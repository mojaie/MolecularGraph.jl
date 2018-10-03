

import YAML

const PERIODIC_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "..", "const", "periodictable.yaml")
))
const H_WEIGHT = PERIODIC_TABLE["H"]["std_weight"]


mutable struct Atom
    symbol::String
    charge::Int
    multiplicity::Int
    mass::Union{Int, Nothing}
    Hcount::Int
    pi::Int
    aromatic::Bool
    Hdonor::Bool
    Hacceptor::Bool
    carbonylC::Bool
    lonepair::Bool
    wctype::Int
    patty::Int
    stereo::Int
    coords::Tuple
    visible::Bool


    function Atom(symbol::String)
        initialize!(new(), symbol)
    end
end


function initialize!(atom::Atom, symbol)
    atom.symbol = symbol
    atom.charge = 0
    atom.multiplicity = 1
    atom.mass = nothing
    atom.Hacceptor = symbol in ("N", "O", "F")
    atom.visible = symbol != "C"
    atom
end


function number(atom::Atom)
    PERIODIC_TABLE[atom.symbol]["number"]
end


function name(atom::Atom)
    PERIODIC_TABLE[atom.symbol]["name"]
end


function color(atom::Atom)
    attr = PERIODIC_TABLE[atom.symbol]
    tuple(get(attr, "color", [0, 192, 192]))
end


function weight(atom::Atom)
    m = PERIODIC_TABLE[atom.symbol]["std_weight"]
    m + H_WEIGHT * atom.Hcount
end


function addH!(atom::Atom, Hs::Int)
    atom.Hcount = Hs
    atom.Hdonor = Hs > 0 && atom.symbol in ("N", "O")
end
