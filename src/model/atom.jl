#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Atom,
    number,
    name,
    color,
    weight,
    addhydrogen!


import YAML

const PERIODIC_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "..", "assets", "const", "periodictable.yaml")
))
const H_WEIGHT = PERIODIC_TABLE["H"]["std_weight"]


mutable struct Atom <: Node
    index::UInt16
    symbol::String
    charge::Int8
    multiplicity::UInt8
    mass::Union{UInt8, Nothing}
    Hcount::UInt8
    pi::UInt8
    aromatic::Bool
    Hdonor::Bool
    Hacceptor::Bool
    carbonylC::Bool
    lonepair::Bool
    wctype::UInt8
    patty::UInt8
    stereo::UInt8
    coords::Tuple
    visible::Bool

    function Atom(symbol::AbstractString)
        atom = new()
        atom.symbol = symbol
        atom.charge = 0
        atom.multiplicity = 1
        atom.mass = nothing
        atom.Hacceptor = symbol in ("N", "O", "F")
        atom.visible = symbol != "C"
        atom
    end
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


function addhydrogen!(atom::Atom, Hs::Integer)
    atom.Hcount = Hs
    atom.Hdonor = Hs > 0 && atom.symbol in ("N", "O")
end
