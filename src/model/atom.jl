#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    PERIODIC_TABLE,
    H_WEIGHT,
    SDFileAtom,
    SmilesAtom,
    SmartsAtom,
    atomsymbol,
    atomnumber,
    atomname,
    atomweight


const PERIODIC_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "..", "assets", "const", "periodictable.yaml")
))
const H_WEIGHT = PERIODIC_TABLE["H"]["std_weight"]


mutable struct SDFileAtom <: Atom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    mass::Union{Float64, Nothing}
    coords::Union{SVector{3}, Nothing}

    function SDFileAtom(sym, chg, multi, mass, coords)
        if !(string(sym) in keys(PERIODIC_TABLE))
            throw(ErrorException("unsupported symbol: $(sym)"))
        end
        new(sym, chg, multi, mass, coords)
    end
end

SDFileAtom(sym) = SDFileAtom(sym, 0, 1, nothing, nothing)
SDFileAtom(sym, chg) = SDFileAtom(sym, chg, 1, nothing, nothing)


mutable struct SmilesAtom <: Atom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    mass::Union{Float64, Nothing}
    isaromatic::Union{Bool, Nothing}
    stereo::Union{Int, Nothing}

    function SmilesAtom(sym, chg, multi, mass, aromatic, stereo)
        if !(string(sym) in keys(PERIODIC_TABLE))
            throw(ErrorException("unsupported symbol: $(sym)"))
        end
        new(sym, chg, multi, mass, aromatic, stereo)
    end
end

SmilesAtom(sym) = SmilesAtom(sym, 0, 1, nothing, nothing, nothing)
SmilesAtom(sym, chg) = SmilesAtom(sym, chg, 1, nothing, nothing, nothing)


mutable struct SmartsAtom <: QueryAtom
    query::Pair
end


function atomsymbol(number::Int)
    for (sym, atom) in PERIODIC_TABLE
        if atom["number"] == number
            return Symbol(sym)
        end
    end
    throw(ErrorException("invalid atomic number $(number)"))
end


atomnumber(symbol::Symbol) = PERIODIC_TABLE[string(symbol)]["number"]
atomnumber(atom::Atom) = PERIODIC_TABLE[string(atom.symbol)]["number"]


atomname(symbol::Symbol) = PERIODIC_TABLE[string(symbol)]["name"]
atomname(atom::Atom) = PERIODIC_TABLE[string(atom.symbol)]["name"]


function atomweight(atom::Atom)
    stdweight = PERIODIC_TABLE[string(atom.symbol)]["std_weight"]
    atom.mass === nothing ? stdweight : atom.mass
end
