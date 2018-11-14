#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    PERIODIC_TABLE,
    H_WEIGHT,
    Atom,
    sdfatom,
    smilesatom,
    atomsymbol,
    atomnumber,
    atomname,
    atomweight


const PERIODIC_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "..", "assets", "const", "periodictable.yaml")
))
const H_WEIGHT = PERIODIC_TABLE["H"]["std_weight"]


struct Atom <: AbstractNode
    symbol::Symbol
    charge::Int
    multiplicity::Int
    mass::Union{Float64, Nothing}
    sdf_coords::Union{SVector{3}, Nothing}
    smiles_aromatic::Union{Bool, Nothing}
    smiles_stereo::Union{Int, Nothing}
end


function sdfatom(symbol, charge, multi, mass, coords)
    if !(string(symbol) in keys(PERIODIC_TABLE))
        throw(IOError("unsupported symbol: $(symbol)"))
    end
    Atom(symbol, charge, multi, mass, coords, nothing, nothing)
end


function smilesatom(symbol, charge, multi, mass, aromatic, stereo)
    if !(string(symbol) in keys(PERIODIC_TABLE))
        throw(IOError("unsupported symbol: $(symbol)"))
    end
    Atom(symbol, charge, multi, mass, nothing, aromatic, stereo)
end


function atomsymbol(number::Int)
    for (sym, atom) in PERIODIC_TABLE
        if atom["number"] == number
            return Symbol(sym)
        end
    end
    throw(OperationError("invalid atomic number $(number)"))
end


atomnumber(symbol::Symbol) = PERIODIC_TABLE[string(symbol)]["number"]
atomnumber(atom::Atom) = PERIODIC_TABLE[string(atom.symbol)]["number"]


atomname(symbol::Symbol) = PERIODIC_TABLE[string(symbol)]["name"]
atomname(atom::Atom) = PERIODIC_TABLE[string(atom.symbol)]["name"]


function atomweight(atom::Atom)
    stdweight = PERIODIC_TABLE[string(atom.symbol)]["std_weight"]
    atom.mass === nothing ? stdweight : atom.mass
end
