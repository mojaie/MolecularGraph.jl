#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    H_WEIGHT,
    Atom,
    sdfatom,
    smilesatom,
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
    mass::Union{Int, Nothing}
    sdf_coords::Union{SVector{3}, Nothing}
    smiles_aromatic::Union{Bool, Nothing}
    smiles_stereo::Union{Int, Nothing}
end


function sdfatom(symbol, charge, multi, mass, coords)
    if !(string(symbol) in keys(PERIODIC_TABLE))
        throw(AnnotationError("Atom '$(symbol)' not supported"))
    end
    Atom(symbol, charge, multi, mass, coords, nothing, nothing)
end


function smilesatom(symbol, charge, multi, mass, aromatic, stereo)
    if !(string(symbol) in keys(PERIODIC_TABLE))
        throw(AnnotationError("Atom '$(symbol)' not supported"))
    end
    Atom(symbol, charge, multi, mass, nothing, aromatic, stereo)
end


function atomnumber(atom::Atom)
    PERIODIC_TABLE[string(atom.symbol)]["number"]
end


function atomname(atom::Atom)
    PERIODIC_TABLE[string(atom.symbol)]["name"]
end


atomweight(atom::Atom) = PERIODIC_TABLE[string(atom.symbol)]["std_weight"]
atomweight(sym::Symbol) = PERIODIC_TABLE[string(sym)]["std_weight"]
