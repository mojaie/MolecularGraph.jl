#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFileAtom, SmilesAtom, SmartsAtom,
    setcharge, setstereo, atomnumber


const ATOMTABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "../../assets/const/atomicweights.yaml")
))

const ATOMSYMBOLMAP = YAML.load(open(
    joinpath(dirname(@__FILE__), "../../assets/const/symboltonumber.yaml")
))



struct SDFileAtom <: Atom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    mass::Union{Int, Nothing}
    coords::Union{Vector{Float64}, Nothing}
    stereo::Symbol

    function SDFileAtom(sym, chg, multi, mass, coords, stereo)
        new(sym, chg, multi, mass, coords, stereo)
    end
end

SDFileAtom() = SDFileAtom(:C, 0, 1, nothing, nothing)
SDFileAtom(sym) = SDFileAtom(sym, 0, 1, nothing, nothing)
SDFileAtom(sym, chg) = SDFileAtom(sym, chg, 1, nothing, nothing)
SDFileAtom(sym, chg, multi, mass, coords
    ) = SDFileAtom(sym, chg, multi, mass, coords, :unspecified)

setcharge(a::SDFileAtom, chg
    ) = SDFileAtom(a.symbol, chg, a.multiplicity, a.mass, a.coords)

setstereo(a::SDFileAtom, direction) = SDFileAtom(
    a.symbol, a.charge, a.multiplicity, a.mass, a.coords, direction)


struct SmilesAtom <: Atom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    mass::Union{Float64, Nothing}
    isaromatic::Union{Bool, Nothing}
    stereo::Symbol

    function SmilesAtom(sym, chg, multi, mass, aromatic, stereo)
        new(sym, chg, multi, mass, aromatic, stereo)
    end
end

SmilesAtom() = SmilesAtom(:C, 0, 1, nothing, nothing, :unspecified)
SmilesAtom(sym) = SmilesAtom(sym, 0, 1, nothing, nothing, :unspecified)
SmilesAtom(sym, chg) = SmilesAtom(sym, chg, 1, nothing, nothing, :unspecified)

setcharge(a::SmilesAtom, chg) = SmilesAtom(
    a.symbol, chg, a.multiplicity, a.mass, a.isaromatic, a.stereo)

setstereo(a::SmilesAtom, direction) = SmilesAtom(
    a.symbol, a.charge, a.multiplicity, a.mass, a.isaromatic, direction)


struct SmartsAtom <: QueryAtom
    query::Pair
end



"""
    atomnumber(atomsymbol::Symbol) -> Int
    atomnumber(atom::Atom) -> Int

Return atom number.
"""
atomnumber(atomsymbol::Symbol) = ATOMSYMBOLMAP[string(atomsymbol)]
atomnumber(atom::Atom) = atomnumber(atom.symbol)



"""
    atomsymbol(n::Int) -> Symbol

Return atom symbol of given atomic number.
"""
atomsymbol(n::Int) = Symbol(ATOMTABLE[n]["Symbol"])
