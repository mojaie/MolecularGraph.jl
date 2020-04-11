#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    monoisotopicmass, nominalmass, exactmass, standardweight


const DEFAULT_WEIGHT_DIGITS = 2
const DEFAULT_MASS_DIGITS = 6


function molecularmass(mol::GraphMol, massfunc, implicithfunc)
    mass = 0.0
    unc = 0.0
    hm, hu = implicithfunc(:H)
    for (sym, hcnt) in zip(atomsymbol(mol), implicithcount(mol))
        m, u = massfunc(sym)
        mass += m + hm * hcnt
        unc += u + hu * hcnt
    end
    return (mass, unc)
end



"""
    monoisotopicmass(atomsymbol::Symbol) -> Tuple{Float64,Float64}
    monoisotopicmass(atom::Atom) -> Tuple{Float64,Float64}
    monoisotopicmass(mol::GraphMol) -> Tuple{Float64,Float64}

Return a tuple of monoisotopic mass of the atom/molecule and its uncertainty.

Monoisotopic mass is the relative atomic mass of the most abundant isotope. Even if there is specific `Atom.mass` value, it will be ignored.
"""
function monoisotopicmass(atomsymbol::Symbol)
    num = atomnumber(atomsymbol)
    mass = ATOMTABLE[num]["Monoisotopic"]
    unc = ATOMTABLE[num]["MonoisotopicUncertainty"]
    return (mass, unc)
end

monoisotopicmass(atom::Atom) = monoisotopicmass(atom.symbol)
monoisotopicmass(mol::GraphMol
    ) = molecularmass(mol, monoisotopicmass, monoisotopicmass)


"""
    monoisotopicmass(::Type{Float64}, atomsymbol::Symbol) -> Float64
    monoisotopicmass(::Type{Float64}, atom::Atom) -> Float64
    monoisotopicmass(::Type{Float64}, mol::GraphMol) -> Float64

Return monoisotopic mass of the atom/molecule rounded to `digits=6`.
"""
function monoisotopicmass(::Type{Float64}, atomsymbol::Symbol)
    mass, unc = monoisotopicmass(atomsymbol)
    mass === NaN && return NaN
    return round(mass, digits=DEFAULT_MASS_DIGITS)
end

monoisotopicmass(::Type{Float64}, atom::Atom
    ) = monoisotopicmass(Float64, atom.symbol)

function monoisotopicmass(::Type{Float64}, mol::GraphMol)
    mass, unc = molecularmass(mol, monoisotopicmass, monoisotopicmass)
    mass === NaN && return NaN
    return round(mass, digits=DEFAULT_MASS_DIGITS)
end


"""
    monoisotopicmass(::Type{T}, atomsymbol::Symbol) where {T<:Integer} -> T
    monoisotopicmass(::Type{T}, atom::Atom) where {T<:Integer} -> T
    monoisotopicmass(::Type{T}, mol::GraphMol) where {T<:Integer} -> T
    nominalmass(atomsymbol::Symbol) -> Int
    nominalmass(atom::Atom) -> Int
    nominalmass(mol::GraphMol) -> Int

Return nominal mass of the atom/molecule.
"""
monoisotopicmass(::Type{T}, atomsymbol::Symbol
    ) where {T<:Integer} = round(T, monoisotopicmass(Float64, atomsymbol))

monoisotopicmass(::Type{T}, atom::Atom
    ) where {T<:Integer} = round(T, monoisotopicmass(Float64, atom.symbol))

monoisotopicmass(::Type{T}, mol::GraphMol
    ) where {T<:Integer} = round(T, monoisotopicmass(Float64, mol))

nominalmass(atomsymbol::Symbol) = monoisotopicmass(Int, atomsymbol)
nominalmass(atom::Atom) = monoisotopicmass(Int, atom)
nominalmass(mol::GraphMol) = monoisotopicmass(Int, mol)



"""
    exactmass(atomsymbol::Symbol) -> Tuple{Float64,Float64}
    exactmass(atomsymbol::Symbol, number::Int) -> Tuple{Float64,Float64}
    exactmass(atom::Atom) -> Tuple{Float64,Float64}
    exactmass(mol::GraphMol) -> Tuple{Float64,Float64}

Return a tuple of calculated exact mass (relative atomic mass) and its uncertainty.

If `number` is not given or `Atom.mass` is not specified, monoisotopic mass will be used instead.
"""
exactmass(atomsymbol::Symbol) = monoisotopicmass(atomsymbol)

function exactmass(atomsymbol::Symbol, number::Int)
    pred = rcd -> rcd["Number"] == number
    iso = ATOMTABLE[atomnumber(atomsymbol)]["Isotopes"]
    k = findfirst(pred.(iso))
    mass = iso[k]["Mass"]
    unc = iso[k]["MassUncertainty"]
    return (mass, unc)
end

function exactmass(atom::Atom)
    atom.mass === nothing || return exactmass(atom.symbol, atom.mass)
    return monoisotopicmass(atom.symbol)
end

exactmass(mol::GraphMol) = molecularmass(mol, exactmass, monoisotopicmass)


"""
    exactmass(::Type{Float64}, atomsymbol::Symbol) -> Float64
    exactmass(::Type{Float64}, atomsymbol::Symbol, number::Int) -> Float64
    exactmass(::Type{Float64}, atom::Atom) -> Float64
    exactmass(::Type{Float64}, mol::GraphMol) -> Float64

Return calculated exact mass (relative atomic mass) rounded to `digit=6`.
"""
exactmass(::Type{Float64}, atomsymbol::Symbol
    ) = monoisotopicmass(Float64, atomsymbol)

function exactmass(::Type{Float64}, atomsymbol::Symbol, number::Int)
    mass, unc = exactmass(atomsymbol, number)
    mass === NaN && return NaN
    return round(mass, digits=DEFAULT_MASS_DIGITS)
end

function exactmass(::Type{Float64}, atom::Atom)
    atom.mass === nothing || return exactmass(Float64, atom.symbol, atom.mass)
    return monoisotopicmass(Float64, atom.symbol)
end

function exactmass(::Type{Float64}, mol::GraphMol)
    mass, unc = molecularmass(mol, exactmass, monoisotopicmass)
    mass === NaN && return NaN
    return round(mass, digits=DEFAULT_MASS_DIGITS)
end



"""
    standardweight(atomsymbol::Symbol) -> Tuple{Float64,Float64}
    standardweight(atom::Atom) -> Tuple{Float64,Float64}
    standardweight(mol::GraphMol) -> Tuple{Float64,Float64}

Return a tuple of standard atomic/molecular weight (average weight based on natural abundance) and its uncertainty.

If `Atom.mass` is specified, relative atomic mass will be used instead of standard weight. 
"""
function standardweight(atomsymbol::Symbol)
    num = atomnumber(atomsymbol)
    wt = ATOMTABLE[num]["Weight"]
    unctype = ATOMTABLE[num]["WeightType"]
    if unctype == "Interval"
        lower = ATOMTABLE[num]["WeightLower"]
        higher = ATOMTABLE[num]["WeightHigher"]
        unc = max(higher - wt, wt - lower)
    elseif unctype == "Uncertainty"
        unc = ATOMTABLE[num]["WeightUncertainty"]
    else
        unc = NaN  # unctype in ("MostStable", "NotAvailable")
    end
    return (wt, unc)
end

function standardweight(atom::Atom)
    atom.mass === nothing || return exactmass(atom.symbol, atom.mass)
    return standardweight(atom.symbol)
end

standardweight(mol::GraphMol
    ) = molecularmass(mol, standardweight, standardweight)



"""
    standardweight(atomsymbol::Symbol) -> Float64
    standardweight(atom::Atom) -> Float64
    standardweight(mol::GraphMol) -> Float64

Return standard atomic/molecular weight (average weight based on natural abundance) rounded to `digit=2`.
"""
function standardweight(::Type{Float64}, atomsymbol::Symbol)
    wt, unc = standardweight(atomsymbol)
    wt === NaN && return NaN
    return round(wt, digits=DEFAULT_WEIGHT_DIGITS)
end

function standardweight(::Type{Float64}, atom::Atom)
    atom.mass === nothing || return exactmass(Float64, atom.symbol, atom.mass)
    return standardweight(Float64, atom.symbol)
end

function standardweight(::Type{Float64}, mol::GraphMol)
    wt, unc = molecularmass(mol, standardweight, standardweight)
    wt === NaN && return NaN
    return round(wt, digits=DEFAULT_WEIGHT_DIGITS)
end
