#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    monoisotopicmass, nominalmass, exactmass, standardweight,
    isotopiccomposition, massspecpeaks, simulatemassspec


const DEFAULT_WEIGHT_DIGITS = 2
const DEFAULT_MASS_DIGITS = 6


function molecularmass(mol::GraphMol, massfunc, implicithfunc)
    mass = 0.0
    unc = 0.0
    hm, hu = implicithfunc(:H)
    for (atom, hcnt) in zip(nodeattrs(mol), implicithcount(mol))
        m, u = massfunc(atom)
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

Return a tuple of calculated exact mass and its uncertainty.

If `number` is not given or `Atom.mass` is not specified, monoisotopic mass will be used instead.
"""
exactmass(atomsymbol::Symbol) = monoisotopicmass(atomsymbol)

function exactmass(atomsymbol::Symbol, number::Int)
    pred = rcd -> rcd["Number"] == number
    iso = ATOMTABLE[atomnumber(atomsymbol)]["Isotopes"]
    k = findfirst(pred.(iso))
    k === nothing && throw(
        ErrorException("No isotope data for $(number)$(atomsymbol)"))
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

Return calculated exact mass rounded to `digit=6`.
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

Return a tuple of standard atomic weight (or molecular weight) and its uncertainty.

If `Atom.mass` is specified, calculated exact mass of the atom will be used instead. 
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
    standardweight(::Type{Float64}, atomsymbol::Symbol) -> Float64
    standardweight(::Type{Float64}, atom::Atom) -> Float64
    standardweight(::Type{Float64}, mol::GraphMol) -> Float64

Return standard atomic weight (or molecular weight) rounded to `digit=2`.
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



"""
    isotopiccomposition(atomsymbol::Symbol, number::Int; threshold=0.001
        ) -> Vector{Tuple{Float64,Float64}}
    isotopiccomposition(mol::GraphMol; threshold=0.001
        ) -> Vector{Tuple{Float64,Float64}}

Return isotopic composition of the atoms/molecule as a vector of tuples of mass and composition.

Records that have lower abundance than the given threshold will be filtered out (default 0.001 = 0.1%)
"""
function isotopiccomposition(atomsymbol::Symbol, number::Int; threshold=0.001)
    z = atomnumber(atomsymbol)
    isotopes = []
    for iso in ATOMTABLE[z]["Isotopes"]
        cmp = iso["Composition"]
        cmp === NaN && continue
        push!(isotopes, iso)
    end
    isocnt = length(isotopes)
    data = Tuple{Float64,Float64}[]
    for seps in combinations(number + isocnt - 1, isocnt - 1)
        nums = Int[]
        prev = 0
        for sep in seps
            push!(nums, sep - prev - 1)
            prev = sep
        end
        push!(nums, number + isocnt - prev - 1)
        mass = 0.0
        dups = [logfactorial(n) for n in nums]
        comp = exp(reduce(-, dups; init=logfactorial(number)))
        for (i, iso) in enumerate(isotopes)
            mass += isotopes[i]["Mass"] * nums[i]
            comp *= isotopes[i]["Composition"] ^ nums[i]
        end
        comp > threshold && push!(data, (mass, comp))
    end
    return data
end


function isotopiccomposition(mol::GraphMol; threshold=0.001)
    data = Tuple{Float64,Float64}[]
    for tup in Iterators.product(
            (isotopiccomposition(sym, cnt; threshold=threshold)
            for (sym, cnt) in countatoms(mol))...)
        mass = 0.0
        comp = 1.0
        for t in tup
            mass += t[1]
            comp *= t[2]
        end
        comp > threshold && push!(data, (mass, comp))
    end
    return sort(data, by=x->x[1])
end



"""
    massspecpeaks(mol::GraphMol; threshold=0.001) -> Matrix{Float64}

Return a vector of tuples of each isotopic masses and their relative intensity in the simulated mass spectrum (base peak intensity = 100).

Records that have lower abundance (not peak intensity) than the given threshold will be filtered out (default 0.001 = 0.1%)
"""
function massspecpeaks(mol::GraphMol; threshold=0.001)
    isocomp = isotopiccomposition(mol; threshold=threshold)
    bp = maximum(x->x[2], isocomp)
    data = Tuple{Float64,Float64}[]
    for (mass, comp) in isocomp
        mass = round(mass, digits=DEFAULT_MASS_DIGITS)
        rint = round(comp / bp * 100, digits=ceil(Int, -log(10, threshold))-2)
        push!(data, (mass, rint))
    end
    return data
end


function gaussianpeak(mass::Float64, intensity::Float64, resolution::Int)
    fwhm = mass / resolution
    s = fwhm / (2 * sqrt(2 * log(2)))
    return x -> exp(- (x - mass)^2 / (2 * s^2)) * intensity
end


"""
    simulatemassspec(peaks::Vector{Tuple{Float64,Float64}};
        resolution=10000, rate=0.01) -> Matrix{Float64}
    simulatemassspec(mol::GraphMol;
        threshold=0.001, resolution=10000, rate=0.01) -> Matrix{Float64}

Return a matrix of simulate mass spectrum (dim 1: datapoints, dim 2: mass and intensity).

# Usage (with Plot.jl)
```julia
using MolecularGraph
using Plots
gr()
Plots.GRBackend()

mol = smilestomol("CCO")
data = simulatemassspec(mol)
plot(
    data[:, 1], data[:, 2],
    leg=false, xlabel = "Mass", ylabel = "Intensity"
)
```
"""
function simulatemassspec(
        peaks::Vector{Tuple{Float64,Float64}}; resolution=10000, rate=0.001)
    fs = Function[]
    for (mass, intensity) in peaks
        push!(fs, gaussianpeak(mass, intensity, resolution))
    end
    superposed(x) = reduce(+, [f(x) for f in fs]; init=0.0)
    xleft = round(Int, minimum(x->x[1], peaks) - 1)
    xright = round(Int, maximum(x->x[1], peaks) + 1)
    arr = collect(xleft:rate:xright)
    return hcat(arr, superposed.(arr))
end

simulatemassspec(mol::GraphMol; kwargs...
    ) =  simulatemassspec(massspecpeaks(mol; kwargs...); kwargs...)
