#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const DEFAULT_WEIGHT_DIGITS = 2
const DEFAULT_MASS_DIGITS = 6


function atom_mass_unc(atom::AbstractAtom, massfunc::F) where F
    return massfunc(atom_symbol(atom), atom_mass(atom))
end


function molecular_mass_unc(mol::SimpleMolGraph, massfunc::F) where F
    mass = 0.0
    unc = 0.0
    hm, hu = massfunc(:H)
    imp_hs = implicit_hydrogens(mol)
    for i in vertices(mol)
        m, u = atom_mass_unc(props(mol, i), massfunc)
        mass += m + hm * imp_hs[i]
        unc += u + hu * imp_hs[i]
    end
    return (mass, unc)
end


function molecular_mass_unc(counter::Dict{Symbol,Int}, massfunc::F) where F
    mass = 0.0
    unc = 0.0
    for (sym, cnt) in counter
        m, u = massfunc(sym)
        mass += m * cnt
        unc += u * cnt
    end
    return (mass, unc)
end



"""
    monoiso_mass_unc(atomsymbol::Symbol) -> Tuple{Float64,Float64}
    monoiso_mass_unc(atom) -> Tuple{Float64,Float64}
    monoiso_mass_unc(mol::MolGraph) -> Tuple{Float64,Float64}

Return a tuple of monoisotopic mass of the atom/molecule and its uncertainty.

Monoisotopic mass is the relative atomic mass of the most abundant isotope. Even if there is specific `Atom.mass` value, it will be ignored.
"""
function monoiso_mass_unc(atomsymbol::Symbol, number::Union{Int, Nothing}=nothing)
    num = atom_number(atomsymbol)
    mass = ATOMTABLE[num]["Monoisotopic"]
    unc = ATOMTABLE[num]["MonoisotopicUncertainty"]
    return (mass, unc)
end

monoiso_mass_unc(atom::AbstractAtom) = atom_mass_unc(atom, monoiso_mass_unc)
monoiso_mass_unc(mol::SimpleMolGraph) = molecular_mass_unc(mol, monoiso_mass_unc)


"""
    monoiso_mass(atomsymbol::Symbol, [digits::Int]) -> Float64
    monoiso_mass(atom, [digits::Int]) -> Float64
    monoiso_mass(mol::MolGraph, [digits::Int]) -> Float64

Return monoisotopic mass of the atom/molecule.
"""
monoiso_mass(atomsymbol::Symbol) = monoiso_mass_unc(atomsymbol)[1]
monoiso_mass(atom::AbstractAtom) = monoiso_mass_unc(atom)[1]
monoiso_mass(mol::SimpleMolGraph) = monoiso_mass_unc(mol)[1]

monoiso_mass(atomsymbol::Symbol, digits::Int) = round(monoiso_mass(atomsymbol), digits=digits)
monoiso_mass(atom::AbstractAtom, digits::Int) = round(monoiso_mass(atom), digits=digits)
monoiso_mass(mol::SimpleMolGraph, digits::Int) = round(monoiso_mass(mol), digits=digits)


"""
    nominal_mass(atomsymbol::Symbol) -> Int
    nominal_mass(atom) -> Int
    nominal_mass(mol::MolGraph) -> Int

Return nominal mass of the atom/molecule.
"""
nominal_mass(atomsymbol::Symbol) = round(Int, monoiso_mass(atomsymbol))
nominal_mass(atom::AbstractAtom) = round(Int, monoiso_mass(atom))
nominal_mass(mol::SimpleMolGraph) = round(Int, monoiso_mass(mol))


"""
    exact_mass_unc(atomsymbol::Symbol, [number::Union{Int, Nothing}]) -> Tuple{Float64,Float64}
    exact_mass_unc(atom) -> Tuple{Float64,Float64}
    exact_mass_unc(mol::MolGraph) -> Tuple{Float64,Float64}

Return a tuple of calculated exact mass and its uncertainty.

If `number` is not given or `Atom.mass` is not specified, monoisotopic mass will be used instead.
"""
function exact_mass_unc(atomsymbol::Symbol, number::Union{Int, Nothing}=nothing)
    number === nothing && return monoiso_mass_unc(atomsymbol)
    pred = rcd -> rcd["Number"] == number
    iso = ATOMTABLE[atom_number(atomsymbol)]["Isotopes"]
    k = findfirst(pred.(iso))
    k === nothing && error("No isotope data for $(number)$(atomsymbol)")
    mass = iso[k]["Mass"]
    unc = iso[k]["MassUncertainty"]
    return (mass, unc)
end

exact_mass_unc(atom::AbstractAtom) = atom_mass_unc(atom, exact_mass_unc)
exact_mass_unc(mol::SimpleMolGraph) = molecular_mass_unc(mol, exact_mass_unc)


"""
    exact_mass(atomsymbol::Symbol, [digits::Int]) -> Float64
    exact_mass(atom, [digits::Int]) -> Float64
    exact_mass(mol::MolGraph, [digits::Int]) -> Float64

Return calculated exact mass.
"""
exact_mass(atomsymbol::Symbol) = exact_mass_unc(atomsymbol)[1]
exact_mass(atom::AbstractAtom) = exact_mass_unc(atom)[1]
exact_mass(mol::SimpleMolGraph) = exact_mass_unc(mol)[1]

exact_mass(atomsymbol::Symbol, digits::Int) = round(exact_mass(atomsymbol), digits=digits)
exact_mass(atom::AbstractAtom, digits::Int) = round(exact_mass(atom), digits=digits)
exact_mass(mol::SimpleMolGraph, digits::Int) = round(exact_mass(mol), digits=digits)


"""
    standard_weight_unc(atomsymbol::Symbol) -> Tuple{Float64,Float64}
    standard_weight_unc(atom) -> Tuple{Float64,Float64}
    standard_weight_unc(mol::MolGraph) -> Tuple{Float64,Float64}

Return a tuple of standard atomic weight (or molecular weight) and its uncertainty.

If `Atom.mass` is specified, calculated exact mass of the atom will be used instead. 
"""
function standard_weight_unc(atomsymbol::Symbol, number::Union{Int, Nothing}=nothing)
    number === nothing || return exact_mass_unc(atomsymbol, number)
    num = atom_number(atomsymbol)
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

standard_weight_unc(atom::AbstractAtom) = atom_mass_unc(atom, standard_weight_unc)
standard_weight_unc(mol::SimpleMolGraph) = molecular_mass_unc(mol, standard_weight_unc)

"""
    standard_weight(atomsymbol::Symbol, [digits::Int]) -> Float64
    standard_weight(atom, [digits::Int]) -> Float64
    standard_weight(mol::MolGraph, [digits::Int]) -> Float64

Return standard atomic weight (or molecular weight).
"""
standard_weight(atomsymbol::Symbol) = standard_weight_unc(atomsymbol)[1]
standard_weight(atom::AbstractAtom) = standard_weight_unc(atom)[1]
standard_weight(mol::SimpleMolGraph) = standard_weight_unc(mol)[1]
standard_weight(atomsymbol::Symbol, digits::Int) = round(standard_weight(atomsymbol), digits=digits)
standard_weight(atom::AbstractAtom, digits::Int) = round(standard_weight(atom), digits=digits)
standard_weight(mol::SimpleMolGraph, digits::Int) = round(standard_weight(mol), digits=digits)


"""
    isotopiccomposition(atomsymbol::Symbol, number::Int; threshold=0.001
        ) -> Vector{Tuple{Float64,Float64}}
    isotopiccomposition(mol::MolGraph; threshold=0.001
        ) -> Vector{Tuple{Float64,Float64}}

Return isotopic composition of the atoms/molecule as a vector of tuples of mass and composition.

Records that have lower abundance than the given threshold will be filtered out (default 0.001 = 0.1%)
"""
function isotopic_composition(atomsymbol::Symbol, number::Int; threshold=0.001)
    z = atom_number(atomsymbol)
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


function isotopic_composition(mol::SimpleMolGraph; threshold=0.001)
    data = Tuple{Float64,Float64}[]
    for tup in Iterators.product(
            (isotopic_composition(sym, cnt; threshold=threshold)
            for (sym, cnt) in atom_counter(mol))...)
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
    massspec_peaks(mol::MolGraph; threshold=0.001) -> Matrix{Float64}

Return a vector of tuples of each isotopic masses and their relative intensity in the simulated mass spectrum (base peak intensity = 100).

Records that have lower abundance (not peak intensity) than the given threshold will be filtered out (default 0.001 = 0.1%)
"""
function massspec_peaks(mol::SimpleMolGraph; threshold=0.001)
    isocomp = isotopic_composition(mol; threshold=threshold)
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
    simulate_massspec(peaks::Vector{Tuple{Float64,Float64}};
        resolution=10000, rate=0.01) -> Matrix{Float64}
    simulate_massspec(mol::MolGraph;
        threshold=0.001, resolution=10000, rate=0.01) -> Matrix{Float64}

Return a matrix of simulate mass spectrum (dim 1: datapoints, dim 2: mass and intensity).

Note that the peaks are just calculated from the isotopic composition of atoms (not intended for simulation of fragmentation).

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
function simulate_massspec(
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

simulate_massspec(mol::SimpleMolGraph; kwargs...
    ) = simulate_massspec(massspec_peaks(mol; kwargs...); kwargs...)



# deprecated names
export
    monoisotopicmass, nominalmass, exactmass,
    isotopiccomposition, massspecpeaks, simulatemassspec

monoisotopicmass = monoiso_mass
nominalmass = nominal_mass
exactmass = exact_mass
isotopiccomposition = isotopic_composition
massspecpeaks = massspec_peaks
simulatemassspec = simulate_massspec