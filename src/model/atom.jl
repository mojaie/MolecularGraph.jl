#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# Constants

@kwdef mutable struct Isotope
    Composition::Float64 = NaN
    Mass::Float64 = NaN
    MassUncertainty::Float64 = NaN
    CompositionUncertainty::Float64 = NaN
    Number::Int = 0
    Symbol::String = ""
end

Base.getindex(iso::Isotope, k::String) = getproperty(iso, Symbol(k))


@kwdef mutable struct AtomTable
    WeightType::String = ""
    Isotopes::Vector{Isotope} = Isotope[]
    WeightHigher::Float64 = NaN
    WeightLower::Float64 = NaN
    WeightUncertainty::Float64 = NaN
    Notes::Vector{String} = String[]
    Monoisotopic::Float64 = NaN
    Number::Int = 0
    Weight::Float64 = NaN
    Symbol::String = ""
    MonoisotopicUncertainty::Float64 = NaN
end

Base.getindex(a::AtomTable, k::String) = getproperty(a, Symbol(k))

"""
    ATOMTABLE

Vector of atomic data (atomic weights, monoisotopic masses and isotopes) indexed by
atomic number (e.g. `ATOMTABLE[6].Weight`).
"""
const ATOMTABLE = let
    weightsfile = joinpath(dirname(@__FILE__), "../../assets/const/atomicweights.yaml")
    include_dependency(weightsfile)
    tbl = AtomTable[]
    for rcd in YAML.load(open(weightsfile))
        a = AtomTable()
        for (k, v) in rcd
            if k == "Isotopes"
                ic = Isotope[]
                for iso in v
                    i = Isotope()
                    for (j, u) in iso
                        setproperty!(i, Symbol(j), u)
                    end
                    push!(ic, i)
                end
                setproperty!(a, Symbol(k), ic)
            elseif k == "Notes" && isnothing(v)
                continue
            else
                setproperty!(a, Symbol(k), v)
            end
        end
        push!(tbl, a)
    end
    tbl
end

"""
    ATOMSYMBOLMAP

Mapping of atom symbols (`Symbol`) to atomic numbers (e.g. `ATOMSYMBOLMAP[:C] == 6`).
"""
const ATOMSYMBOLMAP = let
    symbolfile = joinpath(dirname(@__FILE__), "../../assets/const/symboltonumber.yaml")
    include_dependency(symbolfile)
    YAML.load(open(symbolfile); dicttype=Dict{Symbol,Int})
end

"""
    ATOM_COVALENT_RADII

Mapping of atomic numbers to covalent radii (in Å). Elements with multiple values have
a `Dict{String,Float32}` (`"Csp3"`, `"Csp2"` and `"Csp"` for carbon, low spin `"l.s."` and
high spin `"h.s."` for some metals). See also [`atom_radius`](@ref).

# Reference

1. Cordero, B. et al. Covalent radii revisited. Dalton Trans. 21, 2832–2838 (2008).
"""
const ATOM_COVALENT_RADII = let
    radiifile = joinpath(dirname(@__FILE__), "../../assets/const/covalent_radii.csv")
    include_dependency(radiifile)
    tab, headers = readdlm(radiifile, '\t', String; header=true, comments=true)
    radii = Dict{Int,Union{Float32,Dict{String,Float32}}}()
    for i = 1:size(tab, 1)
        an = parse(Int, tab[i, 1])
        container = if !haskey(radii, an)
            if i < size(tab, 1) && tab[i+1,1] == tab[i,1]  # special handling for elements with multiple options
                radii[an] = Dict{String,Float32}()
            else
                nothing
            end
        else
            radii[an]
        end
        ar = parse(Float32, tab[i, 3])
        if container === nothing
            radii[an] = ar
        else
            container[tab[i, 2]] = ar
        end
    end
    radii
end

"""
    ATOM_VANDERWAALS_RADII

Mapping of atomic numbers to van der Waals radii (in Å). See also [`atom_radius`](@ref).
"""
const ATOM_VANDERWAALS_RADII = let
    radiifile = joinpath(dirname(@__FILE__), "../../assets/const/vanderWaals_radii.csv")
    include_dependency(radiifile)
    tab, headers = readdlm(radiifile, '\t', String; header=true, comments=true)
    radii = Dict{Int,Float32}()
    for i = 1:size(tab, 1)
        an = parse(Int, tab[i, 1])
        ar = parse(Float32, tab[i, 2])
        radii[an] = ar
    end
    radii
end


"""
    AbstractAtom <: AbstractElement

The base class of vertex properties (atom).
"""
abstract type AbstractAtom <: AbstractElement end

"""
    StandardAtom <: AbstractAtom

The base class of atom types that have standard atom properties (e.g. [`SDFAtom`](@ref),
[`SMILESAtom`](@ref)). Built-in property functions (e.g. valence, mass) are implemented
for subtypes of this type.
"""
abstract type StandardAtom <: AbstractAtom end


has_isaromatic(::Type{T}) where T <: AbstractAtom = false
has_mol(::Type{T}) where T <: AbstractAtom = false
has_formula(::Type{T}) where T <: AbstractAtom = false
has_hydrogens(::Type{T}) where T <: AbstractAtom = false
has_label(::Type{T}) where T <: AbstractAtom = false


"""
    atom_number(atom::AbstractAtom) -> Int
    atom_number(atomsymbol::Symbol) -> Int

Return an atomic number of the given atom or the atomic symbol.
"""
atom_number(atom::AbstractAtom) = error("atom_number is not implemented for this atom type")
atom_number(atomsymbol::Symbol) = ATOMSYMBOLMAP[atomsymbol]


"""
    atom_symbol(atom::AbstractAtom) -> Symbol
    atom_symbol(n::Int) -> Symbol

Return an atomic symbol of the given atom or the atomic number.
"""
atom_symbol(atom::AbstractAtom) = error("atom_symbol is not implemented for this atom type")
atom_symbol(n::Int) = Symbol(ATOMTABLE[n].Symbol)


"""
    atom_valence(atom::AbstractAtom) -> Int

Return the valence specified to the atom, or `-1` if not specified (default). If not specified,
the valence is inferred from the typical valence of the element (see [`valence`](@ref)).

Atom types that can have a specified valence (e.g. carbenes and pentavalent nitrogens) should
implement this function. The fallback returns `-1`.
"""
atom_valence(atom::AbstractAtom) = -1


"""
    atom_num_hs(atom::AbstractAtom) -> Int

Return the number of hydrogens specified to the atom (e.g. SMILES bracket atoms and CDXML
`NumHydrogens`), or `-1` if not specified (default). The hydrogens include the hydrogen nodes
attached to the atom, and the atom has no implicit hydrogens other than the specified ones
(see [`implicit_hydrogens`](@ref)).

Atom types that can have a specified number of hydrogens should implement this function.
The fallback returns `-1`.
"""
atom_num_hs(atom::AbstractAtom) = -1


"""
    hydrogen_isotope(atomsymbol::Symbol) -> Int

Return the mass number implied by the atom symbol: `2` for deuterium (`:D`), `3` for
tritium (`:T`) and `0` for any other symbol.
"""
hydrogen_isotope(atomsymbol::Symbol) = atomsymbol === :D ? 2 : atomsymbol === :T ? 3 : 0

"""
    find_isotope(z::Int, num::Int) -> Union{Isotope, Nothing}
    find_isotope(atomsymbol::Symbol, num::Int) -> Union{Isotope, Nothing}

Return the isotope record (or nothing if there is no such isotope).
"""
function find_isotope(z::Int, num::Int)
    i = findfirst(x -> x.Number == num, ATOMTABLE[z].Isotopes)
    return isnothing(i) ? i : ATOMTABLE[z].Isotopes[i]
end

function find_isotope(atomsymbol::Symbol, num::Int)
    i = findfirst(x -> x.Number == num, ATOMTABLE[atom_number(atomsymbol)].Isotopes)
    return isnothing(i) ? i : ATOMTABLE[atom_number(atomsymbol)].Isotopes[i]
end



"""
    SDFAtom

SDFile (CTAB) atom property type.

`valence` is the valence specified to the atom (e.g. `2` for carbenes), or `-1` if not
specified (default). See [`atom_valence`](@ref) and [`valence`](@ref).
"""
struct SDFAtom <: StandardAtom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    isotope::Int
    coords::Union{Vector{Float64},Nothing}
    valence::Int

    function SDFAtom(
            symbol::Union{AbstractString,Symbol},
            charge::Int=0,
            multiplicity::Int=1,
            isotope::Int=0,
            coords::Union{Vector,Nothing}=nothing,
            valence::Int=-1)
        haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("sdfile parse error - unsupported atom symbol $(symbol)")
        if isotope != 0
            isnothing(find_isotope(Symbol(symbol), isotope)) && error("sdfile parse error - unsupported isotope $(isotope)$(symbol)")
        end
        valence >= -1 || error("sdfile parse error - invalid valence $(valence)")
        new(Symbol(symbol), charge, multiplicity, isotope, coords, valence)
    end
end

function SDFAtom(;
        symbol::Union{AbstractString,Symbol}=:C,
        charge::Int=0,
        multiplicity::Int=1,
        isotope::Int=0,
        coords::Union{Vector,Nothing}=nothing,
        valence::Int=-1)
    haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("sdfile parse error - unsupported atom symbol $(symbol)")
    return SDFAtom(Symbol(symbol), charge, multiplicity, isotope, coords, valence)
end

SDFAtom(d::Dict{Symbol,Any}) = SDFAtom(; NamedTuple((k, v) for (k, v) in d)...)

atom_symbol(a::SDFAtom) = a.symbol
atom_number(a::SDFAtom) = atom_number(a.symbol)
atom_charge(a::SDFAtom) = a.charge
multiplicity(a::SDFAtom) = a.multiplicity
isotope(a::SDFAtom) = a.isotope
atom_valence(a::SDFAtom) = a.valence

StructUtils.structlike(::Type{SDFAtom}) = false

# Concrete value types (instead of Any) so that JSON writing is compilable with `juliac --trim`
function StructUtils.lower(x::SDFAtom)
    rcd = Dict{String,Union{String,Int,Vector{Float64}}}()
    x.symbol === :C || setindex!(rcd, string(x.symbol), "symbol")
    x.charge == 0 || setindex!(rcd, x.charge, "charge")
    x.multiplicity == 1 || setindex!(rcd, x.multiplicity, "multiplicity")
    x.isotope == 0 || setindex!(rcd, x.isotope, "isotope")
    isnothing(x.coords) || setindex!(rcd, x.coords, "coords")
    x.valence == -1 || setindex!(rcd, x.valence, "valence")
    return rcd
end

StructUtils.lift(::Type{SDFAtom}, x) = SDFAtom(; NamedTuple((Symbol(k), v) for (k, v) in x)...)



"""
    SMILESAtom

SMILES atom property type.

`num_hs` is the number of hydrogens of a bracket atom (e.g. `0` for `[C]` and `3` for
`[CH3]`, including the hydrogen nodes), or `-1` for the atoms in the organic subset
(implicit hydrogens by the typical valence, default). See [`atom_num_hs`](@ref).
"""
struct SMILESAtom <: StandardAtom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    isotope::Int
    isaromatic::Bool
    stereo::Symbol
    num_hs::Int

    function SMILESAtom(
            symbol::Union{AbstractString,Symbol},
            charge::Int=0,
            multiplicity::Int=1,
            isotope::Int=0,
            isaromatic::Bool=false,
            stereo::Union{String,Symbol}=:unspecified,
            num_hs::Int=-1)
        haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("smiles parse error - unsupported atom symbol $(symbol)")
        if isotope != 0
            isnothing(find_isotope(Symbol(symbol), isotope)) && error("smiles parse error - unsupported isotope $(isotope)$(symbol)")
        end
        num_hs >= -1 || error("smiles parse error - invalid number of hydrogens $(num_hs)")
        new(Symbol(symbol), charge, multiplicity, isotope, isaromatic, Symbol(stereo), num_hs)
    end
end

function SMILESAtom(;
        symbol::Union{AbstractString,Symbol}=:C,
        charge::Int=0,
        multiplicity::Int=1,
        isotope::Int=0,
        isaromatic::Bool=false,
        stereo::Union{String,Symbol}=:unspecified,
        num_hs::Int=-1)
    return SMILESAtom(
        Symbol(symbol), charge, multiplicity, isotope, isaromatic, Symbol(stereo), num_hs)
end

SMILESAtom(d::Dict{Symbol,Any}) = SMILESAtom(; NamedTuple((k, v) for (k, v) in d)...)

has_isaromatic(::Type{SMILESAtom}) = true

atom_symbol(a::SMILESAtom) = a.symbol
atom_number(a::SMILESAtom) = atom_number(a.symbol)
atom_charge(a::SMILESAtom) = a.charge
multiplicity(a::SMILESAtom) = a.multiplicity
isotope(a::SMILESAtom) = a.isotope
isaromatic(a::SMILESAtom) = a.isaromatic
atom_num_hs(a::SMILESAtom) = a.num_hs

StructUtils.structlike(::Type{SMILESAtom}) = false

function StructUtils.lower(x::SMILESAtom)
    rcd = Dict{String,Union{String,Int,Bool}}()
    x.symbol === :C || setindex!(rcd, string(x.symbol), "symbol")
    x.charge == 0 || setindex!(rcd, x.charge, "charge")
    x.multiplicity == 1 || setindex!(rcd, x.multiplicity, "multiplicity")
    x.isotope == 0 || setindex!(rcd, x.isotope, "isotope")
    x.isaromatic === false || setindex!(rcd, x.isaromatic, "isaromatic")
    x.stereo === :unspecified || setindex!(rcd, string(x.stereo), "stereo")
    x.num_hs == -1 || setindex!(rcd, x.num_hs, "num_hs")
    return rcd
end

StructUtils.lift(::Type{SMILESAtom}, x) = SMILESAtom(; NamedTuple((Symbol(k), v) for (k, v) in x)...)



"""
    CommonChemAtom

CommonChem atom property type.
"""
struct CommonChemAtom <: StandardAtom
    z::Int
    chg::Int
    impHs::Int
    isotope::Int
    nRad::Int
    stereo::Symbol

    function CommonChemAtom(
            z::Int,
            chg::Int=0,
            impHs::Int=0,
            isotope::Int=0,
            nRad::Int=0,
            stereo::Union{AbstractString,Symbol}=:unspecified)
        z <= length(ATOMTABLE) || error("commonchem parse error - unsupported atom number $(z)")
        if isotope != 0
            isnothing(find_isotope(z, isotope)) && error("commonchem parse error - unsupported isotope $(isotope) for atom number $(z)")
        end
        new(z, chg, impHs, isotope, nRad, Symbol(stereo))
    end
end

function CommonChemAtom(;
        z::Int=6,
        chg::Int=0,
        impHs::Int=0,
        isotope::Int=0,
        nRad::Int=0,
        stereo::Union{AbstractString,Symbol}=:unspecified)
    z <= length(ATOMTABLE) || error("commonchem parse error - unsupported atom number $(z)")
    return CommonChemAtom(z, chg, impHs, isotope, nRad, Symbol(stereo))
end

CommonChemAtom(d::Dict{Symbol,Any}
    ) = CommonChemAtom(; NamedTuple((k, v) for (k, v) in d)...)
CommonChemAtom(x::SDFAtom) = CommonChemAtom(
    atom_number(x.symbol), x.charge, x.multiplicity - 1, x.isotope)
CommonChemAtom(x::SMILESAtom) = CommonChemAtom(
    atom_number(x.symbol), x.charge, x.multiplicity - 1, x.isotope)

atom_symbol(a::CommonChemAtom) = atom_symbol(a.z)
atom_number(a::CommonChemAtom) = a.z
atom_charge(a::CommonChemAtom) = a.chg
multiplicity(a::CommonChemAtom) = a.nRad + 1
isotope(a::CommonChemAtom) = a.isotope

StructUtils.structlike(::Type{CommonChemAtom}) = false

function StructUtils.lower(x::CommonChemAtom)
    rcd = Dict{String,Union{String,Int}}()
    x.z == 6 || setindex!(rcd, x.z, "z")
    x.chg == 0 || setindex!(rcd, x.chg, "chg")
    x.impHs == 0 || setindex!(rcd, x.impHs, "impHs")
    x.isotope == 0 || setindex!(rcd, x.isotope, "isotope")
    x.nRad == 0 || setindex!(rcd, x.nRad, "nRad")
    x.stereo === :unspecific || setindex!(rcd, string(x.stereo), "stereo")
    return rcd
end

StructUtils.lift(::Type{CommonChemAtom}, x) = CommonChemAtom(; NamedTuple((Symbol(k), v) for (k, v) in x)...)
