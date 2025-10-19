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

const ATOMSYMBOLMAP = let
    symbolfile = joinpath(dirname(@__FILE__), "../../assets/const/symboltonumber.yaml")
    include_dependency(symbolfile)
    YAML.load(open(symbolfile); dicttype=Dict{Symbol,Int})
end

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
atom_symbol(n::Int) = Symbol(ATOMTABLE[n]["Symbol"])


"""
    find_isotope(atomsymbol::Symbol, num::Int) -> Union{Isotope, Nothing}

Return the isotope record (or nothing if there is no such isotope).
"""
function find_isotope(z::Int, num::Int)
    i = findfirst(x -> x["Number"] == num, ATOMTABLE[z]["Isotopes"])
    return isnothing(i) ? i : ATOMTABLE[z]["Isotopes"][i]
end

function find_isotope(atomsymbol::Symbol, num::Int)
    i = findfirst(x -> x["Number"] == num, ATOMTABLE[atom_number(atomsymbol)]["Isotopes"])
    return isnothing(i) ? i : ATOMTABLE[atom_number(atomsymbol)]["Isotopes"][i]
end



"""
    SDFAtom

SDFile (CTAB) atom property type.
"""
struct SDFAtom <: StandardAtom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    isotope::Int
    coords::Union{Vector{Float64},Nothing}

    function SDFAtom(
            symbol::Union{AbstractString,Symbol},
            charge::Int=0,
            multiplicity::Int=1,
            isotope::Int=0,
            coords::Union{Vector,Nothing}=nothing)
        haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("sdfile parse error - unsupported atom symbol $(symbol)")
        if isotope != 0
            isnothing(find_isotope(Symbol(symbol), isotope)) && error("sdfile parse error - unsupported isotope $(isotope)$(symbol)")
        end
        new(Symbol(symbol), charge, multiplicity, isotope, coords)
    end
end

function SDFAtom(;
        symbol::Union{AbstractString,Symbol}=:C,
        charge::Int=0,
        multiplicity::Int=1,
        isotope::Int=0,
        coords::Union{Vector,Nothing}=nothing)
    haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("sdfile parse error - unsupported atom symbol $(symbol)")
    return SDFAtom(Symbol(symbol), charge, multiplicity, isotope, coords)
end

SDFAtom(d::Dict{Symbol,Any}) = SDFAtom(; NamedTuple((k, v) for (k, v) in d)...)

atom_symbol(a::SDFAtom) = a.symbol
atom_number(a::SDFAtom) = atom_number(a.symbol)
atom_charge(a::SDFAtom) = a.charge
multiplicity(a::SDFAtom) = a.multiplicity
isotope(a::SDFAtom) = a.isotope

StructUtils.structlike(::StructUtils.StructStyle, ::Type{SDFAtom}) = false

function JSON.lower(x::SDFAtom)
    rcd = Dict{String,Any}()
    x.symbol === :C || setindex!(rcd, string(x.symbol), "symbol")
    x.charge == 0 || setindex!(rcd, x.charge, "charge")
    x.multiplicity == 1 || setindex!(rcd, x.multiplicity, "multiplicity")
    x.isotope == 0 || setindex!(rcd, x.isotope, "isotope")
    isnothing(x.coords) || setindex!(rcd, x.coords, "coords")
    return rcd
end

JSON.lift(::Type{SDFAtom}, x) = SDFAtom(; NamedTuple((Symbol(k), v) for (k, v) in x)...)



"""
    SMILESAtom

SMILES atom property type.
"""
struct SMILESAtom <: StandardAtom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    isotope::Int
    isaromatic::Bool
    stereo::Symbol

    function SMILESAtom(
            symbol::Union{AbstractString,Symbol},
            charge::Int=0,
            multiplicity::Int=1,
            isotope::Int=0,
            isaromatic::Bool=false,
            stereo::Union{String,Symbol}=:unspecified)
        haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("smiles parse error - unsupported atom symbol $(symbol)")
        if isotope != 0
            isnothing(find_isotope(Symbol(symbol), isotope)) && error("smiles parse error - unsupported isotope $(isotope)$(symbol)")
        end
        new(Symbol(symbol), charge, multiplicity, isotope, isaromatic, Symbol(stereo))
    end
end

function SMILESAtom(;
        symbol::Union{AbstractString,Symbol}=:C,
        charge::Int=0,
        multiplicity::Int=1,
        isotope::Int=0,
        isaromatic::Bool=false,
        stereo::Union{String,Symbol}=:unspecified)
    return SMILESAtom(
        Symbol(symbol), charge, multiplicity, isotope, isaromatic, Symbol(stereo))
end

SMILESAtom(d::Dict{Symbol,Any}) = SMILESAtom(; NamedTuple((k, v) for (k, v) in d)...)

has_isaromatic(::Type{SMILESAtom}) = true

atom_symbol(a::SMILESAtom) = a.symbol
atom_number(a::SMILESAtom) = atom_number(a.symbol)
atom_charge(a::SMILESAtom) = a.charge
multiplicity(a::SMILESAtom) = a.multiplicity
isotope(a::SMILESAtom) = a.isotope
isaromatic(a::SMILESAtom) = a.isaromatic

StructUtils.structlike(::StructUtils.StructStyle, ::Type{SMILESAtom}) = false

function JSON.lower(x::SMILESAtom)
    rcd = Dict{String,Any}()
    x.symbol === :C || setindex!(rcd, string(x.symbol), "symbol")
    x.charge == 0 || setindex!(rcd, x.charge, "charge")
    x.multiplicity == 1 || setindex!(rcd, x.multiplicity, "multiplicity")
    x.isotope == 0 || setindex!(rcd, x.isotope, "isotope")
    x.isaromatic === false || setindex!(rcd, x.isaromatic, "isaromatic")
    x.stereo === :unspecified || setindex!(rcd, string(x.stereo), "stereo")
    return rcd
end

JSON.lift(::Type{SMILESAtom}, x) = SMILESAtom(; NamedTuple((Symbol(k), v) for (k, v) in x)...)



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

StructUtils.structlike(::StructUtils.StructStyle, ::Type{CommonChemAtom}) = false

function JSON.lower(x::CommonChemAtom)
    rcd = Dict{String,Any}()
    x.z == 6 || setindex!(rcd, x.z, "z")
    x.chg == 0 || setindex!(rcd, x.chg, "chg")
    x.impHs == 0 || setindex!(rcd, x.impHs, "impHs")
    x.isotope == 0 || setindex!(rcd, x.isotope, "isotope")
    x.nRad == 0 || setindex!(rcd, x.nRad, "nRad")
    x.stereo === :unspecific || setindex!(rcd, x.stereo, "stereo")
    return rcd
end

JSON.lift(::Type{CommonChemAtom}, x) = CommonChemAtom(; NamedTuple((Symbol(k), v) for (k, v) in x)...)
