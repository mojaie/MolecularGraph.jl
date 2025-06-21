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


# Atom interfaces

"""
    atom_number(atomsymbol::Symbol) -> Int
    atom_number(atomp::AbstractAtom) -> Int

Return an atomic number of the given atom or the atomic symbol.
"""
atom_number(atomsymbol::Symbol) = ATOMSYMBOLMAP[atomsymbol]
atom_number(atom::AbstractAtom) = error("atom_number is not implemented for this atom type")


"""
    atom_symbol(n::Int) -> Symbol
    atom_symbol(atom::AbstractAtom) -> Symbol

Return an atomic symbol of the given atom or the atomic number.
"""
atom_symbol(n::Int) = Symbol(ATOMTABLE[n]["Symbol"])
atom_symbol(atom::AbstractAtom) = error("atom_symbol is not implemented for this atom type")


"""
    atom_charge(atom::AbstractAtom) -> Int

Return atomic charge of the given atom.
"""
atom_charge(atom::AbstractAtom) = error("atom_charge is not implemented for this atom type")


"""
    multiplicity(atom::AbstractAtom) -> Int

Return multiplicity (num of radicals + 1) of the given atom.

This is experimental feature - free radical chemistry is still not introduced to this library.
This does nothing for now, but for example, you can set multiplicity=2 to molecular oxygens manually.
"""
multiplicity(atom::AbstractAtom) = error("multiplicity is not implemented for this atom type")


"""
    atom_mass(atom::AbstractAtom) -> Int

Return specific atomic mass of given atom, or return nothing if the mass is unspecified.
"""
atom_mass(atom::AbstractAtom) = error("atom_mass is not implemented for this atom type")



"""
    SDFAtom

SDFile (CTAB) atom property type.
"""
struct SDFAtom <: AbstractAtom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    mass::Union{Int,Nothing}
    coords::Union{Vector{Float64},Nothing}

    function SDFAtom(
            symbol::Union{AbstractString,Symbol},
            charge::Int=0,
            multiplicity::Int=1,
            mass::Union{Int,Nothing}=nothing,
            coords::Union{Vector,Nothing}=nothing)
        haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("Unsupported atom symbol: $(symbol)")
        new(Symbol(symbol), charge, multiplicity, mass, coords)
    end
end

function SDFAtom(;
        symbol::Union{AbstractString,Symbol}=:C,
        charge::Int=0,
        multiplicity::Int=1,
        mass::Union{Int,Nothing}=nothing,
        coords::Union{Vector,Nothing}=nothing)
    haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("Unsupported atom symbol: $(symbol)")
    return SDFAtom(Symbol(symbol), charge, multiplicity, mass, coords)
end

SDFAtom(d::Dict{String,Any}
    ) = SDFAtom(; NamedTuple((Symbol(k), v) for (k, v) in d)...)
SDFAtom(d::Dict{Symbol,Any}
    ) = SDFAtom(; NamedTuple((k, v) for (k, v) in d)...)

atom_symbol(a::SDFAtom) = a.symbol
atom_number(a::SDFAtom) = atom_number(a.symbol)
atom_charge(a::SDFAtom) = a.charge
multiplicity(a::SDFAtom) = a.multiplicity
atom_mass(a::SDFAtom) = a.mass

function to_dict(::Val{:default}, a::SDFAtom)
    rcd = Dict{String,Any}()
    a.symbol === :C || setindex!(rcd, string(a.symbol), "symbol")
    a.charge == 0 || setindex!(rcd, a.charge, "charge")
    a.multiplicity == 1 || setindex!(rcd, a.multiplicity, "multiplicity")
    isnothing(a.mass) || setindex!(rcd, a.mass, "mass")
    isnothing(a.coords) || setindex!(rcd, a.coords, "coords")
    return rcd
end

function to_dict(::Val{:rdkit}, a::SDFAtom)
    rcd = Dict{String,Any}()
    a.symbol === :C || setindex!(rcd, atom_number(a.symbol), "z")
    a.charge == 0 || setindex!(rcd, a.charge, "chg")
    a.multiplicity == 1 || setindex!(rcd, a.multiplicity - 1, "nRad")
    isnothing(a.mass) || setindex!(rcd, a.mass, "isotope")
    return rcd
end


"""
    SMILESAtom

SMILES atom property type.
"""
struct SMILESAtom <: AbstractAtom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    mass::Union{Int,Nothing}
    isaromatic::Bool
    stereo::Symbol

    function SMILESAtom(
            symbol::Union{AbstractString,Symbol},
            charge::Int=0,
            multiplicity::Int=1,
            mass::Union{Int,Nothing}=nothing,
            isaromatic::Bool=false,
            stereo::Union{String,Symbol}=:unspecified)
        new(Symbol(symbol), charge, multiplicity, mass, isaromatic, Symbol(stereo))
    end
end

function SMILESAtom(;
        symbol::Union{AbstractString,Symbol}=:C,
        charge::Int=0,
        multiplicity::Int=1,
        mass::Union{Int,Nothing}=nothing,
        isaromatic::Bool=false,
        stereo::Union{String,Symbol}=:unspecified)
    return SMILESAtom(
        Symbol(symbol), charge, multiplicity, mass, isaromatic, Symbol(stereo))
end

SMILESAtom(d::Dict{String,Any}
    ) = SMILESAtom(; NamedTuple((Symbol(k), v) for (k, v) in d)...)
SMILESAtom(d::Dict{Symbol,Any}
    ) = SMILESAtom(; NamedTuple((k, v) for (k, v) in d)...)

atom_symbol(a::SMILESAtom) = a.symbol
atom_number(a::SMILESAtom) = atom_number(a.symbol)
atom_charge(a::SMILESAtom) = a.charge
multiplicity(a::SMILESAtom) = a.multiplicity
atom_mass(a::SMILESAtom) = a.mass

function to_dict(::Val{:default}, a::SMILESAtom)
    rcd = Dict{String,Any}()
    a.symbol === :C || setindex!(rcd, string(a.symbol), "symbol")
    a.charge == 0 || setindex!(rcd, a.charge, "charge")
    a.multiplicity == 1 || setindex!(rcd, a.multiplicity, "multiplicity")
    isnothing(a.mass) || setindex!(rcd, a.mass, "mass")
    a.isaromatic === false || setindex!(rcd, a.isaromatic, "isaromatic")
    a.stereo === :unspecified || setindex!(rcd, string(a.stereo), "stereo")
    return rcd
end

function to_dict(::Val{:rdkit}, a::SMILESAtom)
    rcd = Dict{String,Any}()
    a.symbol === :C || setindex!(rcd, atom_number(a.symbol), "z")
    a.charge == 0 || setindex!(rcd, a.charge, "chg")
    a.multiplicity == 1 || setindex!(rcd, a.multiplicity - 1, "nRad")
    isnothing(a.mass) || setindex!(rcd, a.mass, "isotope")
    return rcd
end


"""
    CommonChemAtom

CommonChem atom property type.
"""
struct CommonChemAtom <: AbstractAtom
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
        z <= length(ATOMTABLE) || error("Unsupported atom number: $(z)")
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
    z <= length(ATOMTABLE) || error("Unsupported atom number: $(z)")
    return CommonChemAtom(z, chg, impHs, isotope, nRad, Symbol(stereo))
end

CommonChemAtom(d::Dict{String,Any}
    ) = CommonChemAtom(; NamedTuple((Symbol(k), v) for (k, v) in d)...)
CommonChemAtom(d::Dict{Symbol,Any}
    ) = CommonChemAtom(; NamedTuple((k, v) for (k, v) in d)...)

atom_symbol(a::CommonChemAtom) = atom_symbol(a.z)
atom_number(a::CommonChemAtom) = a.z
atom_charge(a::CommonChemAtom) = a.chg
multiplicity(a::CommonChemAtom) = a.nRad + 1
atom_mass(a::CommonChemAtom) = a.isotope == 0 ? nothing : a.mass

function to_dict(::Val, a::CommonChemAtom)
    rcd = Dict{String,Any}()
    a.z == 6 || setindex!(rcd, a.z, "z")
    a.chg == 0 || setindex!(rcd, a.chg, "chg")
    a.impHs == 0 || setindex!(rcd, a.impHs, "impHs")
    a.isotope == 0 || setindex!(rcd, a.isotope, "isotope")
    a.nRad == 0 || setindex!(rcd, a.nRad, "nRad")
    a.stereo === :unspecific || setindex!(rcd, a.stereo, "stereo")
    return rcd
end
