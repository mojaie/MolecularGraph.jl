#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

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

atom_symbol(a::AbstractDict) = a[:symbol]
atom_number(a::AbstractDict) = atom_number(a[:symbol])
atom_charge(a::AbstractDict) = a[:charge]
multiplicity(a::AbstractDict) = a[:multiplicity]
atom_mass(a::AbstractDict) = a[:mass]


"""
    atom_number(atomsymbol::Symbol) -> Int
    atom_number(atomprop::SDFAtom) -> Int
    atom_number(atomprop::SMILESAtom) -> Int
    atom_number(atomprop::CommonChemAtom) -> Int

Return atom number.
"""
atom_number(atomsymbol::Symbol) = ATOMSYMBOLMAP[atomsymbol]


"""
    atom_symbol(n::Int) -> Symbol
    atom_symbol(atomprop::SDFAtom) -> Symbol
    atom_symbol(atomprop::SMILESAtom) -> Symbol
    atom_symbol(atomprop::CommonChemAtom) -> Symbol

Return atom symbol of given atomic number.
"""
atom_symbol(n::Int) = Symbol(ATOMTABLE[n]["Symbol"])


"""
    SDFAtom

SDFile (CTAB) atom property type.
"""
struct SDFAtom
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

Base.getindex(a::SDFAtom, prop::Symbol) = getproperty(a, prop)
Base.:(==)(a1::SDFAtom, a2::SDFAtom) = all(
    [getfield(a1, f1) == getfield(a2, f2) for (f1, f2) in zip(fieldnames(typeof(a1)), fieldnames(typeof(a2)))])
Base.hash(a::SDFAtom, h::UInt
    ) = hash(a.symbol, hash(a.charge, hash(a.multiplicity, hash(a.mass, hash(a.coords, h)))))

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
struct SMILESAtom
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

Base.getindex(a::SMILESAtom, prop::Symbol) = getproperty(a, prop)
Base.:(==)(a1::SMILESAtom, a2::SMILESAtom) = all(
    [getfield(a1, f1) == getfield(a2, f2) for (f1, f2) in zip(fieldnames(typeof(a1)), fieldnames(typeof(a2)))])
Base.hash(a::SMILESAtom, h::UInt
    ) = hash(a.symbol, hash(a.charge, hash(a.multiplicity, hash(a.mass, hash(a.isaromatic, hash(a.stereo, h))))))

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
struct CommonChemAtom
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

Base.getindex(a::CommonChemAtom, prop::Symbol) = getproperty(a, prop)
Base.:(==)(a1::CommonChemAtom, a2::CommonChemAtom) = all(
    [getfield(a1, f1) == getfield(a2, f2) for (f1, f2) in zip(fieldnames(typeof(a1)), fieldnames(typeof(a2)))])
Base.hash(a::CommonChemAtom, h::UInt
    ) = hash(a.z, hash(a.chg, hash(a.impHs, hash(a.mass, hash(a.nRad, hash(a.stereo, h))))))

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