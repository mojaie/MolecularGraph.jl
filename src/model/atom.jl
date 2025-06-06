#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const ATOMTABLE = let
    weightsfile = joinpath(dirname(@__FILE__), "../../assets/const/atomicweights.yaml")
    include_dependency(weightsfile)
    YAML.load(open(weightsfile))
end

const ATOMSYMBOLMAP = let
    symbolfile = joinpath(dirname(@__FILE__), "../../assets/const/symboltonumber.yaml")
    include_dependency(symbolfile)
    YAML.load(open(symbolfile); dicttype=Dict{String,Int})
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
atom_number(atomsymbol::Symbol) = ATOMSYMBOLMAP[string(atomsymbol)]


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

    function SDFAtom(sym=:C, chg=0, mul=1, ms=nothing, coords=nothing)
        haskey(ATOMSYMBOLMAP, string(sym)) || error("Unsupported atom symbol: $(sym)")
        new(sym, chg, mul, ms, coords)
    end
end

SDFAtom(d::Dict{T,Any}) where T <: Union{AbstractString,Symbol} = SDFAtom(
    Symbol(d[T("symbol")]), d[T("charge")], d[T("multiplicity")],
    d[T("mass")], d[T("coords")])
SDFAtom(arr::Vector) = SDFAtom(Symbol(arr[1]), arr[2], arr[3], arr[4], arr[5])

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

ELEMENT_TYPE_REGISTRY["SDFAtom"] = SDFAtom
to_dict(::Val{:default}, a::SDFAtom) = Any[a.symbol, a.charge, a.multiplicity, a.mass, a.coords]

function to_dict(::Val{:rdkit}, a::SDFAtom)
    rcd = Dict{String,Any}()
    a.symbol === :C || (rcd["z"] = atom_number(a.symbol))
    a.charge == 0 || (rcd["chg"] = a.charge)
    a.multiplicity == 1 || (rcd["nRad"] = a.multiplicity - 1)
    isnothing(a.mass) || (rcd["isotope"] = a.mass)
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
    isaromatic::Union{Bool,Nothing}
    stereo::Symbol

    function SMILESAtom(sym=:C, chg=0, mul=1, ms=nothing, isarom=false, stereo=:unspecified)
        new(sym, chg, mul, ms, isarom, stereo)
    end
end

# U may be necessary to determine whether T is String (for seriarization) or Symbol (for SMILES parser)
SMILESAtom(d::Dict{T,U}) where {T<:Union{AbstractString,Symbol},U} = SMILESAtom(
    Symbol(get(d, T("symbol"), :C)), get(d, T("charge"), 0),
    get(d, T("multiplicity"), 1), get(d, T("mass"), nothing),
    get(d, T("isaromatic"), false), Symbol(get(d, T("stereo"), :unspecified))
)
SMILESAtom(arr::Vector) = SMILESAtom(Symbol(arr[1]), arr[2], arr[3], arr[4], arr[5], Symbol(arr[6]))

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

ELEMENT_TYPE_REGISTRY["SMILESAtom"] = SMILESAtom
to_dict(::Val{:default}, a::SMILESAtom) = Any[a.symbol, a.charge, a.multiplicity, a.mass, a.isaromatic, a.stereo]

function to_dict(::Val{:rdkit}, a::SMILESAtom)
    rcd = Dict{String,Any}()
    a.symbol === :C || (rcd["z"] = atom_number(a.symbol))
    a.charge == 0 || (rcd["chg"] = a.charge)
    a.multiplicity == 1 || (rcd["nRad"] = a.multiplicity - 1)
    isnothing(a.mass) || (rcd["isotope"] = a.mass)
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
    mass::Int
    nRad::Int
    stereo::Symbol

    function CommonChemAtom(z=6, chg=0, impHs=0, mass=0, nRad=0, stereo=:unspecified)
        z <= length(ATOMTABLE) || error("Unsupported atom number: $(z)")
        new(z, chg, impHs, mass, nRad, stereo)
    end
end

CommonChemAtom(d::Dict{T,Any}) where T <: Union{AbstractString,Symbol} = CommonChemAtom(
    get(d, T("z"), 6), get(d, T("chg"), 0), get(d, T("impHs"), 0),
    get(d, T("mass"), 0), get(d, T("nRad"), 0), Symbol(get(d, T("stereo"), :unspecified))
)
CommonChemAtom(arr::Vector) = CommonChemAtom(
    arr[1], arr[2], arr[3], arr[4], arr[5], Symbol(arr[6]))

Base.getindex(a::CommonChemAtom, prop::Symbol) = getproperty(a, prop)
Base.:(==)(a1::CommonChemAtom, a2::CommonChemAtom) = all(
    [getfield(a1, f1) == getfield(a2, f2) for (f1, f2) in zip(fieldnames(typeof(a1)), fieldnames(typeof(a2)))])
Base.hash(a::CommonChemAtom, h::UInt
    ) = hash(a.z, hash(a.chg, hash(a.impHs, hash(a.mass, hash(a.nRad, hash(a.stereo, h))))))

atom_symbol(a::CommonChemAtom) = atom_symbol(a.z)
atom_number(a::CommonChemAtom) = a.z
atom_charge(a::CommonChemAtom) = a.chg
multiplicity(a::CommonChemAtom) = a.nRad + 1
atom_mass(a::CommonChemAtom) = a.mass == 0 ? nothing : a.mass

ELEMENT_TYPE_REGISTRY["CommonChemAtom"] = CommonChemAtom
to_dict(::Val{:default}, a::CommonChemAtom) = Any[a.z, a.chg, a.impHs, a.mass, a.nRad, a.stereo]

function to_dict(::Val{:rdkit}, a::CommonChemAtom)
    rcd = Dict{String,Any}()
    a.z == 6 || (rcd["z"] = a.z)
    a.chg == 0 || (rcd["chg"] = a.chg)
    a.impHs == 0 || (rcd["impHs"] = a.impHs)
    a.mass == 0 || (rcd["mass"] = a.mass)
    a.nRad == 0 || (rcd["nRad"] = a.nRad)
    a.stereo === :unspecific || (rcd["stereo"] = a.stereo)
    return rcd
end