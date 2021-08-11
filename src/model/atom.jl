#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFileAtom, SmilesAtom,
    setcharge, setstereo, atomnumber, todict


const ATOMTABLE = let
    weightsfile = joinpath(dirname(@__FILE__), "../../assets/const/atomicweights.yaml")
    include_dependency(weightsfile)
    YAML.load(open(weightsfile))
end

const ATOMSYMBOLMAP = let
    symbolfile = joinpath(dirname(@__FILE__), "../../assets/const/symboltonumber.yaml")
    include_dependency(symbolfile)
    YAML.load(open(symbolfile))
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

SDFileAtom(data::Dict{String,Any}) = SDFileAtom(
    Symbol(data["symbol"]),
    data["charge"],
    data["multiplicity"],
    data["mass"],
    data["coords"],
    Symbol(data["stereo"])
)

function todict(a::SDFileAtom)
    data = Dict{String,Any}()
    for field in fieldnames(SDFileAtom)
        data[string(field)] = getfield(a, field)
    end
    return data
end

setcharge(a::SDFileAtom, chg
    ) = SDFileAtom(a.symbol, chg, a.multiplicity, a.mass, a.coords)

setstereo(a::SDFileAtom, direction) = SDFileAtom(
    a.symbol, a.charge, a.multiplicity, a.mass, a.coords, direction)



struct SmilesAtom <: Atom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    mass::Union{Int, Nothing}
    isaromatic::Union{Bool, Nothing}
    stereo::Symbol

    function SmilesAtom(sym, chg, multi, mass, aromatic, stereo)
        new(sym, chg, multi, mass, aromatic, stereo)
    end
end

SmilesAtom() = SmilesAtom(:C, 0, 1, nothing, nothing, :unspecified)
SmilesAtom(sym) = SmilesAtom(sym, 0, 1, nothing, nothing, :unspecified)
SmilesAtom(sym, chg) = SmilesAtom(sym, chg, 1, nothing, nothing, :unspecified)

SmilesAtom(data::Dict{String,Any}) = SmilesAtom(
    Symbol(data["symbol"]),
    data["charge"],
    data["multiplicity"],
    data["mass"],
    data["isaromatic"],
    Symbol(data["stereo"])
)

function todict(a::SmilesAtom)
    data = Dict{String,Any}()
    for field in fieldnames(SmilesAtom)
        data[string(field)] = getfield(a, field)
    end
    return data
end

setcharge(a::SmilesAtom, chg) = SmilesAtom(
    a.symbol, chg, a.multiplicity, a.mass, a.isaromatic, a.stereo)

setstereo(a::SmilesAtom, direction) = SmilesAtom(
    a.symbol, a.charge, a.multiplicity, a.mass, a.isaromatic, direction)


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
