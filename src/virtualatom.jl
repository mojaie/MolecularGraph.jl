#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


const GeneralMolGraph = MolGraph{Int,AbstractAtom,AbstractBond}


struct VirtualAtom <: AbstractAtom
    label::Vector{Vector{Tuple{Symbol,String}}}

    function VirtualAtom(markup::Vector{Vector{Tuple{Symbol,String}}})
        new(sanitize_markup(markup))
    end
end

VirtualAtom(label::String) = VirtualAtom([[(:default, label)]])

has_label(::Type{VirtualAtom}) = true

atom_number(atom::VirtualAtom) = -1
atom_symbol(atom::VirtualAtom) = Symbol(write_formula(atom.label))
atom_charge(atom::VirtualAtom) = 0
multiplicity(atom::VirtualAtom) = 1
isotope(atom::VirtualAtom) = 0


struct HydrogenatedAtom{T<:StandardAtom} <: AbstractAtom
    center::T
    hydrogens::Int
    coords::Union{Vector{Float64},Nothing}
end

HydrogenatedAtom(atom::T, hs::Int; coords=nothing
    ) where T <: AbstractAtom = HydrogenatedAtom{T}(atom, hs, coords)

has_hydrogens(::Type{T}) where T <: HydrogenatedAtom = true
has_isaromatic(::Type{T}) where T <: HydrogenatedAtom = true

atom_number(atom::HydrogenatedAtom) = atom_number(atom.center)
atom_symbol(atom::HydrogenatedAtom) = atom_symbol(atom.center)
atom_charge(atom::HydrogenatedAtom) = atom_charge(atom.center)
multiplicity(atom::HydrogenatedAtom) = multiplicity(atom.center)
isotope(atom::HydrogenatedAtom) = isotope(atom.center)
isaromatic(atom::HydrogenatedAtom) = has_isaromatic(typeof(atom.center)) ? isaromatic(atom.center) : false



struct FormulaGroup <: AbstractAtom
    formula::Dict{Symbol,Int}
    charge::Int
    label::Vector{Vector{Tuple{Symbol,String}}}

    function FormulaGroup(
            formula::Dict{Symbol,Int},
            charge::Int,
            markup::Vector{Vector{Tuple{Symbol,String}}})
        new(formula, charge, sanitize_markup(markup))
    end
end

has_formula(::Type{FormulaGroup}) = true
has_label(::Type{FormulaGroup}) = true

atom_number(group::FormulaGroup) = -1
atom_symbol(group::FormulaGroup) = Symbol(write_formula(atom_markup(group)))
atom_charge(group::FormulaGroup) = group.charge
multiplicity(group::FormulaGroup) = 1
isotope(group::FormulaGroup) = 0



struct StructGroup{T<:Integer,V<:StandardAtom,E<:StandardBond} <: AbstractAtom
    mol::MolGraph{T,V,E}
    label::Vector{Vector{Tuple{Symbol,String}}}
    term::Dict{Int,T}  # terminal No. => vertex in StructGroup

    function StructGroup{T,V,E}(
            mol::MolGraph{T,V,E},
            label::Vector{Vector{Tuple{Symbol,String}}}) where {T,V,E}
        new(mol, sanitize_markup(label), Dict{Int,T}())
    end
end

StructGroup(mol::MolGraph{T,V,E}, label::String
    ) where {T,V,E} = StructGroup{T,V,E}(mol, [[(:default, label)]])

has_mol(::Type{<:StructGroup}) = true
has_label(::Type{<:StructGroup}) = true

atom_number(group::StructGroup) = -1
atom_symbol(group::StructGroup) = Symbol(write_formula(group.label))
atom_charge(group::StructGroup) = sum(atom_charge(group.mol))
multiplicity(group::StructGroup) = 1
isotope(group::StructGroup) = 0



struct StructGroupBond{T<:StandardBond} <: AbstractBond
    bond::T
    src::Pair{Int,Int}  # term label (-1 if source is not StructGroup) => term vertex in StructGroup
    dst::Pair{Int,Int}
end

StructGroupBond(bond::T; src::Pair{Int,Int}=-1=>-1, dst::Pair{Int,Int}=-1=>-1
    ) where T = StructGroupBond{T}(bond, src, dst)

has_submap(::Type{<:StructGroupBond}) = true

bond_order(bond::StructGroupBond) = bond_order(bond.bond)

function Graphs.add_edge!(mol::SimpleMolGraph, u::Integer, v::Integer, prop::StructGroupBond)
    e = u_edge(mol, u, v)
    add_u_edge!(mol, e, prop)
    src = mol[e.src]
    if has_mol(typeof(src))
        prop.src[1] > 0 || error("src vertex should be specified")
        src.term[prop.src[1]] = prop.src[2]
    end
    dst = mol[e.dst]
    if has_mol(typeof(dst))
        prop.dst[1] > 0 || error("dst vertex should be specified")
        dst.term[prop.dst[1]] = prop.dst[2]
    end
end


# TODO: functional groups (Ph, Bz, Ts ...)
# TODO: sugar
methyl_group() = HydrogenatedAtom(SMILESAtom(:C), 3, "Me")
ethyl_group() = StructGroup(smilestomol("[H]CC"), "Et")  # term => 1
tbutyl_group() = StructGroup(smilestomol("[H]C(C)(C)C"), "tBu")  # term => 1

gly() = StructGroup(smilestomol("[H]NC(C=O)O"), "Gly")  # Cterm => 1, Nterm => 6
ala() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)C"), "Ala")  # Cterm => 1, Nterm => 7
ser() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CO"), "Ser")
cys() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CS"), "Cys")
met() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CCSC"), "Met")
lys() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CCCCN"), "Lys")
val() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)C(C)C"), "Val")
thr() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)C(O)C"), "Thr")
ile() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)C(C)CC"), "Ile")
leu() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CC(C)C"), "Leu")
pro() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)C1NCCC1"), "Pro")
asn() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CC(=O)N"), "Asn")
asp() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CC(=O)O"), "Asp")
gln() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CCC(=O)N"), "Gln")
glu() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CCC(=O)O"), "Glu")
phe() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)Cc1ccccc1"), "Phe")
arg() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)CCCNC(=N)N"), "Arg")
his() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)Cc1nc[nH]c1"), "His")
tyr() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)Cc1ccc(O)cc1"), "Tyr")
trp() = StructGroup(smilestomol("[H]N[C@H](C(=O)O)Cc1cnc2ccccc12"), "Trp")
