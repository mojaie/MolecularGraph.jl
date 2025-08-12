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

    function StructGroup{T,V,E}(
                mol::MolGraph{T,V,E},
                label::Vector{Vector{Tuple{Symbol,String}}}) where {T,V,E}
        new(mol, sanitize_markup(label))
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



struct StructGroupBond{T<:Integer,E<:StandardBond} <: AbstractBond
    bond::E
    src::T  # -1 if source vertex is not StructGroup
    dst::T
end

StructGroupBond(bond::E; src::T=-1, dst::T=-1
    ) where {T,E} = StructGroupBond{T,E}(bond, src, dst)

has_submap(::Type{<:StructGroupBond}) = true

bond_order(bond::StructGroupBond) = bond_order(bond.bond)


# TODO: functional groups (Ph, Bz, Ts ...)
# TODO: sugar
methyl_group() = HydrogenatedAtom(SMILESAtom(:C), 3, "Me")
ethyl_group() = StructGroup(smilestomol("CC"), "Et")
tbutyl_group() = StructGroup(smilestomol("C(C)(C)C"), "tBu")

gly() = StructGroup(smilestomol("NCC=O"), "Gly")
ala() = StructGroup(smilestomol("N[C@H](C=O)C"), "Ala")
ser() = StructGroup(smilestomol("N[C@H](C=O)CO"), "Ser")
cys() = StructGroup(smilestomol("N[C@H](C=O)CS"), "Cys")
met() = StructGroup(smilestomol("N[C@H](C=O)CCSC"), "Met")
lys() = StructGroup(smilestomol("N[C@H](C=O)CCCCN"), "Lys")
val() = StructGroup(smilestomol("N[C@H](C=O)C(C)C"), "Val")
thr() = StructGroup(smilestomol("N[C@H](C=O)C(O)C"), "Thr")
ile() = StructGroup(smilestomol("N[C@H](C=O)C(C)CC"), "Ile")
leu() = StructGroup(smilestomol("N[C@H](C=O)CC(C)C"), "Leu")
pro() = StructGroup(smilestomol("N[C@H](C=O)C1NCCC1"), "Pro")
asn() = StructGroup(smilestomol("N[C@H](C=O)CC(=O)N"), "Asn")
asp() = StructGroup(smilestomol("N[C@H](C=O)CC(=O)O"), "Asp")
gln() = StructGroup(smilestomol("N[C@H](C=O)CCC(=O)N"), "Gln")
glu() = StructGroup(smilestomol("N[C@H](C=O)CCC(=O)O"), "Glu")
phe() = StructGroup(smilestomol("N[C@H](C=O)Cc1ccccc1"), "Phe")
arg() = StructGroup(smilestomol("N[C@H](C=O)CCCNC(=N)N"), "Arg")
his() = StructGroup(smilestomol("N[C@H](C=O)Cc1nc[nH]c1"), "His")
tyr() = StructGroup(smilestomol("N[C@H](C=O)Cc1ccc(O)cc1"), "Tyr")
trp() = StructGroup(smilestomol("N[C@H](C=O)Cc1cnc2ccccc12"), "Trp")
