#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


function sanitize_html(s::AbstractString)
    s = replace(s, "&" => "&amp;")
    s = replace(s, "<" => "&lt;")
    s = replace(s, ">" => "&gt;")
    s = replace(s, "\"" => "&quot;")
    s = replace(s, "'" => "&#39;")
    return s
end

function sanitize_markup(markup::Vector{Vector{Tuple{Symbol,String}}})
    san = Vector{Tuple{Symbol,String}}[]
    for gp in markup
        push!(san, [(sym, sanitize_html(unsafe_str)) for (sym, unsafe_str) in gp])
    end
    return san
end



struct VirtualAtom <: AbstractAtom
    label::Vector{Vector{Tuple{Symbol,String}}}

    function VirtualAtom(markup::Vector{Vector{Tuple{Symbol,String}}})
        new(sanitize_markup(markup))
    end
end

has_label(::Type{VirtualAtom}) = true

VirtualAtom(label::String) = VirtualAtom([[(:default, label)]])

atom_number(atom::VirtualAtom) = -1
atom_symbol(atom::VirtualAtom) = Symbol(write_formula(atom.label))
atom_charge(atom::VirtualAtom) = 0
multiplicity(atom::VirtualAtom) = 1
atom_color(atom::VirtualAtom; color_theme=DEFAULT_ATOM_COLOR, kwargs...) = color_theme[:C]
atom_markup(atom::VirtualAtom) = atom.label
atom_mass_unc(atom::VirtualAtom, f::F) where F = (0.0, 0.0)



struct HydrogenatedAtom{T<:StandardAtom} <: AbstractAtom
    center::T
    hydrogens::Int
    coords::Union{Vector{Float64},Nothing}
end

HydrogenatedAtom(atom::T, hs::Int; coords=nothing
    ) where T <: AbstractAtom = HydrogenatedAtom{T}(atom, hs, coords)

is_group(::Type{HydrogenatedAtom}) = true
has_hydrogens(::Type{HydrogenatedAtom}) = true
has_isaromatic(::Type{HydrogenatedAtom}) = true

atom_number(atom::HydrogenatedAtom) = atom_number(atom.center)
atom_symbol(atom::HydrogenatedAtom) = atom_symbol(atom.center)
atom_charge(atom::HydrogenatedAtom) = atom_charge(atom.center)
multiplicity(atom::HydrogenatedAtom) = multiplicity(atom.center)
isaromatic(atom::HydrogenatedAtom) = has_isaromatic(typeof(atom.center)) ? isaromatic(atom.center) : false
atom_color(atom::HydrogenatedAtom; color_theme=DEFAULT_ATOM_COLOR, kwargs...) = atom_color(atom.center)
atom_markup(atom::HydrogenatedAtom
    ) = [[(:default, string(atom_symbol(atom.center))), (:sub, string(atom.hydrogens))]]
atom_mass_unc(atom::HydrogenatedAtom, f::F
    ) where F = molecular_mass_unc(Dict(atom_symbol(atom.center) => 1, :H => atom.hydrogens), f)

atom_counter(atom::HydrogenatedAtom) = Dict(atom_symbol(atom.center) => 1, :H => atom.hydrogens)



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

is_group(::Type{FormulaGroup}) = true
has_label(::Type{FormulaGroup}) = true

atom_number(group::FormulaGroup) = -1
atom_symbol(group::FormulaGroup) = Symbol(write_formula(atom_markup(group)))
atom_charge(group::FormulaGroup) = group.charge
multiplicity(group::FormulaGroup) = 1
atom_color(group::FormulaGroup; color_theme=DEFAULT_ATOM_COLOR, kwargs...) = color_theme[:C]
atom_markup(group::FormulaGroup) = isempty(group.label) ? markup_formula(group.formula, charge=group.charge) : group.label
atom_mass_unc(group::FormulaGroup, f::F) where F = atom_mass_unc(group.formula, f)

atom_counter(group::FormulaGroup) = group.formula



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


is_group(::Type{<:StructGroup}) = true
has_label(::Type{<:StructGroup}) = true

atom_number(group::StructGroup) = -1
atom_symbol(group::StructGroup) = Symbol(write_formula(group.label))
atom_charge(group::StructGroup) = sum(atom_charge(group.mol))
multiplicity(group::StructGroup) = sum(multiplicity(group.mol))
atom_color(group::StructGroup; color_theme=DEFAULT_ATOM_COLOR, kwargs...) = color_theme[:C]
atom_markup(group::StructGroup) = group.label
atom_mass_unc(group::StructGroup, f::F) where F = molecular_mass_unc(group.mol, f)  # TODO: minus num conn * mass(:H)

atom_counter(group::StructGroup) = atom_counter(group.mol)



struct StructGroupBond{T<:Integer,E<:StandardBond} <: AbstractBond
    bond::E
    src::T  # -1 if source vertex is not StructGroup
    dst::T
    srcsub::Bool
    dstsub::Bool
end

StructGroupBond(
    bond::E; src::T=-1, dst::T=-1, srcsub::Bool=false, dstsub::Bool=false
) where {T,E} = StructGroupBond{T,E}(bond, src, dst, srcsub, dstsub)

is_group(::Type{<:StructGroupBond}) = true

bond_order(bond::StructGroupBond) = bond_order(bond.bond)

const GeneralMolGraph = MolGraph{Int,AbstractAtom,AbstractBond}


# TODO: functional groups (Ph, Bz, Ts ...)
methyl_group() = HydrogenatedAtom(SMILESAtom(:C), 3, "Me")
ethyl_group() = StructGroup(smilestomol("CC"), "Et")
tbutyl_group() = StructGroup(smilestomol("C(C)(C)C"), "tBu")

gly() = StructGroup(smilestomol("NCC(=O)O"), "Gly")
ala() = StructGroup(smilestomol("N[C@H](C(=O)O)C"), "Ala")
ser() = StructGroup(smilestomol("N[C@H](C(=O)O)CO"), "Ser")
cys() = StructGroup(smilestomol("N[C@H](C(=O)O)CS"), "Cys")
met() = StructGroup(smilestomol("N[C@H](C(=O)O)CCSC"), "Met")
lys() = StructGroup(smilestomol("N[C@H](C(=O)O)CCCCN"), "Lys")
val() = StructGroup(smilestomol("N[C@H](C(=O)O)C(C)C"), "Val")
thr() = StructGroup(smilestomol("N[C@H](C(=O)O)C(O)C"), "Thr")
ile() = StructGroup(smilestomol("N[C@H](C(=O)O)C(C)CC"), "Ile")
leu() = StructGroup(smilestomol("N[C@H](C(=O)O)CC(C)C"), "Leu")
pro() = StructGroup(smilestomol("N[C@H](C(=O)O)C1NCCC1"), "Pro")
asn() = StructGroup(smilestomol("N[C@H](C(=O)O)CC(=O)N"), "Asn")
asp() = StructGroup(smilestomol("N[C@H](C(=O)O)CC(=O)O"), "Asp")
gln() = StructGroup(smilestomol("N[C@H](C(=O)O)CCC(=O)N"), "Gln")
glu() = StructGroup(smilestomol("N[C@H](C(=O)O)CCC(=O)O"), "Glu")
phe() = StructGroup(smilestomol("N[C@H](C(=O)O)Cc1ccccc1"), "Phe")
arg() = StructGroup(smilestomol("N[C@H](C(=O)O)CCCNC(=N)N"), "Arg")
his() = StructGroup(smilestomol("N[C@H](C(=O)O)Cc1nc[nH]c1"), "His")
tyr() = StructGroup(smilestomol("N[C@H](C(=O)O)Cc1ccc(O)cc1"), "Tyr")
trp() = StructGroup(smilestomol("N[C@H](C(=O)O)Cc1cnc2ccccc12"), "Trp")


# TODO: sugar