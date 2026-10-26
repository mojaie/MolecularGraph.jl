```@meta
CurrentModule = MolecularGraph
```

# Customize atom and bond properties

Atoms and bonds are the vertex and edge properties of `MolGraph{T,V,E}` (`V` and `E` are the
atom and bond types). The built-in readers use the following types.

| Reader | Atom type | Bond type |
|:--|:--|:--|
| [`sdfilereader`](@ref), [`sdftomol`](@ref) | [`SDFAtom`](@ref) | [`SDFBond`](@ref) |
| [`smilestomol`](@ref) | [`SMILESAtom`](@ref) | [`SMILESBond`](@ref) |
| [`mol_from_json`](@ref) (RDKit CommonChem) | [`CommonChemAtom`](@ref) | [`CommonChemBond`](@ref) |
| [`cdxmltomol`](@ref) | [`CDXMLAtom`](@ref) | [`CDXMLBond`](@ref) |
| [`smartstomol`](@ref) | [`QueryAtom`](@ref) | [`QueryBond`](@ref) |

Custom atom and bond types can be used to keep additional information on atoms and bonds.
`GeneralMolGraph` with virtual atoms (experimental) can represent abbreviated groups and
placeholders (e.g. R groups) that are not single atoms.


## Interfaces

Atom types are subtypes of [`AbstractAtom`](@ref) and should implement the following functions.

| Function | Description |
|:--|:--|
| [`atom_symbol`](@ref)`(atom)` | atomic symbol (`Symbol`, e.g. `:C`) |
| [`atom_number`](@ref)`(atom)` | atomic number |
| [`atom_charge`](@ref)`(atom)` | formal charge |
| [`multiplicity`](@ref)`(atom)` | spin multiplicity (1: no unpaired electron, 2: radical) |
| [`isotope`](@ref)`(atom)` | mass number (0: natural abundance) |
| [`atom_valence`](@ref)`(atom)` | (optional) specified valence, e.g. 2 for carbenes (the fallback returns `-1`: inferred from the element) |
| [`atom_num_hs`](@ref)`(atom)` | (optional) specified number of hydrogens including hydrogen nodes, e.g. 2 for `[CH2]` (the fallback returns `-1`: implicit hydrogens inferred from the valence) |

Bond types are subtypes of [`AbstractBond`](@ref) and should implement [`bond_order`](@ref)`(bond)`.

As subtypes of [`AbstractElement`](@ref), atoms and bonds are compared by all fields
(`Base.:(==)` and `Base.hash`), and fields can be accessed by `atom[:symbol]`. Atom and bond
types should be immutable. To change an atom or bond of a molecule, replace it by
`mol[i] = newatom` so that the descriptors are updated (see [MolecularGraph interfaces](@ref)).

Most of the property functions (e.g. [`valence`](@ref), [`molecular_formula`](@ref),
[`standard_weight`](@ref)) work with atom types that implement the functions above. Some functions
are implemented only for molecules whose atom and bond types are subtypes of
[`StandardAtom`](@ref) and [`StandardBond`](@ref) (e.g. [`printv2mol`](@ref), [`inchi`](@ref),
[`wclogp`](@ref), [`add_hydrogens!`](@ref)). Custom types of single atoms and ordinary bonds
should subtype them.

Optional features are declared by the following traits (non-exported functions of the type).

| Trait | Required function or fields | Used by |
|:--|:--|:--|
| `MolecularGraph.has_isaromatic(::Type) = true` | `MolecularGraph.isaromatic(atom)` | aromatic atoms in SMILES (kekulization) |
| `MolecularGraph.has_hydrogens(::Type) = true` | fields `center::StandardAtom`, `hydrogens::Int` | [`HydrogenatedAtom`](@ref) |
| `MolecularGraph.has_formula(::Type) = true` | field `formula::Dict{Symbol,Int}` | [`FormulaGroup`](@ref) |
| `MolecularGraph.has_label(::Type) = true` | field `label` (atom markup) | labels of virtual atoms in drawing |
| `MolecularGraph.has_mol(::Type) = true` | fields `mol`, `term` | [`StructGroup`](@ref) |
| `MolecularGraph.has_submap(::Type) = true` | fields `src`, `dst` | [`StructGroupBond`](@ref) |


## Custom atom and bond types with readers

Readers create atoms and bonds by the constructors of the atom and bond types of the molecule
type given as the first argument, and the default callbacks (see
[MolecularGraph interfaces](@ref)) read some fields of them. Therefore, custom types used with
the readers should have the following constructors and fields.

| Reader | Constructors | Fields read by the default callbacks |
|:--|:--|:--|
| SDFile | `V(symbol, charge, multiplicity, isotope, coords)`, `V(symbol, charge, multiplicity, isotope, coords, valence)` (only atoms with a specified valence; hydrogen counts `hhh`/`HCOUNT` are ignored except for [`CDXMLAtom`](@ref)), `E(order, notation, isordered)`, `E(; order, isordered)` (V3000) | atom: `coords`, bond: `order`, `notation`, `isordered` |
| SMILES | `V(; symbol, charge, isotope, isaromatic, stereo, num_hs)` (bracket atoms), `V(; symbol, isaromatic)` (organic subset), `V(; symbol=:H)`, `E(; order, isaromatic, direction)`, `E()` | atom: `stereo` (and `isaromatic` by the trait), bond: `order`, `direction` |

For example, an atom type that has an additional label can be used with `smilestomol` as follows.

```julia
using MolecularGraph, Graphs
import MolecularGraph: atom_symbol, atom_number, atom_charge, multiplicity, isotope

struct LabeledAtom <: MolecularGraph.StandardAtom
    symbol::Symbol
    charge::Int
    isotope::Int
    isaromatic::Bool
    stereo::Symbol
    label::String
end

# Constructor called by the SMILES parser (unsupported keywords are ignored)
LabeledAtom(; symbol=:C, charge=0, isotope=0, isaromatic=false, stereo=:unspecified,
    label="", kwargs...) = LabeledAtom(Symbol(symbol), charge, isotope, isaromatic, Symbol(stereo), label)

atom_symbol(a::LabeledAtom) = a.symbol
atom_number(a::LabeledAtom) = atom_number(a.symbol)
atom_charge(a::LabeledAtom) = a.charge
multiplicity(a::LabeledAtom) = 1
isotope(a::LabeledAtom) = a.isotope
MolecularGraph.has_isaromatic(::Type{LabeledAtom}) = true
MolecularGraph.isaromatic(a::LabeledAtom) = a.isaromatic

mol = smilestomol(MolGraph{Int,LabeledAtom,SMILESBond}, "CC(=O)[O-]")
mol[1] = LabeledAtom(symbol=:C, label="methyl")
mol[1].label  # "methyl"
molecular_formula(mol)  # "C2H3O2"
```

Molecules can also be constructed directly from an edge list and atom and bond lists. Give the
callbacks to calculate the descriptors that the property functions need.

```julia
mol = MolGraph(
    Edge.([(1, 2)]), [LabeledAtom(symbol=:C, label="a"), LabeledAtom(symbol=:O)], [SMILESBond(1)];
    on_init=smiles_on_init!, on_update=smiles_on_update!)
```

!!! note
    `JSON.json(mol)` writes molecules with custom types by the generic method of JSON.jl,
    but [`mol_from_json`](@ref) cannot read them (only the built-in molecule types are supported).


## Virtual atoms (experimental)

!!! warning "Experimental"
    `GeneralMolGraph` and the virtual atom types are experimental. The API and behavior may
    change without deprecation, and some functions do not support them yet (see the note below).

`GeneralMolGraph` (`MolGraph{Int,AbstractAtom,AbstractBond}`) can have different atom and bond
types in a molecule. The readers create the default atom types of the reader
(e.g. `smilestomol(GeneralMolGraph, smiles)` creates `SMILESAtom`s), and then the atoms can be
replaced by the following virtual atom types.

- [`VirtualAtom`](@ref): a placeholder such as an R group (e.g. `VirtualAtom("R")`). It has no
  element and is shown as the label in the drawing and the molecular formula.
- [`HydrogenatedAtom`](@ref): an atom with hydrogens as one vertex
  (e.g. `HydrogenatedAtom(SMILESAtom(:N), 1)` for NH).
- [`FormulaGroup`](@ref): a group given by the molecular formula and the charge
  (e.g. `FormulaGroup(Dict(:C => 2, :H => 7, :N => 1), 1, markup)`).
- [`StructGroup`](@ref): a group given by the structure with a label (e.g. amino acid
  residues). Bonds to `StructGroup`s should be [`StructGroupBond`](@ref)s that specify the
  terminal atoms of the group (`src` and `dst` are pairs of the terminal number and the vertex
  in the group).

```julia
mol = smilestomol(GeneralMolGraph, "c1ccccc1C")
mol[7] = VirtualAtom("R")
molecular_formula(mol)  # "C6H5R"

# Glu-Cys-Gly
mol = smilestomol(GeneralMolGraph, "N[C@H](C(=O)O)CCC(=O)")
add_vertex!(mol, MolecularGraph.cys())
add_vertex!(mol, MolecularGraph.gly())
add_edge!(mol, 9, 11, StructGroupBond(SMILESBond(1), dst=1=>1))
add_edge!(mol, 11, 12, StructGroupBond(SMILESBond(1), src=2=>7, dst=1=>1))
molecular_formula(mol)  # "C10H17N3O6S"
```

!!! note
    Functions that require `StandardAtom` and `StandardBond` (e.g. [`printv2mol`](@ref),
    [`inchi`](@ref), [`wclogp`](@ref)) are not available for `GeneralMolGraph`, and mass
    calculation fails if the molecule has `VirtualAtom`s. `GeneralMolGraph` has abstract
    atom and bond types, so it is slower than `MolGraph` with concrete types and is not
    compilable with `juliac --trim`.
