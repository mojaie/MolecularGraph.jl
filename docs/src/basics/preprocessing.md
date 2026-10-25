# Preprocessing

This page describes the following preprocessing strategies.

- Remove hydrogen vertices
- Extract molecules of interest
- Standardize charges
- Dealing with resonance structures
- Custom property updater

```@example preprocessing
using MolecularGraph, Graphs
```

Public databases (e.g. PubChem, ChEMBL) and flat file databases (e.g. SDFile) have different
conventions, and the molecules may not always be suitable for your analysis as they are.
For example,

- whether hydrogens are explicitly written or omitted
- whether salt and water molecules are included in the molecular graph or provided as metadata
- representation of resonance structures (e.g. diazo group; `[C-]-[N+]#N` <-> `C=[N+]=[N-]`)
- charges depend on the condition - powder, dissolved or in physiological conditions


## Remove hydrogen vertices

SDFiles downloaded from PubChem have hydrogen vertices. The SDFile used here can be obtained
from PubChem (e.g. `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6437877/SDF`
for cefditoren pivoxil).

```@example preprocessing
datadir = joinpath(pkgdir(MolecularGraph), "docs", "src", "assets", "data")
mol = sdftomol(joinpath(datadir, "Cefditoren Pivoxil.mol"))
```

In practice, hydrogens that are not important are often removed from molecular graphs for
simplicity.

- [`remove_hydrogens!`](@ref)`(mol)` removes hydrogen vertices that are not important (no
  charge, no unpaired electron, no specific isotope and not involved in stereochemistry).
- [`remove_all_hydrogens!`](@ref)`(mol)` removes all hydrogen vertices.

Both functions return the vmap (the original indices of the remaining vertices).

```@example preprocessing
vmap = remove_hydrogens!(mol)
nv(mol)
```

```@example preprocessing
mol
```

[`add_hydrogens!`](@ref) does the opposite; it adds all implicit hydrogens as vertices.


## Extract molecules of interest

`connected_components(mol)` (from Graphs.jl) returns the connected components, the sets of
vertices of the individual molecules in the molecular graph object.

```@example preprocessing
mol2 = smilestomol("CC(=O)OC1=CC=CC=C1C(=O)O.O.CCO")
html_fixed_size(mol2, 250, 250, atomindex=true)
```

```@example preprocessing
connected_components(mol2)
```

To extract the molecule of interest, you can iterate over the connected components, select
the component by your criteria and apply `induced_subgraph(mol, vertices)` to extract the
molecule. `induced_subgraph` returns the subgraph molecule and the vmap. For example, the
following code extracts the component that has the most oxygen atoms.

```@example preprocessing
symbols = atom_symbol(mol2)
comp = argmax(c -> count(==(:O), symbols[c]), connected_components(mol2))
sub, vmap = induced_subgraph(mol2, comp)
sub
```

Or simply [`extract_largest_component!`](@ref)`(mol)` can be used. This removes the vertices
that do not belong to the largest component (the connected component which has the largest
number of vertices) from the graph.

```@example preprocessing
let mol = copy(mol2)
    extract_largest_component!(mol)
    mol
end
```


## Standardize charges

- [`protonate_acids!`](@ref)`(mol)` removes charges on oxo/thio acid anions
- [`deprotonate_oniums!`](@ref)`(mol)` removes charges on ammonium/oxonium cations

```@example preprocessing
charged = smilestomol("CCCC(=O)[O-].[N+]CCCC[N+]")
```

```@example preprocessing
let mol = copy(charged)
    protonate_acids!(mol)
    deprotonate_oniums!(mol)
    mol
end
```

These functions do not modify the atom objects; they store the standardized charges as a
descriptor (`atom_charge(mol)` returns them). The non-mutating versions
([`protonate_acids`](@ref) and [`deprotonate_oniums`](@ref)) just return the new charge
vectors.

```@example preprocessing
println(protonate_acids(charged))
```

Note that these descriptors are recalculated from the atom objects when the molecule is
edited (see [Custom property updater](@ref) below for how to apply such preprocessing
automatically).


## Dealing with resonance structures

Substructure match methods in this library compare atom symbols and the number of ``\pi``
electrons, so in many cases you don't have to care about the variation of resonance
structures.

```@example preprocessing
quinoline1 = smilestomol("N1=CC=CC2=C1C=CC=C2")
quinoline2 = smilestomol("N=1C=CC=C2C1C=CC=C2")
html_grid([quinoline1, quinoline2], 2, 150)
```

```@example preprocessing
has_exact_match(quinoline1, quinoline2)
```

```@example preprocessing
println(pi_electron(quinoline1))
println(pi_electron(quinoline2))
```

```@example preprocessing
inchikey(quinoline1) == inchikey(quinoline2)
```

For functional groups that have multiple representations, the following functions convert
them into one of the representations. As with the charge standardization, the non-mutating
versions return a tuple of the atom charges and the bond orders, and the `!` versions
store them as descriptors.

- [`polarize!`](@ref): polarizes dipole double bonds (e.g. `S=O` -> `[S+]-[O-]`)
- [`depolarize!`](@ref): depolarizes dipole double bonds (e.g. `[C+]-[O-]` -> `C=O`)
- [`to_triple_bond!`](@ref): 1,3-dipoles as triple bonds (e.g. `C=[N+]=[N-]` -> `[C-]-[N+]#N`)
- [`to_allene_like!`](@ref): 1,3-dipoles as allene-like double bonds (e.g. `[C-]-[N+]#N` ->
  `C=[N+]=[N-]`)

```@example preprocessing
dmso = smilestomol("CS(=O)C")
diazo = smilestomol("CC(=O)C=[N+]=[N-]")
dmso_p = copy(dmso)
polarize!(dmso_p)
diazo_t = copy(diazo)
to_triple_bond!(diazo_t)
html_grid([dmso, dmso_p, diazo, diazo_t], 4, 150)
```

Kekulization assigns double bonds and single bonds to aromatic rings which consist of SMILES
lowercase atoms. [`kekulize`](@ref) returns a tuple of the bond orders and the pyrrole-like
aromatic atoms (e.g. `[nH]`), and [`kekulize!`](@ref) stores them in the molecule.
Kekulization is applied automatically to molecules parsed from SMILES by the property
auto-update mechanism (see below). To see how kekulization works, here we construct a
molecule with a custom updater that skips kekulization.

```@example preprocessing
function custom_updater!(mol)
    # preprocessing
    default_atom_charge!(mol)
    default_bond_order!(mol)
    # kekulize!(mol)  # skip kekulize!
    # recalculate bottleneck descriptors
    sssr!(mol)
    apparent_valence!(mol)
    valence!(mol)
    lone_pair!(mol)
    is_ring_aromatic!(mol)
end

pyrene = smilestomol("c1cc2cccc3c2c4c1cccc4cc3"; on_update=custom_updater!)
println(bond_order(pyrene))
```

```@example preprocessing
kekulize!(pyrene)
println(bond_order(pyrene))
```

```@example preprocessing
pyrene
```


## Custom property updater

As described in [Properties and descriptors](@ref), the descriptors of a molecule are
recalculated by the `on_update` callback when the molecule has been edited. The following
code is the default property updater for SMILES molecules ([`smiles_on_update!`](@ref)).

```julia
function smiles_on_update!(mol)
    # preprocessing
    default_atom_charge!(mol)
    default_bond_order!(mol)
    kekulize!(mol)  # also sets apparent_valence
    # recalculate bottleneck descriptors
    sssr!(mol)
    valence!(mol)
    lone_pair!(mol)
    is_ring_aromatic!(mol)
end
```

- `default_atom_charge!` and `default_bond_order!` are necessary to set atom
  charges and bond orders as descriptors. You can add custom preprocessing functions
  (e.g. [`protonate_acids!`](@ref)) after them to standardize charges and notation of
  resonance structures.
- The lines after that cache bottleneck descriptors. These are essential for most of the
  other descriptors, but you can remove them if you don't use descriptors at all
  (e.g. in molecular graph generation tasks).

The callbacks (`on_init` and `on_update`) are part of the molecule type and fixed when the
molecule is constructed, so they cannot be replaced by assignment (`mol.state.on_update = f`
throws an error). Give the custom updater to the molecule readers, or create a molecule with
other callbacks by `MolecularGraph.with_callbacks`.

```julia
mol = smilestomol(smiles; on_update=my_updater!)
mol = sdftomol(file; config=(on_init=sdf_on_init!, on_update=my_updater!))
mol = mol_from_json(json; on_update=my_updater!)
mol2 = MolecularGraph.with_callbacks(mol; on_update=my_updater!)
```

For example, the following updater always neutralizes acids and oniums, even after the
molecule is edited.

```@example preprocessing
function neutralizing_updater!(mol)
    smiles_on_update!(mol)
    protonate_acids!(mol)
    deprotonate_oniums!(mol)
end

neutral = MolecularGraph.with_callbacks(copy(charged); on_update=neutralizing_updater!)
add_vertex!(neutral, SMILESAtom(symbol=:C))  # edit the molecule to trigger the update
add_edge!(neutral, 7, nv(neutral), SMILESBond())
neutral
```

See [MolecularGraph interfaces](@ref) for details of the callback mechanism.
