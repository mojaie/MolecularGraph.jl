# Properties and descriptors

This page describes the built-in molecular properties and descriptors, and how they are
kept consistent with the molecular graph.

- Built-in molecule properties and descriptors
  - Lipinski's rule of five (RO5)
  - Molecular formula
  - Atom and bond properties
  - Graph topology (rings and fused rings)
- Auto-update mechanism of properties

```@example properties
using MolecularGraph, Graphs
```


## Built-in molecule properties

Here are examples of built-in properties often used in molecular data mining.

### Lipinski's rule of five (RO5)

Methods that return a `Float64` value accept the number of `digits` (after the decimal point)
as an optional argument.

- [`standard_weight`](@ref)`(mol[, digits])`: molecular weight
- [`hydrogen_acceptor_count`](@ref)`(mol)`: number of hydrogen bond acceptors
- [`hydrogen_donor_count`](@ref)`(mol)`: number of hydrogen bond donors
- [`wclogp`](@ref)`(mol[, digits])`: predicted logP by Wildman and Crippen (1999),
  [https://doi.org/10.1021/ci990307l](https://doi.org/10.1021/ci990307l)
- [`rotatable_count`](@ref)`(mol)`: number of rotatable bonds (used in some RO5 variants)

```@example properties
mol = smilestomol("CC(=O)OC1=CC=CC=C1C(=O)O")
```

```@example properties
standard_weight(mol, 2)
```

```@example properties
hydrogen_acceptor_count(mol)
```

```@example properties
hydrogen_donor_count(mol)
```

```@example properties
wclogp(mol, 2)
```

```@example properties
rotatable_count(mol)
```

See [Mass and isotopes](@ref) for other mass-related functions (e.g. monoisotopic mass and
exact mass).

### Molecular formula

- [`molecular_formula`](@ref)`(mol)`: molecular formula in the Hill system
- [`empirical_formula`](@ref)`(mol)`: empirical formula in the Hill system

```@example properties
oleic_acid = smilestomol(raw"CCCCCCCC\C=C/CCCCCCCC(O)=O")
```

```@example properties
molecular_formula(oleic_acid)
```

```@example properties
empirical_formula(oleic_acid)
```

[`atom_counter`](@ref) returns the number of atoms of each element (including implicit
hydrogens), and [`heavy_atom_count`](@ref) returns the number of non-hydrogen atoms.

```@example properties
atom_counter(oleic_acid)
```

```@example properties
heavy_atom_count(oleic_acid)
```

### Atom and bond properties

The following methods return a vector of atom properties in vertex index order.

- [`atom_symbol`](@ref)`(mol)`: atom letters as symbols, e.g. `:C`, `:O` and `:N`
- [`atom_charge`](@ref)`(mol)`: electric charge of the atom. Only integer charges are
  allowed in the model.
- [`multiplicity`](@ref)`(mol)`: 1: no unpaired electron (default), 2: radical, 3: biradical
- [`lone_pair`](@ref)`(mol)`: number of lone pairs on the atom
- [`implicit_hydrogens`](@ref)`(mol)`: number of implicit hydrogens, that do not appear as
  graph vertices but are automatically calculated, drawn in images and used for the
  calculation of other descriptors. [`explicit_hydrogens`](@ref) and [`total_hydrogens`](@ref)
  are also available.
- [`valence`](@ref)`(mol)`: atom valence, specific to each atom species and considering the
  electric charge. The number of implicit hydrogens is obtained by subtracting the degree of
  the vertex from the valence.
- [`is_aromatic`](@ref)`(mol)`: whether the atom is aromatic or not. Only binary aromaticity
  is allowed in the model.
- [`pi_electron`](@ref)`(mol)`: number of pi electrons
- [`hybridization`](@ref)`(mol)`: orbital hybridization, e.g. `:sp`, `:sp2` and `:sp3`

```@example properties
paraquat = smilestomol("C[n+]1ccc(cc1)c2cc[n+](cc2)C.[Cl-].[Cl-]")
html_fixed_size(paraquat, 250, 250, atomindex=true)
```

Each of the methods returns a vector. Here they are printed in one line each for brevity.

```@example properties
for f in (atom_symbol, atom_charge, multiplicity, lone_pair, implicit_hydrogens,
          valence, is_aromatic, pi_electron, hybridization)
    println(rpad(f, 20), f(paraquat))
end
```

Each atom object can also be accessed by indexing (e.g. `paraquat[2]`), and the functions
above (e.g. `atom_symbol(paraquat[2])`) are also available for the atom object. Note that
the atom object holds only the primary properties given by the input (e.g. SMILES), so the
molecule-level functions (e.g. `atom_charge(paraquat)`) should be used to obtain the values
after the preprocessing and standardization described below.

```@example properties
paraquat[2]
```

The following methods return a vector of bond properties in lexicographic order of edges
(called edge rank, the same order as `edges(mol)` iterates).

- [`bond_order`](@ref)`(mol)`: only integer bond orders are allowed in the model
- [`is_edge_aromatic`](@ref)`(mol)`: whether the bond is aromatic or not
- [`is_rotatable`](@ref)`(mol)`: whether the bond is rotatable or not

```@example properties
fluconazole = smilestomol("OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F")
html_fixed_size(fluconazole, 250, 250, atomindex=true)
```

```@example properties
for f in (bond_order, is_edge_aromatic, is_rotatable)
    println(rpad(f, 18), f(fluconazole))
end
```

[`edge_rank`](@ref MolecularGraph.edge_rank)`(mol)` returns a dictionary of edges and their positions in the bond
property vectors. Edges are normalized so that the source index is smaller than the
destination index ([`u_edge`](@ref) does this for you).

```@example properties
ernk = MolecularGraph.edge_rank(fluconazole)
is_rotatable(fluconazole)[ernk[u_edge(fluconazole, 9, 2)]]
```

### Graph topology

- [`sssr`](@ref)`(mol)`: smallest set of smallest rings
- `degree(mol)`: the degree of each vertex (from Graphs.jl)
- `connected_components(mol)`: connected components (from Graphs.jl)
- [`fused_rings`](@ref)`(mol)`: fused ring systems (2-edge connected components)
- [`ring_count`](@ref), [`which_ring`](@ref), [`smallest_ring`](@ref) and
  [`is_in_ring`](@ref) return the ring membership of each atom.

```@example properties
stannsoporfin = smilestomol(raw"[Cl-].[Cl-].[Sn+4].CCC1=C2[N-]C(C=C3N=C(C=C4[N-]C(=CC5=NC(=C2)C(C)=C5CC)C(C)=C4CCC(O)=O)C(CCC(O)=O)=C3C)=C1C")
html_fixed_size(stannsoporfin, 300, 300, atomindex=true)
```

```@example properties
sssr(stannsoporfin)
```

```@example properties
println(degree(stannsoporfin))
```

```@example properties
connected_components(stannsoporfin)
```

```@example properties
fused_rings(stannsoporfin)
```

```@example properties
println(ring_count(stannsoporfin))
```


## Auto-update mechanism of properties

The following are examples of properties often considered and implemented in molecular
graph models.

- atom (vertex) properties
  - atom symbol, e.g. C, O, N, ...
  - atomic charge
  - whether aromatic or not
- bond (edge) properties
  - bond order
  - whether aromatic or not
- graph topology
  - connected components (a molecule object with multiple molecules)
  - smallest set of smallest rings (SSSR)
- graph properties
  - metadata, e.g. name, compound ID

What is important is that these properties interact with each other. In other words, many of
the vertex/edge properties depend on other properties of the vertex/edge or of adjacent
vertices/edges. This also means that removal or addition of vertices/edges, or changes in the
graph topology, can affect individual properties.

To reduce this kind of side effect, keep consistency and obtain reproducible results,
molecule properties in this library are implemented based on the following categories.

- primary property, or simply 'property': information on atom symbols, charges, bond
  orders, etc. that defines the molecule obtained from SMILES, SDFile or databases.
- secondary property, or 'descriptor': properties that depend on primary properties and the
  graph topology. Methods to generate secondary properties are required to be 'pure
  functions'. That means these methods should take only the graph object (`SimpleGraph`) and
  primary property vectors as arguments and should not alter objects outside the method.

(Note that these terms are not common technical terms. They are categorized as above just for
implementation convenience.)

Primary properties are stored in atom/bond property objects (e.g. [`SMILESAtom`](@ref)).
On the other hand, secondary properties are calculated ad hoc, or retrieved from the cache
(stored in the `descriptors` field of the graph-level property `mol.gprops`) if they are
considered expensive and frequently used.

```@example properties
MolecularGraph.get_descriptor(mol, :sssr)
```

If primary properties or the graph topology are changed (e.g. vertices are removed), the
molecule object is marked as having updates (`mol.state.has_updates`), and the descriptors are
recalculated by a pre-defined routine (the `on_update` callback; see [Preprocessing](preprocessing.md)
for details) the next time a descriptor is accessed.

```@example properties
html_fixed_size(mol, 250, 250, atomindex=true)
```

```@example properties
println(is_aromatic(mol))
```

Removing a vertex of the benzene ring makes the remaining atoms non-aromatic.

```@example properties
let mol = copy(mol)
    rem_vertex!(mol, 6)
    println("has updates: ", mol.state.has_updates)
    arom = is_aromatic(mol)  # descriptors are recalculated here
    println("has updates: ", mol.state.has_updates)
    println(arom)
end
```

See [MolecularGraph interfaces](@ref) for details of the update mechanism.
