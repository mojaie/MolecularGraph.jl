# Getting Started

This tutorial walks through the fundamental operations of molecular data mining with MolecularGraph.jl:

- Chemical structure drawing
- Removing hydrogens
- Calculating molecular properties
- Substructure match and molecule query
- Molecular fingerprints for similarity search
- Maximum common substructure (MCS)


## Installation

MolecularGraph.jl is registered in the General registry. Enter the Pkg mode of the Julia REPL (press `]`) and type:

```
pkg> add MolecularGraph
```

As you will often use fundamental features related to the graph model, it is also recommended to add the Graphs.jl package (`add Graphs`) so that the Graphs.jl interface functions (e.g. `nv`, `ne` and `neighbors`) are available.

[RDKitMinimalLib.jl](https://github.com/eloyfelix/RDKitMinimalLib.jl) can optionally be installed (`add RDKitMinimalLib`) if you need a SMILES writer and similarity fingerprints that are not implemented in MolecularGraph.jl. These features are provided as a package extension, so `using RDKitMinimalLib` enables them automatically.

```@example gettingstarted
using MolecularGraph, Graphs
```


## Chemical structure drawing

### From SDFile

A molecule object can be created from an SDFile (.mol or .sdf) by using [`sdftomol`](@ref). Structure data files for this tutorial were downloaded from [PubChem](https://pubchem.ncbi.nlm.nih.gov/) (e.g. `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6437877/SDF` for Cefditoren Pivoxil, CID 6437877) and are bundled with the documentation sources.

```@example gettingstarted
datadir = joinpath(pkgdir(MolecularGraph), "docs", "src", "assets", "data")
molfile = joinpath(datadir, "Cefditoren Pivoxil.mol")
mol = sdftomol(molfile)
```

In environments that support HTML output (e.g. Jupyter, Pluto.jl and this documentation), a molecule object is displayed as a structure image (250 x 250 px by default). See [Drawing molecules](@ref) for drawing options.


### From SMILES

A molecule object can also be created from a SMILES string by [`smilestomol`](@ref). As SMILES has no atom coordinates, 2D coordinates are calculated automatically. This library internally uses Schrodinger's [coordgenlibs](https://github.com/schrodinger/coordgenlibs) for 2D coordinate generation.

```@example gettingstarted
smilestomol(raw"O=C3N2/C(=C(/C=C\c1scnc1C)CS[C@@H]2[C@@H]3NC(=O)C(=N\OC)/c4nc(sc4)N)C(=O)O")
```


## Removing hydrogens

SDFiles downloaded from PubChem have explicit hydrogen nodes. In practice, hydrogens that carry no important information are removed from molecular graphs for simplicity.

- [`remove_hydrogens!`](@ref) removes "trivial" hydrogen vertices (no charge, no unpaired electron, no specific isotope, and not involved in stereochemistry).
- [`remove_all_hydrogens!`](@ref) removes all hydrogen vertices.

```@example gettingstarted
remove_hydrogens!(mol)
mol
```

Even if these hydrogen nodes are removed from the graph, the actual number of hydrogens can be inferred from the atom valence (so-called "implicit hydrogens"). Most of the molecular property and descriptor calculations in this library take such implicit hydrogens into account. For example, the molecular weight does not change although the number of vertices does.

```@example gettingstarted
original = sdftomol(molfile)
nv(original), standard_weight(original, 2)  # number of vertices and molecular weight
```

```@example gettingstarted
nv(mol), standard_weight(mol, 2)
```

See [Preprocessing](basics/preprocessing.md) for more details on hydrogen removal and other preprocessing operations.


## Calculating molecular properties

As already shown above, `standard_weight(mol, digits)` returns the standard molecular weight. Other fundamental molecular properties can be calculated as follows (see [Properties and descriptors](@ref) for the full list).

```@example gettingstarted
hydrogen_acceptor_count(mol)  # F, O and N
```

```@example gettingstarted
hydrogen_donor_count(mol)  # O or N attached to at least one hydrogen
```

```@example gettingstarted
rotatable_count(mol)  # rotatable bonds
```

```@example gettingstarted
wclogp(mol, 3)  # predicted logP by the Wildman and Crippen method
```

```@example gettingstarted
molecular_formula(mol)  # formula in the Hill system
```


## Substructure match and molecule query

[`has_substruct_match`](@ref)`(mol, query)` returns whether `query` is a substructure of `mol`.

```@example gettingstarted
substr = smilestomol("C1=CCSC2N1C(=O)C2")  # cephem skeleton
has_substruct_match(mol, substr)
```

Similar to SMILES, a SMARTS query can be parsed by [`smartstomol`](@ref).

```@example gettingstarted
# non-aromatic ring nitrogen adjacent to a double bond
query = smartstomol(raw"[$([N;r]A=A)]")
has_substruct_match(mol, query)
```

[`substruct_matches`](@ref) returns an iterator of all matches, each of which is a mapping from the vertices of `mol` to the vertices of the query. The matched atoms and bonds can be highlighted in the structure image.

```@example gettingstarted
substr = smilestomol("s1cncc1")  # thiazole
# matched vertex indices
matched = unique(vcat([collect(keys(m)) for m in substruct_matches(mol, substr)]...))
# matched edges can be obtained from the node-induced subgraph
matched_edges = induced_subgraph_edges(mol.graph, matched)
# highlight matched vertices and edges
html_fixed_size(mol, 300, 300, atomhighlight=matched, bondhighlight=matched_edges)
```

See [Substructure and query](@ref) for more detailed substructure matching and visualization tasks.


## Molecular fingerprints for similarity search

The RDKit extension (requires `using RDKitMinimalLib`) provides some of the fingerprint methods implemented in RDKit. The Morgan fingerprint with radius 2 (almost equivalent to ECFP4) is often used for similarity search. The following code is not executed in this documentation.

```julia
using RDKitMinimalLib

cefditoren = sdftomol(joinpath(datadir, "Cefditoren Pivoxil.mol"))
remove_all_hydrogens!(cefditoren)
ceftazidime = sdftomol(joinpath(datadir, "Ceftazidime.mol"))
remove_all_hydrogens!(ceftazidime)

a = morgan_fp_vector(cefditoren)  # BitVector, default radius=2
b = morgan_fp_vector(ceftazidime)
# Jaccard index (0 = totally dissimilar <-> 1 = perfect match)
sum(a .& b) / sum(a .| b)
```


## Maximum common substructure (MCS)

MolecularGraph.jl provides the following fundamental MCS calculation methods for cheminformatics applications:

- Node-induced (MCIS) or edge-induced (MCES)
- Connected or disconnected
- With topological constraint (known as TD-MCS)

```@example gettingstarted
cefditoren = sdftomol(joinpath(datadir, "Cefditoren Pivoxil.mol"))
remove_all_hydrogens!(cefditoren)
ceftazidime = sdftomol(joinpath(datadir, "Ceftazidime.mol"))
remove_all_hydrogens!(ceftazidime)
html_grid([cefditoren, ceftazidime], 2, 250)
```

[`connected_mces`](@ref) calculates the connected maximum common edge-induced subgraph. The result is a tuple of the edge mapping and the calculation status (`:done` if the calculation was completed, or `:timedout`, etc.).

```@example gettingstarted
mapping, status = connected_mces(cefditoren, ceftazidime)
length(mapping), status
```

```@example gettingstarted
# atoms and bonds highlighting
edges1 = collect(keys(mapping))
nodes1 = unique(vcat([[src(e), dst(e)] for e in edges1]...))
svg1 = html_fixed_size(cefditoren, 250, 250,
    atomhighlight=nodes1, bondhighlight=edges1)

edges2 = collect(values(mapping))
nodes2 = unique(vcat([[src(e), dst(e)] for e in edges2]...))
svg2 = html_fixed_size(ceftazidime, 250, 250,
    atomhighlight=nodes2, bondhighlight=edges2)

HTML("""<div style="display:flex">$(svg1.content)$(svg2.content)</div>""")
```

[`tdmcis`](@ref) calculates the maximum common node-induced subgraph with topological constraint. The `tolerance` keyword argument specifies the allowed mismatch of topological distances between matched atoms.

```@example gettingstarted
mapping, status = tdmcis(cefditoren, ceftazidime, tolerance=1)

# atoms and bonds highlighting
nodes1 = collect(keys(mapping))
edges1 = induced_subgraph_edges(cefditoren.graph, nodes1)
svg1 = html_fixed_size(cefditoren, 250, 250,
    atomhighlight=nodes1, bondhighlight=edges1)

nodes2 = collect(values(mapping))
edges2 = induced_subgraph_edges(ceftazidime.graph, nodes2)
svg2 = html_fixed_size(ceftazidime, 250, 250,
    atomhighlight=nodes2, bondhighlight=edges2)

HTML("""<div style="display:flex">$(svg1.content)$(svg2.content)</div>""")
```

See [Maximum common substructure](@ref) for more details.
