# Stereochemistry

This page describes the following functionalities related to stereochemistry.

- Stereochemistry as a molecular graph property
- Stereospecific implicit hydrogens
- Limitations

```@example stereo
using MolecularGraph, Graphs
```


## Stereochemistry as a molecular graph property

Stereochemistry properties are often held as primary properties of vertices and edges in
chemical data formats such as SDFile (wedges of bonds) and SMILES (`@`, `@@`, `/` and `\`).

Stereochemistry is an interesting property in the context of graph theory. Despite being a
property about topology, it cannot be represented as a general simple graph. In order to
handle such stereochemistry features, `MolGraph` holds stereochemistry information objects
as graph-level properties, which are obtained from the input when the molecule is
constructed (by the `on_init` callback, see [MolecularGraph interfaces](@ref)).

```@example stereo
tsa = smilestomol(raw"O=C(NO)\C=C\C(=C\[C@H](C(=O)c1ccc(N(C)C)cc1)C)C")
html_fixed_size(tsa, 250, 250, atomindex=true)
```

`mol[:stereocenter]` returns a `StereocenterMap`, a mapping of the stereocenter vertex
to the `Stereocenter` object that has the fields `lookingFrom`, `first`, `second` and
`isclockwise`. Looking from the `lookingFrom` vertex toward the center, `first`, `second`
and the remaining neighbor (or the implicit hydrogen) are arranged clockwise if
`isclockwise` is true, otherwise counterclockwise (the same convention as SMILES `@@`
and `@`).

```@example stereo
for (center, stereo) in tsa[:stereocenter]
    println(center, " => ", Tuple(stereo))  # (lookingFrom, first, second, isclockwise)
end
```

`mol[:stereobond]` returns a `StereobondMap`, a mapping of the double bond (edge) to the
`Stereobond` object that has the fields `first`, `second` and `is_cis`. `first` and
`second` are the neighbors of the source and the destination vertex of the edge,
respectively, and `is_cis` indicates whether they are on the same side of the double bond.

```@example stereo
for (bond, stereo) in tsa[:stereobond]
    println(bond, " => ", Tuple(stereo))  # (first, second, is_cis)
end
```

[`set_stereocenter!`](@ref) and [`set_stereobond!`](@ref) can be used to set the
stereochemistry information manually.


## Stereospecific implicit hydrogens

When a molecule is constructed, stereocenters are rearranged so that the hydrogen vertex at
the stereocenter (if any; see `stereo_hydrogen` and [`safe_stereo_hydrogen!`](@ref))
does not appear in the `Stereocenter` object. Therefore, hydrogen vertices attached to
stereocenters can be safely removed by [`remove_all_hydrogens!`](@ref) (or `rem_vertex!`).
[`remove_hydrogens!`](@ref) keeps such stereospecific hydrogens as vertices.

Here we use an SDFile of estradiol, which can be obtained from PubChem
(e.g. `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/5757/SDF`).

```@example stereo
datadir = joinpath(pkgdir(MolecularGraph), "docs", "src", "assets", "data")
molfile = joinpath(datadir, "Estradiol.mol")
mol1 = sdftomol(molfile)
remove_hydrogens!(mol1)
html_fixed_size(mol1, 250, 250, atomindex=true)
```

```@example stereo
inchikey(mol1)
```

```@example stereo
mol2 = sdftomol(molfile)
remove_all_hydrogens!(mol2)
html_fixed_size(mol2, 250, 250, atomindex=true)
```

The stereochemistry information in `mol[:stereocenter]` is preserved after the removal of
the hydrogens.

```@example stereo
mol1[:stereocenter] == mol2[:stereocenter]
```

`mol[:stereocenter]` is the reference of the stereochemistry, and the wedges in the 2D
drawing are its notation. The wedges on the removed hydrogens are lost with the bonds, but
the wedges are automatically updated to represent `mol[:stereocenter]` (e.g. a wedge on
another bond of the stereocenter; see [`MolecularGraph.assign_stereocenter_wedges!`](@ref)),
while the 2D coordinates taken from the SDFile are kept. Therefore, the drawing above and the
InChI (which is generated from the 2D coordinates and wedges) are consistent with the
stereochemistry.

```@example stereo
inchikey(mol2) == inchikey(mol1)
```

[`coordgen!`](@ref) can be used to regenerate the 2D coordinates and wedges (wedges of
stereocenters not displayed by coordgenlibs are also complemented in the same way).

```@example stereo
coordgen!(mol2)
html_fixed_size(mol2, 250, 250, atomindex=true)
```

```@example stereo
inchikey(mol2) == inchikey(mol1)
```

The 2D rendering of fused rings above is still a bit hard to read, because the
stereospecific hydrogens at the ring fusions are implicit. To give the hydrogens back to the
molecule, simply add hydrogen vertices as follows.

```@example stereo
let mol = copy(mol2)
    ringcount = ring_count(mol)
    for center in keys(mol[:stereocenter])
        if degree(mol, center) == 3 && ringcount[center] > 1
            add_vertex!(mol, SDFAtom(symbol=:H))
            add_edge!(mol, center, nv(mol), SDFBond())
        end
    end
    coordgen!(mol)
    html_fixed_size(mol, 250, 250, atomindex=true)
end
```

The SDFile writer ([`printv2mol`](@ref), [`sdfilewriter`](@ref)) does this automatically
(`givebackhydrogen=true` by default). If the molecule has 2D coordinates, they are kept and
the hydrogens are placed with wedges, so `coordgen!` is not called.

```@example stereo
let dump = printv2mol(mol2)
    mol = sdftomol(IOBuffer(dump))
    inchikey(mol) == inchikey(mol1)
end
```


## Limitations

- As the stereochemistry can be affected by changes of adjacent vertices and edges,
  operations that alter the graph topology of molecules that have stereochemistry should
  be avoided (except for the removal of hydrogen vertices shown above).
- Interconversion of stereochemical data and 2D coordinate generation still have many
  challenges.
