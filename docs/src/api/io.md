
# Molecule I/O


## SDFile reader

```@autodocs
Modules = [MolecularGraph]
Pages   = ["sdfilereader.jl"]
Private = false
Order   = [:constant, :function, :type]
```


## SDFile writer

```@autodocs
Modules = [MolecularGraph]
Pages   = ["sdfilewriter.jl"]
Private = false
Order   = [:constant, :function, :type]
```


## SMILES/SMARTS reader

```@autodocs
Modules = [MolecularGraph]
Pages   = ["smarts/base.jl"]
Private = false
Order   = [:constant, :function, :type]
```


## CDXML reader

```@autodocs
Modules = [MolecularGraph]
Pages   = ["model/cdxml.jl"]
Private = false
Order   = [:constant, :function, :type]
```


## JSON

Molecules can be serialized into JSON text by `JSON.json(mol)` and restored by
[`mol_from_json`](@ref). `mol_from_json` also reads RDKit CommonChem JSON.

```@docs
MolecularGraph.JSON.json(::MolecularGraph.JSONMolGraph)
mol_from_json
```
