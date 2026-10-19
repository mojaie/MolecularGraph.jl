# RDKit interface

These functions require [RDKitMinimalLib.jl](https://github.com/eloyfelix/RDKitMinimalLib.jl) to be loaded.

!!! note
    RDKitMinimalLib.jl does not work on Windows currently: the RDKit library of RDKit_jll
    for Windows does not export the functions of the RDKit C API (e.g. `could not load symbol
    "get_mol"`).

```@autodocs
Modules = [MolecularGraph]
Pages   = ["rdkit.jl"]
Private = false
Order   = [:constant, :function, :type]
```
