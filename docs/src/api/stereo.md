
# Stereochemistry

```@autodocs
Modules = [MolecularGraph]
Pages   = ["stereo.jl"]
Private = false
Order   = [:constant, :function, :type]
```

## Stereochemistry and graph modification

Non-exported functions called inside `rem_edge!`, `rem_vertex!` and `rem_vertices!`.

```@docs
MolecularGraph.stereo_on_edge_removal!
MolecularGraph.stereo_on_vertex_removal!
```

## 2D notation of stereocenters

Non-exported functions to keep 2D wedges consistent with the stereocenters.

```@docs
MolecularGraph.assign_stereocenter_wedges!
MolecularGraph.stereocenter_from_draw2d
```

## Stereochemistry types

```@docs
Stereocenter
StereocenterMap
Stereobond
StereobondMap
```
