
# Stereochemistry

```@autodocs
Modules = [MolecularGraph]
Pages   = ["stereo.jl"]
Private = false
Order   = [:constant, :function, :type]
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
