
# Molecular graph models

## Basic molecular graph type and interfaces

```@autodocs
Modules = [MolecularGraph]
Pages   = ["./model/interface.jl", "./model/molgraph.jl"]
Private = false
Order   = [:constant, :function, :type]
```

## Default atom and bond types

```@autodocs
Modules = [MolecularGraph]
Pages   = ["./model/atom.jl", "./model/bond.jl"]
Private = false
Order   = [:constant, :function, :type]
```

## Interfaces for customized models

Non-exported types and functions to implement or use in customized molecular graph models
(see [MolecularGraph interfaces](@ref)).

```@docs
MolecularGraph.AbstractMolGraph
MolecularGraph.SimpleMolGraph
MolecularGraph.edge_rank
MolecularGraph.get_descriptor
MolecularGraph.has_descriptor
MolecularGraph.set_descriptor!
MolecularGraph.add_u_edge!
MolecularGraph.rem_u_edge!
MolecularGraph.AbstractState
MolecularGraph.MolState
MolecularGraph.dispatch_update!
MolecularGraph.notify_updates!
MolecularGraph.notify_new_edges!
MolecularGraph.AbstractProperty
MolecularGraph.SimpleMolProperty
MolecularGraph.remap!
MolecularGraph.MolProperty
MolecularGraph.MolDescriptor
MolecularGraph.QueryMolProperty
MolecularGraph.AbstractElement
MolecularGraph.StandardAtom
MolecularGraph.StandardBond
```
