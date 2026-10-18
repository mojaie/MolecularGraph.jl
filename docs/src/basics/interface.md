
# MolecularGraph interfaces

This section describes API functions that should be implemented in a customized molecular graph model type.


## Common graph interface

Molecular graphs that is a subtype of [`AbstractMolGraph`](@ref) implements the following standard API functions of `Graphs.AbstractGraph`.

- `Graphs.edges`
- `Graphs.edgetype`
- `Graphs.has_edge`
- `Graphs.has_vertex`
- `Graphs.inneighbors`
- `Graphs.ne`
- `Graphs.nv`
- `Graphs.outneighbors`
- `Graphs.vertices`
- `Graphs.is_directed`: As typical molecular graph models are undirected, `is_directed(::Type{<:AbstractMolGraph})` returns `false`.

Following `Base` functions are also recommended to be implemented

- `Base.eltype`: returns graph element type (`Int` by default)
- `Base.copy`: returns a shallow copy
- `Base.zero`: returns a null graph
- `Base.:(==)`: check if two graphs has the exactly same graph topology and properties.
- `Base.show`: for pretty printing


### SimpleMolGraph interface

If we uses `Graphs.Edge` (`Graphs.SimpleEdge`),


[`SimpleMolGraph`](@ref) is the base abstract type of molecular graph models that have `Graphs.SimpleGraph`.

- [`u_edge`](@ref)
- [`edge_rank`](@ref)
- [`edge_neighbors`](@ref)
- `ordered_neighbors`
- `ordered_edge_neighbors`

- `Graphs.induced_subgraph`


### Working with properties

- `Base.getindex`
- [`get_descriptor`](@ref)
- [`has_descriptor`](@ref)
- [`vproptype`](@ref)
- [`eproptype`](@ref)

editing molecules (recomended to use `ReactiveMolGraph`)

- `Base.setindex!`
- [`set_descriptor!`](@ref)
- [`add_u_edge!`](@ref): called by `Graphs.add_edge!` to store undirected edge
- [`rem_u_edge!`](@ref): called by `Graphs.rem_edge!` to remove undirected edge
- `Graphs.add_vertex!`
- `Graphs.rem_vertex!`
- `Graphs.rem_vertices!`


`ReactiveMolGraph` interface

- [`AbstractState`](@ref)
  - [`dispatch_update!`](@ref)
  - [`notify_updates!`](@ref)
  - [`notify_new_edges!`](@ref)
  - [`MolState`](@ref)

initialize!
remap_gprops
default_on_init!
default_on_update!


## Graph property interfaces

- `Base.:(==)`
- `Base.hash`
- `Base.copy`

### Graph-level property (gprop)

- [`AbstractMolProperty`](@ref)
  - [`to_dict`](@ref)
  - [`SimpleMolProperty`](@ref)
    - [`reconstruct`](@ref)
    - [`remap!`](@ref): for property update (ReactiveMolGraph)
    - [`MolProperty`](@ref)
    - [`MolDescriptor`](@ref)
    - [`QueryMolProperty`](@ref)


### Vertex and edge property (vprop and eprop)

- [`AbstractElement`](@ref)
  - [`to_dict`](@ref)
  - [`AbstractAtom`](@ref)
    - [`atom_number`](@ref)
    - [`atom_symbol`](@ref)
    - [`atom_charge`](@ref)
    - [`multiplicity`](@ref)
    - [`isotope`](@ref)
    - [`StandardAtom`](@ref)
  - [`AbstractBond`](@ref)
    - [`bond_order`](@ref)
    - [`StandardBond`](@ref)
  - [`QueryTree`](@ref)
    - [`QueryAtom`](@ref)
    - [`QueryBond`](@ref)