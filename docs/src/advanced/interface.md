```@meta
CurrentModule = MolecularGraph
```


# MolecularGraph interfaces

This section describes the type hierarchy of molecular graph models and the API functions
that a customized molecular graph model should implement (or can use).

```
Graphs.AbstractGraph{T}
└─ AbstractMolGraph{T}
   └─ SimpleMolGraph{T}                    topology in Graphs.SimpleGraph{T}
      └─ ReactiveMolGraph{T,V,E}           auto-update of descriptors
         ├─ MolGraph{T,V,E,S}              V <: AbstractAtom, E <: AbstractBond
         └─ QueryMolGraph{T,V,E,S}         V, E <: QueryTree
```

Built-in molecule readers return `MolGraph` (e.g. `SDFMolGraph` =
`MolGraph{Int,SDFAtom,SDFBond}`) or `QueryMolGraph` (SMARTS). A molecular graph model
has the following fields:

- `graph`: graph topology (`Graphs.SimpleGraph{T}` for `SimpleMolGraph`)
- `vprops`: vertex properties (atoms), keyed by vertex
- `eprops`: edge properties (bonds), keyed by edge
- `gprops`: graph-level properties and stored descriptors
- `state`: update flags and callbacks (`ReactiveMolGraph` only)


## Common graph interface

[`AbstractMolGraph`](@ref) implements the following functions of the `Graphs.AbstractGraph`
interface by forwarding them to the `graph` field, so most of the graph algorithms in
Graphs.jl can be applied to molecules directly.

- `Graphs.edges`
- `Graphs.edgetype`
- `Graphs.has_edge`
- `Graphs.has_vertex`
- `Graphs.inneighbors`
- `Graphs.ne`
- `Graphs.nv`
- `Graphs.outneighbors`
- `Graphs.vertices`
- `Graphs.is_directed`: As typical molecular graph models are undirected,
  `is_directed(::Type{<:AbstractMolGraph})` returns `false`.

The following `Base` functions are also implemented (and recommended for custom models).

- `Base.eltype`: returns the graph element type (`Int` by default)
- `Base.copy`: returns a copy of the molecule
- `Base.zero`: returns a null graph
- `Base.:(==)`: checks if two graphs have exactly the same graph topology and properties
- `Base.show`: for pretty printing


## SimpleMolGraph interface

[`SimpleMolGraph`](@ref) is the base abstract type of molecular graph models that have
`Graphs.SimpleGraph` as the topology.

`SimpleGraph` does not have an undirected edge type, so an undirected edge is represented
by `Graphs.Edge` whose source index is smaller than the destination index. Edges given to
the API functions are normalized in this form, and edge properties are stored with the
normalized edges as keys.

- [`u_edge`](@ref): returns the normalized edge (`u_edge(mol, 3, 1) == Edge(1, 3)`)
- [`edge_rank`](@ref): returns the indices of edges in `edges(mol)` order. Edge descriptors
  (e.g. [`bond_order`](@ref)) are vectors in this order.
- [`edge_neighbors`](@ref): returns the neighbors of both ends of the edge
- `ordered_neighbors`, `ordered_edge_neighbors`: aliases of `neighbors` and `edge_neighbors`,
  used where the vertices should be in ascending order (e.g. stereochemistry).
- `Graphs.induced_subgraph`: returns the subgraph molecule with the properties, and the vmap.


### Working with properties

- `Base.getindex`: `mol[v]` returns the atom of the vertex `v`, `mol[u, v]` or `mol[e]` returns
  the bond of the edge, and `mol[:stereocenter]` returns the graph-level property.
- `Base.haskey`: `haskey(mol, :stereocenter)` checks the graph-level property,
  `haskey(mol, "compound_id")` checks the metadata (e.g. SDFile option fields).
- [`vproptype`](@ref), [`eproptype`](@ref): return the atom and bond types.
- [`get_descriptor`](@ref), [`has_descriptor`](@ref): access to the stored descriptors.
  Descriptor functions (e.g. [`sssr`](@ref), [`valence`](@ref)) return the stored values if
  available, otherwise calculate them. Use the `Val` form (e.g. `get_descriptor(mol, Val(:sssr))`)
  inside performance-critical code to keep the return type concrete.


### Editing molecules

The following functions edit the molecule. For models other than `ReactiveMolGraph`, the stored
descriptors are not updated automatically.

- `Base.setindex!`: `mol[v] = atom`, `mol[u, v] = bond` or `mol[e] = bond` replaces the atom or bond.
- `Graphs.add_vertex!(mol, atom)`: adds a vertex with the atom.
- `Graphs.add_edge!(mol, u, v, bond)`: adds an edge with the bond. This calls
  [`add_u_edge!`](@ref) with the normalized edge.
- `Graphs.rem_edge!(mol, u, v)`: removes the edge. This calls [`rem_u_edge!`](@ref).
- `Graphs.rem_vertex!`, `Graphs.rem_vertices!`: remove vertices. As with `Graphs.SimpleGraph`,
  the vertex indices are changed by the removal (`rem_vertices!` returns the vmap).
- [`set_descriptor!`](@ref): stores a descriptor (typically called in update callbacks).

A custom `SimpleMolGraph` model should implement [`add_u_edge!`](@ref), [`rem_u_edge!`](@ref),
`Graphs.add_vertex!`, `Graphs.rem_vertex!` and `Graphs.rem_vertices!` to keep the properties
consistent with the topology.


## ReactiveMolGraph interface

[`ReactiveMolGraph`](@ref) is a molecular graph model that automatically updates the stored
descriptors (e.g. standardized charges and bond orders, rings, valence) when the molecule is
edited. The update mechanism works as follows.

1. When a molecule is constructed, [`initialize!`](@ref) calls the `on_init` callback once
   (e.g. stereochemistry from the input), and then calls the `on_update` callback.
1. Editing functions (`setindex!`, `add_vertex!`, `add_edge!`, `rem_vertex!`, ...) do not
   recalculate descriptors but set the update flag by [`notify_updates!`](@ref).
   `add_edge!` also calls [`notify_new_edges!`](@ref) because rings (`sssr`) should be
   recalculated only if new edges are added.
1. Vertex removal remaps the vertex indices in the graph-level properties by
   [`remap_gprops!`](@ref) (see [`remap!`](@ref) below).
1. When a descriptor is accessed by [`get_descriptor`](@ref), [`dispatch_update!`](@ref) calls
   the `on_update` callback if the update flag is set, and then clears the flag. Therefore,
   multiple edits are followed by only one recalculation.

The flags and callbacks are stored in the `state` field ([`MolState`](@ref), a subtype of
[`AbstractState`](@ref)). The callbacks are type parameters of `MolState`, so they are fixed
when the molecule is constructed. Use [`with_callbacks`](@ref) to get a molecule with
other callbacks.

### Callbacks

The default callbacks of the built-in readers are:

| Reader | `on_init` | `on_update` |
|:--|:--|:--|
| [`sdfilereader`](@ref), [`sdftomol`](@ref) | [`sdf_on_init!`](@ref) | [`sdf_on_update!`](@ref) |
| [`smilestomol`](@ref) | [`smiles_on_init!`](@ref) | [`smiles_on_update!`](@ref) |
| [`smartstomol`](@ref) | [`smarts_on_init!`](@ref) | [`default_on_update!`](@ref) |

[`default_on_init!`](@ref) and [`default_on_update!`](@ref) do nothing. Custom callbacks
can be given to the readers to customize the preprocessing and descriptor calculation
(e.g. to skip unnecessary calculations for performance).

```julia
function my_on_update!(mol::MolecularGraph.SimpleMolGraph)
    smiles_on_update!(mol)
    # additional descriptors or preprocessing
end

mol = smilestomol("CC(=O)O"; on_update=my_on_update!)
mol = sdftomol("aspirin.mol"; config=(on_init=sdf_on_init!, on_update=my_on_update!))
mol = mol_from_json(json; on_update=my_on_update!)
mol2 = MolecularGraph.with_callbacks(mol; on_update=smiles_on_update!)
```

Inside callbacks, use `set_descriptor!` to store descriptors. The update is disabled while
the callback is running, so descriptor functions called in the callback do not trigger
another update.


## Graph property interfaces

Graph-level property types implement `Base.:(==)` (all fields are compared), `Base.copy`
and [`remap!`](@ref).

### Graph-level property (gprop)

- [`AbstractProperty`](@ref): the base type of graph-level properties
  - [`SimpleMolProperty`](@ref): the base type of graph-level properties of `SimpleMolGraph`
    - [`MolProperty`](@ref): graph-level properties of `MolGraph` (stereochemistry, metadata,
      logs and descriptors)
    - [`MolDescriptor`](@ref): stored descriptors (the `descriptors` field of `MolProperty`)
    - [`QueryMolProperty`](@ref): graph-level properties of `QueryMolGraph`

When vertices are removed, [`remap!`](@ref) is called for each field of the property as
`remap!(Val(:field), gprop, vmap, edges)` to update the vertex indices (e.g. stereocenters).
Fields that do not depend on vertex indices need no method (the fallback does nothing).
Descriptors that are recalculated by `on_update` also need no remapping, but cached
descriptors that are not recalculated (e.g. `sssr` without new edges, coordinates) are remapped.


### Vertex and edge property (vprop and eprop)

Atom and bond types are subtypes of [`AbstractElement`](@ref), which provides `Base.:(==)` and
`Base.hash` comparing all fields, and `atom[:field]` access to the fields. Atoms and bonds
should be immutable; use `setindex!` on the molecule to replace them.

- [`AbstractElement`](@ref)
  - [`AbstractAtom`](@ref): atom types should implement the following functions
    (descriptor functions of molecules such as [`atom_symbol`](@ref)`(mol)` use them).
    - [`atom_number`](@ref)
    - [`atom_symbol`](@ref)
    - [`atom_charge`](@ref)
    - [`multiplicity`](@ref)
    - [`isotope`](@ref)
    - [`StandardAtom`](@ref): the base type of atoms that have all of the above
      (e.g. [`SDFAtom`](@ref), [`SMILESAtom`](@ref)). Most of the built-in property functions
      are implemented for this type.
  - [`AbstractBond`](@ref): bond types should implement [`bond_order`](@ref).
    - [`StandardBond`](@ref): the base type of standard bonds (e.g. [`SDFBond`](@ref),
      [`SMILESBond`](@ref)).
  - [`QueryTree`](@ref): atom and bond queries of `QueryMolGraph`
    - [`QueryAtom`](@ref)
    - [`QueryBond`](@ref)

See [Customize atom and bond properties](@ref) for custom atom and bond types.

!!! note
    `JSON.json(mol)` and [`mol_from_json`](@ref) support only the built-in molecule types
    (SDFile, SMILES, SMARTS and CommonChem). Molecules with custom atom or bond types are
    serialized by the generic method of JSON.jl, which is not compilable with `juliac --trim`.
