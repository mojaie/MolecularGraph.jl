```@meta
CurrentModule = MolecularGraph
```


# Core functions


- [`vproptype`](@ref)
- [`eproptype`](@ref)


## SimpleGraph support function

- [`edge_rank`](@ref)
- [`edge_neighbors`](@ref)
- `ordered_neighbors`
- `ordered_edge_neighbors`

- [`u_edge`](@ref)
- [`add_u_edge!`](@ref)
- [`rem_u_edge!`](@ref)

- [`get_descriptor`](@ref)
- [`has_descriptor`](@ref)
- [`set_descriptor!`](@ref)


- [`atom_number`](@ref)
- [`atom_symbol`](@ref)
- [`atom_charge`](@ref)
- [`multiplicity`](@ref)
- [`isotope`](@ref)

- [`bond_order`](@ref)


## Serialization

- [`JSON.json`](@ref MolecularGraph.JSON.json(::MolecularGraph.JSONMolGraph)): serialize a molecule into JSON text
- [`mol_from_json`](@ref): read a molecule from JSON text
