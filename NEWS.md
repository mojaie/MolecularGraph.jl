# NEWS

## v0.17.3

- Fixed `mincycles` (#119)
- Fixed package compiler scripts

## v0.17.2

- Fixed `rem_vertex` (#114)
- The default destination path of the compiled package is now `/usr/local/moleculargraphjl`
- The behavior of stereocenters with an unnecessarily high number of wedges has been clarified.
- Fixed serialization of molecules with invalid stereochemistry.

## v0.17.1

- The library builder (built with PackageCompiler) is now Julia v1.11 compatible
- Fixed inappropriate size of atom indices and highlights when exporting very large PNG structure images.

## v0.17

- valence and aromaticity
  - `pi_delocalized` and `hybridization_delocalized` are deprecated. `pi_electron` and `hybridization` always consider contribution of N, O, S lone pairs to the adjacent conjugation system.
  - `valence` works with non-organic atoms and hypervalent atoms, and `hybridization` may better explain molecular geometry (e.g. -SO2-).
  - Aromaticity determination algorithm has been improved, allowing more reasonable detection in many cases even for fused rings (e.g. azulene).
  - (sub)structure match functions matches `atom_symbol` and `hybridization` of each atoms by default. `hybridization` may better explain molecular geometry.
- metadata behavior (e.g. option block in SDFile)
  - `metadata` function would be deprecated (too general and easy to conflict with other packages).
  - `get_prop(mol::MolGraph, key::String)` returns the value of the metadata field.
  - `set_prop!` should work similarly.
  - Still `get_prop(mol::MolGraph, key::Symbol)` returns other graph properties (e.g. stereochem) but in many cases these should be automatically generated and should not be modified.
  - Molecules built from SMILES should be able to have manually curated metadata fields.
  - convenient metadata setter/getters (e.g. mol["compound_id"] = "ABC00001")

## v0.16

Toward version 0.16.0, I have been working on 2D structure images export to PNG (with Cairo) and building binary packages for C and Python.

- `html_fixed_size` and `html_grid` now takes MolGraph objects, not SVG.
- Some molecule parameter methods for 2D drawing (e. g. `is_atom_visible`, `double_bond_style`) are no longer exposed.
- Added Cairo.jl to the dependencies

## v0.15

This release contains significant performance improvements in maximum common substructure (MCS) methods.
Here are some API changes:

- `tcmcis`/`tcmces` are renamed to `tdmcis`/`tdmces` (Abbreviation of Topologically-constrained Disconnected MCS in the original literature)
- libinchi was updated based on InChI version 1.06. Now `inchi` and `inchikey` can take InChI options as optional arguments.
- `pi_electron` and `hybridization` now do not consider delocalization of lone pairs on N, O and S. Default subgraph isomorphism matcher functions
  use `pi_electron` to consider bond orders, so the matching behavior would change a bit.
- `pi_delocalized` and `hybridization_delocalized`, which consider the delocalization, are also available (maybe slightly expensive).
- MCS calculation methods (`disconnected_mcis`/`mces`, `connected_mcis`/`mces` and `tdmcis`/`mces`), and clique detection methods called from them
  now return `Tuple{Dict{T,T},Symbol}`. The dict is the mapping of matched vertices/edges, and the symbol is the status (:done, :targetreached or :timedout)
