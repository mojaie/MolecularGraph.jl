# NEWS

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
