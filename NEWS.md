# NEWS

## v0.22.0

- Breaking changes
  - Moved to JSON.jl v1.
    - `StructUtils` has been added as a direct dependency.
    - Introduced new internal data types for graph-level properties (e.g. `StereocenterMap`, `Coords2d`, etc.).
    - Removed `to_dict` and `to_json`. Use `JSON.json` and `JSON.parse` instead. `mol_from_json(json::AbstractString)` is a convenient interface to auto-detect JSON format.
    - Slightly improved (de)serialization performance.
- Implement generic inchitomol() with stereo support (#138) (by @hhaensel).
- MolGraph default initializers records unusual valence warning to gprops.log.
- Improved error handling in compiled package to prevent segfault.
- Properly scaled coordinates in SDFile output (adapted conventional bond length, 0.825).
- Added tests
- Fixed some descriptor functions and improved SDFile reader performance.
- Fixed some wrong stereochemistry.

## v0.21.1

- Implement `Base.setindex!` molecule property accessor which was actually not implemented in v0.21.0.
- Fixed regression in 2D molecule drawing.
- Fixed a bug in serialization of molecules.

## v0.21.0

- Added `VirtualAtom` as an experimental feature. This enables treating atom placeholders (R-groups) and molecular fragments (e.g., amino acid residues) as atom vertices.
- Enhanced `has_exact_match` and `has_substruct_match` to support stereospecific queries with the option `stereo=true` (#133).
- Updated molecule property accessors (#130). `get_props` and `set_props!` are now deprecated; use `Base.getindex` and `Base.setindex!` instead.
- Fixed `on_init` and `on_update` callbacks for automatic property recalculation to allow adding custom preprocessing methods (e.g., `remove_all_hydrogen!`, `extract_largest_component!`, `protonate_acids!`).
- Fixed SMILES parser to correctly recognize the order of atom properties (#124).
- Improved text placement in 2D drawings.

## v0.20.2

- Reverted some undesirable changes in v0.20
  - `svgcolor` in hex form
  - It is still challenging to run juliac.jl in Docker environment. PackageCompiler.jl will be a primary option for a while longer.
  - Newly introduced `remove_all_hydrogen!` in preprocessing method for substructure match can cause mismatch in isomorphism mappings (#136). This has been removed and alternative workflow will be documented in the tutorial soon.
- Text positioning in 2D drawings has been slightly improved.

## v0.20.1

- The shortest path algorithm used inside `sssr` was specialized for typical molecular graph tasks and get more efficient.
- Implemented `Base.copy` for `MolGraph` and `MolProperty` types (#131). Similar to `copy(graph::SimpleGraph)`, this returns a deepcopy in effect, but far faster than `deepcopy`.
- Fixed RDKit fingerprint options and added some method aliases.
- `MolProperty` types are no longer mutable.
- Fixed sdfilereader (#135)

## v0.20.0

This version introduces substantial breaking changes to internal data structures and methods, while the public APIs remain largely unchanged.

- `RDKitMinimalLib` and `Cairo` will be weak dependencies (Added extensions `RDKitExt` and `CairoExt`). Cairo is a large library but has been used only for PNG image export. RDKitMinimalLib is not so large but there are known OS compatibility issues.
- The graph-level properties of `MolGraph` have been refactored from a `Dict`-based structure to a dedicated data class, `MolProperty`. This may improve type stability and enables a more sophisticated auto-update mechanism for molecular descriptors.
- Molecular query structures (parsed from SMARTS) are no longer recursive. The new query representation is based on `SimpleDiGraph`, allowing for more straightforward topological searches and equality checks.
- SMILES/SMARTS parsing now accepts the special case `[HH]` to represent molecular hydrogen (#124).
- Rewrite many functions to increase type-stablitly (hopefully this improves JIT compile time). 

### API changes

- Atom and bond constructors now accept keyword arguments (e.g., `SMILESAtom(symbol=:N, isaromatic=true)`).
- Options in molecule readers (e.g., SMILES, SDFile) have been changed from a `Dict`-based interface to keyword arguments (e.g., `smilestomol("CCC", on_update=custom_updater!)`).
- Graph-level properties (e.g., SDFile metadata and stereochemistry) are now accessed via fields (e.g., mol.gprops.stereocenter) instead of by index (e.g., `mol.gprops[:stereocenter]`). Metadata accessors remain available (e.g., `mol["compound_id"]` is equivalent to `mol.gprops.metadata["compound_id"]`).
- MolGraph can now store multiple sets of 2D/3D Cartesian coordinates (e.g., generated 3D conformers). As there is currently no API for managing multiple coordinate sets, drawing methods will use the first entry by default.


## v0.19.1

- Fixed SDFilereader (#126)
- Added `make build` command to build ilbrary package instead of the PackageCompiler.jl script. This was enabled by `juliac.jl` script, but still not compatible with `--trim` option.

## v0.19.0

There may be many breaking changes. See updated tutorials for details.

- All imports and exports are aggregated to the package entrypoint file (src/MolecularGraph.jl).
  - Removed unnecessarily exposed APIs, that are internally used or have too general name (e.g. `metadata`) to avoid confliction.
- [RDKitMinimalLib.jl](https://github.com/eloyfelix/RDKitMinimalLib.jl) has been added as a direct dependency.
  - `smiles(mol)` generates SMILES from MolGraph (#67)
  - RDKit fingerprints (Morgan, RDKit, etc.) (#72)
- Geometry features in molecular drawing and stereochemistry are now based on [GeometryBasics.jl](https://github.com/JuliaGeometry/GeometryBasics.jl)
- Color features in molecular drawings are now based on [Colors.jl](https://github.com/JuliaGraphics/Colors.jl)
  - Color options in `drawsvg` and `drawpng` accepts rgb(), hex color codes and any other parsable representations in Colors.jl
- Improved serialization and deserialization
  - Safer JSON deserialization
  - RDKit CommonChem format reader/writer
  - `MolGraph(json::String)` now automatically detect element types and JSON formats.
- Stereochemistry
  - Fixed wrong stereocenter recognition in SMILES
  - Fixed wrong stereobond recognition in <8-membered rings
- Fixed Makie errors in 3D drawing (#112)
  - Axis visibility setting seems to be removed from MakieCore, but now we can set up `LScene` with `show_axis=false` options or apply `hidedecorations!` to `Axis` when we use (See tutorials for example usage)

## v0.18.0

- The minimum required Julia version will be updated from 1.6 to 1.8 (due to OrderedCollections requiring Julia ≥ 1.7.1).
- OrderedCollections has been added as a direct dependency.
- The container for SDFile options (metadata) will be changed from `Dict` to `OrderedDict`, improving the consistency of metadata field order.
- A new function `sdfilescanner` has been added. It works similarly to `sdfilereader` but does not parse the CTAB block. Instead, it returns a dictionary of metadata and the raw CTAB block as a string.

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
