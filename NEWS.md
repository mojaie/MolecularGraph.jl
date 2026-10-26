# NEWS

## Unreleased

- Atoms can have a specified valence (`valence` field of `SDFAtom` and `SMILESAtom`,
  `atom_valence`), e.g. carbenes and pentavalent nitrogens. The specified valence is used by
  `valence`, `implicit_hydrogens` and `lone_pair` (`-1`: inferred from the element as before).
  It is read from and written to molfiles (V2000 `vvv` column, V3000 `VAL`) and molecule JSON.
- SMILES bracket atoms have no implicit hydrogens, following the SMILES specification
  (e.g. `[CH2]` is a carbene and `[CH3]` is a methyl radical; before, `[CH2]` and `[CH3]` were
  CH4). `[Si]`, `[Se]` etc. without hydrogen count are now atoms with valence 0.
  Bracket atoms of elements other than typical organic atoms (e.g. `[Sn]`, `[Na+]`) always
  have the specified valence (written as `vvv` to molfiles, e.g. 15 for `[Sn]`), and so do
  CDXML atoms of these elements with the number of hydrogens.
  The valence of aromatic bracket atoms is also determined from the SMILES (the same as
  `apparent_valence` before kekulization, see below), e.g. `c1cc[c]cc1` is a phenyl radical
  (C6H5; before, C6H6).
- `apparent_valence` of aromatic atoms (`isaromatic` of the atom type) before kekulization
  counts aromatic bonds as single bonds plus one pi bond (e.g. 3 for benzene `c`; before, 2),
  so that molecules that are not kekulized have the correct implicit hydrogens (the valence
  of atoms without pi bonds such as `[nH]` and `o` is overestimated until kekulization).
  `kekulize!` overwrites `apparent_valence` with the kekulized bond orders, and
  `smiles_on_update!` and `cdxml_on_update!` no longer call `apparent_valence!` after
  `kekulize!` (custom callbacks should not call it after `kekulize!`).
- `kekulize` assigns double bonds to the aromatic atoms without double bonds whose total bond
  order is lower than the typical valence (or the specified valence, `atom_valence`), so
  `[se]1cccc1`, `[cH-]1cccc1`, `[cH+]1cccccc1` and `[o+]1ccccc1` can be kekulized (failed
  before).
- `is_ring_aromatic`: a ring atom with an exocyclic double bond to oxygen keeps the pi
  electrons other than the double bond, so pyridine N-oxide in the hypervalent form
  (`O=n1ccccc1`) is aromatic.
- Molfile writer: `vvv` is written for atoms with the specified valence (`atom_valence`) and
  atoms whose `valence` descriptor is not inferred from the element and bonds (e.g. set by
  `set_descriptor!`). Molfiles of non-organic atoms are written back as read (`vvv` 0 or 15).
- `lone_pair` is calculated from the valence electrons, charge, valence and multiplicity
  (`lone_pair(symbol_arr, charge_arr, multiplicity_arr, valence_arr)`); the former
  `lone_pair(symbol_arr, charge_arr, connectivity_arr)` method is removed. Lone pairs of
  hypervalent atoms with double bonds are corrected (e.g. 0 for the nitrogen of `CN(=O)=O`,
  1 before), which also changes their `hybridization` (`:sp2` instead of `:sp3`).
- CDXML: the number of hydrogens (`NumHydrogens`, and `hhh`/`HCOUNT` of molfiles read as
  `CDXMLMolGraph`) is converted to the specified valence of `CDXMLAtom` (`valence` field) after
  the bonds are parsed. The `num_hs` field/keyword of `CDXMLAtom` (added after v1.0.0) is
  replaced by `valence`, and `printv2mol` writes the valence (`vvv`) instead of the hydrogen
  count (`hhh`) for CDXML molecules. Aromatic atoms (bond order 1.5) ignore the number of
  hydrogens.
- Fixed: V3000 molfiles could not be read as `CDXMLMolGraph` (no keyword constructor of
  `CDXMLBond`), and V2000 property blocks (`M  CHG` etc.) failed for `CDXMLAtom`.

## v1.0.0

The API is considered stable from this version. The documentation now covers all exported
functions and types (API reference) and the former Pluto notebook tutorials are merged into
the documentation (the notebook repository will be archived).

- Breaking changes
  - Callbacks are now part of the molecule type for type stability: `MolState{T,F1,F2}` has
    constant `on_init`/`on_update` fields, and `MolGraph` and `QueryMolGraph` have an additional
    type parameter (`MolGraph{T,V,E,S}`). Callbacks of an existing molecule cannot be replaced
    in place (e.g. `mol.state.on_update = f`); use `MolecularGraph.with_callbacks(mol; on_update=f)`
    instead. `MolGraph{T,V,E}` is no longer a concrete type (constructors work as before).
  - Stereochemistry is updated when the molecular graph is modified.
    - `rem_edge!` of a bond at a stereocenter or a stereobond (double bond) removes the
      stereochemistry with a warning. `rem_vertex!`/`rem_vertices!` of a heavy atom adjacent to
      a stereocenter or a stereobond also removes it with a warning.
    - Hydrogens referred by stereocenters/stereobonds can be removed safely; the references are
      replaced with the remaining neighbors.
    - 2D wedges are kept consistent with `mol[:stereocenter]` (the reference of stereochemistry)
      after graph modification and coordinate generation. This affects SDFile output, InChI and
      drawing (e.g. wedges on removed hydrogens are regenerated on other bonds).
  - Three-coordinate stereocenters (with an implicit hydrogen) in 2D molfiles are interpreted
    geometrically, consistently with InChI. The result differs from previous versions when the
    neighbors are in a half-plane with the wedge on the middle bond, or there are two wedges in
    the same direction.
  - `printv2mol`/`sdfilewriter` with `givebackhydrogen=true` keep the existing 2D coordinates
    (the hydrogens are placed with wedges instead of calling coordgen), and do not add
    hydrogens to molecules with 3D coordinates.
  - `bondhighlight` of `drawsvg`/`drawpng` accepts only existing bonds (`ArgumentError`).
  - `isotopic_composition(atomsymbol, number)` returns records sorted by mass in ascending order,
    like `isotopic_composition(mol)`.
  - `set_stereocenter!` and `set_stereobond!` return `nothing`.
- Compatible with `juliac --trim=safe` (JuliaC.jl). An experimental C shared library built by
  JuliaC is available in `build/juliac` (the PackageCompiler version is kept until the trim
  option becomes stable).
- Complicated query containment tasks are now more efficient with a built-in SAT solver.
- `Statistics` is no longer a dependency.
- Molecule JSON is written and read by own code (the format is unchanged, about 2x faster).
  `mol_from_json(f, json)` can be used to construct molecules with a type-stable callback.
- Fixed SMILES/SMARTS starting with a chiral atom (e.g. `[C@@H](F)(Cl)Br`).
- Fixed `coordgen!` not displaying wedges of some stereocenters.
- Fixed shifted 2D wedges after `add_edge!`/`rem_edge!`, and SSSR not recalculated after the
  removal of a ring bond.
- RDKitMinimalLib.jl extension is not available on Windows (RDKit_jll for Windows does not
  export the RDKit C API).

## v0.23.0

- Breaking changes
  - The minimum required Julia version will be updated from 1.9 to 1.10.
  - `maximum_clique` was renamed to `maximum_clique_mg` to avoid name confliction (#146).
  - Now compatible with Makie 0.24 (#147). Removed `MakieCore` and now `Makie` is available as an extention.
- `Aqua` was added as a test dependency (#146).
- CDXML reader (#145)

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
