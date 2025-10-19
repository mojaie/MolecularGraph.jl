
# MolecularGraph.jl

**MolecularGraph.jl** is a graph-based molecule modeling and chemoinformatics analysis toolkit fully implemented in Julia.

[README.md on GitHub](https://github.com/mojaie/MolecularGraph.jl)


## Installation

```
(@v1.8) pkg> add MolecularGraph
```


## Usage

See [Pluto.jl notebook tutorials](https://github.com/mojaie/MolecularGraph.jl_notebook)

- [Getting started](https://mojaie.github.io/MolecularGraph.jl_notebook/getting_started.jl.html)
- [Molecular graph basics](https://mojaie.github.io/MolecularGraph.jl_notebook/molecular_graph_basics.jl.html)
  - Scope of MolecularGraph.jl
  - Considerations in molecular graph implementation
  - Basic operations provided by Graphs.jl interface
  - MolGraph type and atom/bond properties
- [Properties and descriptors](https://mojaie.github.io/MolecularGraph.jl_notebook/properties_and_descriptors.jl.html)
  - Built-in molecule properties and descriptors
    - Lipinski's Rule of five (RO5)
    - Molecular formula
    - Atom and bond properties
    - Graph topology (ring and fused ring)
  - Auto-update mechanism of properties
- [Preprocessing](https://mojaie.github.io/MolecularGraph.jl_notebook/preprocessing.jl.html)
  - Remove hydrogen vertices
  - Extract molecules of interest
  - Standardize charges
  - Dealing with resonance structure
  - Customize property updater
- [Mass and isotopes](https://mojaie.github.io/MolecularGraph.jl_notebook/mass_and_isotopes.jl.html)
  - Molecular weight and exact mass
  - Uncertainty
  - Isotopic composition
  - Simulate mass spectrum
- [Substructure and query](https://mojaie.github.io/MolecularGraph.jl_notebook/substructure_and_query.jl.html)
  - Substructure match
  - InChI and InChIKey
  - SMARTS query
  - Structural alerts (e.g. PAINS)
  - Functional group analysis
  - Query containment
- [Maximum common substructure (MCS)](https://mojaie.github.io/MolecularGraph.jl_notebook/maximum_common_substructure.jl.html)
  - Maximum common induced substructure (MCIS)
  - Maximum common edge-induced substructure (MCES)
  - Connected or disconnected MCS
  - Working with larger molecules
  - Topological constraint (tdMCS)
- [Drawing molecule](https://mojaie.github.io/MolecularGraph.jl_notebook/drawing_molecule.jl.html)
  - Settings of 2D structure images
    - Change image size
    - Layout for web and Pluto notebook
  - Regenerate 2D coordinates
  - 3D molecule rendering using Makie.jl


## License

[MIT license](http://opensource.org/licenses/MIT)  
See [Assets/README.md](https://github.com/mojaie/MolecularGraph.jl/tree/master/assets) for details of external data sets and their licenses.


## Copyright

(C) 2018-2025 Seiji Matsuoka and contributors