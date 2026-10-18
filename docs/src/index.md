
# MolecularGraph.jl

**MolecularGraph.jl** is a graph-based molecule modeling and chemoinformatics analysis toolkit fully implemented in Julia.

[README.md on GitHub](https://github.com/mojaie/MolecularGraph.jl)


## Installation

```
pkg> add MolecularGraph
```


## Usage

- [Getting Started](@ref): read molecules, draw structures, calculate properties and search substructures
- Basics
  - [Concept of molecular graph models](@ref): scope, graph operations, atom and bond properties
  - [Core functions](@ref)
  - [Properties and descriptors](@ref): built-in properties and the auto-update mechanism
  - [Drawing molecules](@ref): 2D structure images (SVG/PNG), 2D coordinates and 3D rendering
  - [Preprocessing](basics/preprocessing.md): hydrogens, components, charges, resonance structures and custom updaters
  - [Stereochemistry](basics/stereo.md)
  - [Mass and isotopes](@ref): molecular weight, exact mass, isotopic composition and mass spectra
  - [Substructure and query](@ref): substructure match, InChI, SMARTS queries and structural alerts
  - [Maximum common substructure](@ref): MCIS, MCES and topological constraints (tdMCS)
- Advanced topics
  - [MolecularGraph interfaces](@ref)
  - [Customize atom and bond properties](@ref)

Interactive [Pluto.jl notebook tutorials](https://github.com/mojaie/MolecularGraph.jl_notebook)
are also available (some of them are based on older versions).


## License

[MIT license](http://opensource.org/licenses/MIT)  
See [Assets/README.md](https://github.com/mojaie/MolecularGraph.jl/tree/master/assets) for details of external data sets and their licenses.


## Copyright

(C) 2018-2025 Seiji Matsuoka and contributors