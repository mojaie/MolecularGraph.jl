
MolecularGraph.jl
===================================================

[![DOI](https://zenodo.org/badge/151080560.svg)](https://zenodo.org/badge/latestdoi/151080560)
[![Build Status](https://travis-ci.org/mojaie/MolecularGraph.jl.svg?branch=master)](https://travis-ci.org/mojaie/MolecularGraph.jl)
![TagBot](https://github.com/mojaie/MolecularGraph.jl/workflows/TagBot/badge.svg)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mojaie.github.io/MolecularGraph.jl/dev)


**MolecularGraph.jl** is a graph-based molecule modeling and chemoinformatics analysis toolkit fully implemented in Julia.


<img src="./assets/image/demo.svg" width="200"/><img src="./assets/image/Acetohexamide.svg" width="200"/><img src="./assets/image/Bivalirudin.svg" width="200"/><img src="./assets/image/Cefmenoxime.svg" width="200"/><img src="./assets/image/Succinic acid.svg" width="200"/><img src="./assets/image/Quinacrine.svg" width="200"/><img src="./assets/image/Docetaxel.svg" width="200"/><img src="./assets/image/FerrocenylethylMaleimide.svg" width="200"/>

<img src="./assets/image/mcsdemo.png" width="400"/><img src="./assets/image/funcgroupdemo.png" width="400"/><img src="./assets/image/massspecdemo.png" width="400"/>



## Usage

- [Try it with Jupyter notebook tutorials](https://github.com/mojaie/MolecularGraph.jl_notebook)
- [Documentation and API Reference](https://mojaie.github.io/MolecularGraph.jl/dev)


## Features

- Chemical structure file I/O
  - Structure image drawing and export to SVG
  - SDFile import/export (.sdf, .mol)
  - SMILES/SMARTS parser

- Basic descriptors
  - H-bond donor/acceptor
  - rotatable bonds
  - Aromaticity
  - Wildman-Crippen logP

- Atomic mass
  - standard atomic/molecular weight
  - relative atomic/molecular mass
  - isotopic composition

- Molecular graph topology
  - Ring, scaffold, connectivity
  - Graph traversal

- 2D geometry
  - Stereochemistry drawing
  - Coordinates generation ([coordgenlibs](https://github.com/schrodinger/coordgenlibs))

- Sub(super)structure
  - Library search by using SMARTS query
  - Subgraph isomorphism detection with VF2 algorithm
  - Node-induced and edge-induced
  - Constraints (mandatory/forbidden mapping)

- SMARTS and terminology graph-based functional group analysis

- Maximum common substructure (MCS)
  - By clique detection algorithm
  - Node-induced (MCIS) and edge-induced (MCES)
  - Connected and disconnected
  - Topological constraint (known as tdMCS)
  - Diameter restriction (MCS-DR) and graph-based local similarity (GLS)


## License

[MIT license](http://opensource.org/licenses/MIT)


## Copyright

(C) 2018-2020 Seiji Matsuoka
