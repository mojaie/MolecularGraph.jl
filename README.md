
MolecularGraph.jl
===================================================

[![DOI](https://zenodo.org/badge/151080560.svg)](https://zenodo.org/badge/latestdoi/151080560)
[![Build Status](https://travis-ci.org/mojaie/MolecularGraph.jl.svg?branch=master)](https://travis-ci.org/mojaie/MolecularGraph.jl)
![TagBot](https://github.com/mojaie/MolecularGraph.jl/workflows/TagBot/badge.svg)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mojaie.github.io/MolecularGraph.jl/dev)


**MolecularGraph.jl** is a graph-based molecule modeling and chemoinformatics analysis toolkit fully implemented in Julia.


<img src="./assets/image/demo.svg" width="200"/><img src="./assets/image/Bivalirudin.svg" width="200"/><img src="./assets/image/Succinic acid.svg" width="200"/><img src="./assets/image/Docetaxel.svg" width="200"/><img src="./assets/image/FerrocenylethylMaleimide.svg" width="200"/><img src="./assets/image/spacefilling.png" width="200"/>

<img src="./assets/image/mcsdemo.png" width="400"/><img src="./assets/image/funcgroupdemo.png" width="400"/><img src="./assets/image/massspecdemo.png" width="400"/>



## Usage

- [Try it with Jupyter notebook tutorials](https://github.com/mojaie/MolecularGraph.jl_notebook)
- [Documentation and API Reference](https://mojaie.github.io/MolecularGraph.jl/dev)



## Features

- Chemical structure file I/O
  - 2D structure image drawing and export to SVG
  - 3D structure drawing
  - SDFile import/export (.sdf, .mol)
  - SMILES/SMARTS parser

- Database
  - InChI ([InChI](https://www.inchi-trust.org/))
  - Serialization (molecule object <-> JSON)

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
  - Ring, scaffold, connected components
  - Minimum cycle basis (de Pina algorithm)
    - Smallest set of smallest rings (SSSR)
  - Planarity (left-right planarity test)
  - Maximum matching
    - Kekulization
  - Graph traversal

- 2D geometry
  - Stereochemistry drawing
  - Coordinates generation ([coordgenlibs](https://github.com/schrodinger/coordgenlibs))

- Sub(super)structure
  - Library search by using SMARTS query
  - Subgraph isomorphism detection with VF2 algorithm
  - Monomorphism, node-induced and edge-induced
  - Constraints (mandatory/forbidden mapping)

- SMARTS query-based substructure analysis
  - functional group mining
  - structural alerts (by using ChEMBL dataset)

- Maximum common substructure (MCS)
  - By clique detection algorithm
  - Node-induced (MCIS) and edge-induced (MCES)
  - Connected and disconnected
  - Topological constraint (known as tdMCS)
  - Diameter restriction (MCS-DR) and graph-based local similarity (GLS)



## License

[MIT license](http://opensource.org/licenses/MIT)  
See [Assets/README.md](https://github.com/mojaie/MolecularGraph.jl/tree/master/assets) for details of external data sets and their licenses.



## Copyright

(C) 2018-2021 Seiji Matsuoka and contributors
