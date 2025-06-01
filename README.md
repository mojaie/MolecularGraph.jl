
MolecularGraph.jl
===================================================

[![DOI](https://zenodo.org/badge/151080560.svg)](https://zenodo.org/badge/latestdoi/151080560)
![TagBot](https://github.com/mojaie/MolecularGraph.jl/workflows/TagBot/badge.svg)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mojaie.github.io/MolecularGraph.jl/dev)


**MolecularGraph.jl** is a graph-based molecule modeling and chemoinformatics analysis toolkit fully implemented in Julia.


<img src="./assets/image/demo.svg" width="200"/><img src="./assets/image/Bivalirudin.svg" width="200"/><img src="./assets/image/Succinic acid.svg" width="200"/><img src="./assets/image/Docetaxel.svg" width="200"/><img src="./assets/image/FerrocenylethylMaleimide.svg" width="200"/><img src="./assets/image/spacefilling.png" width="200"/>

<img src="./assets/image/mcsdemo.png" width="400"/><img src="./assets/image/funcgroupdemo.png" width="400"/><img src="./assets/image/massspecdemo.png" width="400"/>

https://github.com/user-attachments/assets/9844a6f7-d172-4278-82b8-d4f307dff13f


Documentation and API Reference
------------------------------------

[Documentation and API Reference](https://mojaie.github.io/MolecularGraph.jl/dev)


How to use
------------------------------------

[Pluto.jl notebook tutorials](https://github.com/mojaie/MolecularGraph.jl_notebook)

To run codes in your environment, see `Edit or run this notebook` instruction shown in the top-right of the tutorial pages below.

- [Getting started](https://mojaie.github.io/MolecularGraph.jl_notebook/getting_started.html)
- [Molecular graph basics](https://mojaie.github.io/MolecularGraph.jl_notebook/molecular_graph_basics.html)
  - Scope of MolecularGraph.jl
  - Considerations in molecular graph implementation
  - Basic operations provided by Graphs.jl interface
  - MolGraph type and atom/bond properties
- [Properties and descriptors](https://mojaie.github.io/MolecularGraph.jl_notebook/properties_and_descriptors.html)
  - Built-in molecule properties and descriptors
    - Lipinski's Rule of five (RO5)
    - Molecular formula
    - Atom and bond properties
    - Graph topology (ring and fused ring)
  - Auto-update mechanism of properties
- [Preprocessing](https://mojaie.github.io/MolecularGraph.jl_notebook/preprocessing.html)
  - Remove hydrogen vertices
  - Extract molecules of interest
  - Standardize charges
  - Dealing with resonance structure
  - Customize property updater
- [Mass and isotopes](https://mojaie.github.io/MolecularGraph.jl_notebook/mass_and_isotopes.html)
  - Molecular weight and exact mass
  - Uncertainty
  - Isotopic composition
  - Simulate mass spectrum
- [Substructure and query](https://mojaie.github.io/MolecularGraph.jl_notebook/substructure_and_query.html)
  - Substructure match
  - InChI and InChIKey
  - SMARTS query
  - Structural alerts (e.g. PAINS)
  - Functional group analysis
  - Query containment
- [Maximum common substructure (MCS)](https://mojaie.github.io/MolecularGraph.jl_notebook/maximum_common_substructure.html)
  - Maximum common induced substructure (MCIS)
  - Maximum common edge-induced substructure (MCES)
  - Connected or disconnected MCS
  - Working with larger molecules
  - Topological constraint (tdMCS)
- [Drawing molecule](https://mojaie.github.io/MolecularGraph.jl_notebook/drawing_molecule.html)
  - 2D structure images
    - Display and export in SVG format
    - Display and export in PNG format
    - Layout for Pluto notebook
    - Regenerate coordinates
  - 3D molecule rendering using Makie.jl
- [Stereochemistry](https://mojaie.github.io/MolecularGraph.jl_notebook/stereochemistry.html)
  - Stereochemistry as a molecular graph property
  - Stereospecific implicit hydrogens


Features
------------------

- Chemical structure file I/O
  - 2D structure image drawing and export to SVG
  - 3D structure drawing ([Makie.jl](https://github.com/MakieOrg/Makie.jl))
  - SDFile reader/writer (.sdf, .mol)
  - SMILES/SMARTS reader
  - SMILES writer ([RDKitMinimalLib.jl](https://github.com/eloyfelix/RDKitMinimalLib.jl))
  - Coordinates generation ([coordgenlibs](https://github.com/schrodinger/coordgenlibs))

- Properties and descriptors
  - H-bond donor/acceptor
  - rotatable bonds
  - Aromaticity
  - Wildman-Crippen logP
  - Fingerprints ([RDKitMinimalLib.jl](https://github.com/eloyfelix/RDKitMinimalLib.jl))

- Substructure and query
  - InChI ([InChI](https://www.inchi-trust.org/))
  - Serialization (molecule object <-> loss-less format or CommonChem JSON)
  - Subgraph isomorphism detection with VF2 algorithm
    - SMARTS query match
    - Monomorphism, node-induced and edge-induced match
    - Constraints (mandatory/forbidden vertex mapping)
  - Functional group query set
  - Structural alerts detection with ChEMBL dataset
  - SMARTS query containment

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

- Maximum common substructure (MCS)
  - By clique detection algorithm
  - Node-induced (MCIS) and edge-induced (MCES)
  - Connected and disconnected
  - Topological constraint (known as tdMCS)
  - Diameter restriction (MCS-DR) and graph-based local similarity (GLS)


License
-----------------

[MIT license](http://opensource.org/licenses/MIT)  
See [Assets/README.md](https://github.com/mojaie/MolecularGraph.jl/tree/master/assets) for details of external data sets and their licenses.


Copyright
-----------------

(C) 2018-2025 Seiji Matsuoka and contributors
