
# Concept of molecular graph models


### Scope of MolecularGraph.jl

MolecularGraph.jl mainly targets on small organic molecules that we deal with in medicinal chemistry. In particular, MolecularGraph.jl deals with a subset of graphs with the following characteristics:

- Undirected graph
- Simple graph, that means no self-loops and multi-edges
- Typically <100 vertices and <100 edges
- Typically degree <=4
- Almost all are planar and many of them are outerplanar
- Vertices, edges and the graph itself have multiple properties (attributes) e.g. atom symbol, bond order and metadata of the molecule
- Properties have interactions. That means changes in graph topology and vertex/edge properties can affect other vertex/edge properties.

Therefore, following molecules are not supported and may be better to use other appropriate models.

- Inorganic molecules
- Bio/synthetic polymers


### Considerations in molecular graph implementation

According to the characteristics shown above,

- Some of the graph algorithm functions in this library do not support some of general graph structure. For example, self-loops and multi-edges are not considered.
- Graph algorithm implementations for general graphs are not always optimal for molecular graphs as well. Even well-known implementations such as deposited to Graph.jl may be reimplemented depending on benchmark results.
- The properties of nodes and edges in molecular graphs are very different from the attributes like distances and labels in general graphs. At this time, MolecularGraph.jl does not depend on MetaGraphs.jl or MetaGraphsNext.jl and has its own data structure to deal with propagation and interaction of molecular graph properties.


### Auto-update mechanism of properties

The following are examples of properties often considered and implemented in molecular graph models

- atom (vertex) properties
  - atom symbol e.g. C, O, N, ...
  - atomic charge
  - whether aromatic or not
- bond (edge) properties
  - bond order
  - whether aromatic or not
- graph topology
  - connected components (a molecule object with multiple molecules)
  - smallest set of smallest rings (SSSR)
- graph properties
  - metadata e.g. name, compound ID

What is important is that these properties can have interactions. In other words, many of the vertex/edge properties depend on other properties in the vertex/edge or adjacent vertexes/edges. This also means that removal or addition of vertices/edges or changes in graph topology can affect individual properties.

To reduce this kind of side-effects, keep consistency and obtain reproducible results, molecule properties in this library is implemented based on following categories.

- primary property, or simply 'property': Information on atom symbol, charges, bond orders, etc. that define the molecule obtained from SMILES, SDFile or databases.
- secondary property, or 'descriptor': Properties that depend on primary properties and graph topology. Methods to generate secondary properties are required to be 'pure function'. That means these methods should take only graph object (SimpleGraph) and primary property vector as arguments and should not alter objects outside the method.

(Note that these terms are not common technical terms. Just for implementation convenience, they are categorized as above.)

Primary properties are stored in atom/bond property objects (e.g. SMILESAtom). On the other hand, secondary properties are calculated ad-hoc, or called from caches (stored in MolGraph.state field) if it is considered to be expensive and frequently called.

If primary properties or graph topology were changed (e.g. remove vertices), the molecule object is marked as `:has_updates`, and recalculated with a pre-defined routine (see preprocessing tutorial for details) when the next time the secondary property is called.