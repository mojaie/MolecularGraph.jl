
# Design of molecular graph models


## Graph type hierarchy

- (AbstractGraph)
  - (UndirectedGraph)
    - VectorGraph
    - MapGraph(deprecated)
    - (GraphView)
      - (GraphMol)
      - SubgraphView
      - ComplementGraphView
  - (DirectedGraph)
    - DiGraph
    - (DiGraphView)
      - DiSubgraphView
      - ComplementDiGraphView
      - ReverseGraphView
  - MultiGraph
  - MultiDiGraph
  - HyperGraph
  - (GView) (Union{GraphView,DiGraphView})


## Molecule type hierarchy

- (GraphMol)
  - (GeneralMol)
    - VectorMol{A<:Atom, B<:Bond}
      - SDFile (Alias of VectorMol{SDFileAtom, SDFileBond})
      - SMILES (Alias of VectorMol{SmilesAtom, SmilesBond})
    - (GeneralMolView)
      - SubstructureView{T<:VectorGraph}
  - MapMol(deprecated)
  - QueryMol
    - QueryMol{A<:QueryAtom, B<:QueryBond}
      - SMARTS (Alias of QueryMol{SmartsAtom, SmartsBond})



#### Methods

- getatom
- getbond
- neighbors
- neighborcount (or degree)
- atomcount
- bondcount
- updateatom!
- updatebond!
- unlinkatom!
- unlinkbond!


### VectorMol

`VectorMol` is vector(array)-based molecular model which is specialized for element-wise fast computation of molecular properties.

`VectorMol` can iterate over atom properties faster than `MapMol`and can store calculated molecular properties and annotation arrays which are suitable for vector computation. On the other hand, `VectorMol` does not have abilities to modify its graph structure (adding or removing elements). `VectorMol` can be converted to `MapMol` but the calculated properties and annotations will be lost.


### MapMol

`MapMol` is used as a molecular model builder for general purpose.

This type inherits `AbstractMapMol`, a molecular graph model which have map(dict)-based structure. The map-based molecular graph can insert and delete elements (atoms and bonds). This can be easily converted to `VectorMol` object by using `vectormol` method


### QueryMol

`QueryMol` consists of `QueryAtom`s and `QueryBond`s that represents molecular query (ex. atom symbol is 'O' and charge is -1, bond order is 1 and not in rings, ...). This type of objects typically built from SMARTS query.



## Atom

- AbstractNode
  - AbstractAtom
    - Atom
      - SDFileAtom
      - SmilesAtom
    - QueryAtom
      - SmartsAtom



## Bond

- AbstractEdge
  - AbstractBond
    - Bond
      - SDFileBond
      - SmilesBond
    - QueryBond
      - SmartsBond
