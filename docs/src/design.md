
# Design of molecular graph models


## Graph type hierarchy

- (AbstractGraph)
  - (UndirectedGraph)
    - (OrderedGraph)
      - PlainGraph
      - ImmutablePlainGraph
      - GraphMol{Atom,Bond}
        - SDFile (Alias of GraphMol{SDFileAtom,SDFileBond})
        - SMILES (Alias of GraphMol{SmilesAtom,SmilesBond})
      - QueryMol{QueryAtom,QueryBond}
        - SMARTS (Alias of QueryMol{SmartsAtom,SmartsBond})
      - LineGraph
      - CartesianProductGraph
      - ModularProductGraph
    - SubgraphView{UndirectedGraph}
  - (DirectedGraph)
    - (OrderedDiGraph)
      - PlainDiGraph
      - FunctionalGroupClassGraph
    - DiSubgraphView{DirectedGraph}
  - HyperGraph?


### AbstractGraph methods

  - getnode, getedge, hasedge
  - `neighbors` and its derivatives
  - nodecount
  - edgecount
  - nodeset
  - edgeset


### DirectedGraph methods

  - `outneighbors` and `inneighbors`


### OrderedGraph methods

  - nodesiter
  - edgesiter
  - nodeattrs
  - edgeattrs


### OrderedGraph

`OrderedGraph` consists of vectors of neighborhood map (incident edge => adjacent node) and edge (tuple of node index pair) vector.



### QueryMol

`QueryMol` consists of `QueryAtom` and `QueryBond` that represent molecular query (ex. atom symbol is 'O' and charge is -1, bond order is 1 and not in rings, ...). This type of objects typically built from SMARTS query.



## Node type hierarchy

- (AbstractNode)
  - (Atom)
    - SDFileAtom
    - SmilesAtom
  - (QueryAtom)
    - SmartsAtom



## Edge type hierarchy

- (AbstractEdge)
  - (UndirectedEdge)
    - Edge
    - (Bond)
      - SDFileBond
      - SmilesBond
    - (QueryBond)
      - SmartsBond
  - (DirectedEdge)
    - Arrow
