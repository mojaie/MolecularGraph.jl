
# Design of molecular graph models


## Graph type hierarchy

- (AbstractGraph)
  - (UndirectedGraph)
    - (OrderedGraph)
      - PlainGraph
      - GraphMol{Atom,Bond}
        - SDFile (Alias of GraphMol{SDFileAtom,SDFileBond})
        - SMILES (Alias of GraphMol{SmilesAtom,SmilesBond})
      - QueryMol{QueryAtom,QueryBond}
        - SMARTS (Alias of QueryMol{SmartsAtom,SmartsBond})
      - LineGraph{OrderedGraph}
      - ModularProductGraph{OrderedGraph}
    - SubgraphView{UndirectedGraph}
  - (DirectedGraph)
    - (OrderedDiGraph)
      - DiGraph
      - FunctionalGroupClassGraph
    - DiSubgraphView{DirectedGraph}
  - HyperGraph


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

  - nodekeys(values)
  - edgekeys(values)
  - nodesiter
  - edgesiter


### OrderedGraph

`GraphMol` is vector(array)-based molecular model which can be used for element-wise fast computation of molecular properties. This is a subtype of `OrderedGraph`



### QueryMol

`QueryMol` consists of `QueryAtom`s and `QueryBond`s that represents molecular query (ex. atom symbol is 'O' and charge is -1, bond order is 1 and not in rings, ...). This type of objects typically built from SMARTS query.



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
