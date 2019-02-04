
# Design of molecular graph models


## Molecule


- MolGraph
  - MapMolGraph
    - GeneralMapMol{A<:Atom, B<:Bond}
      - SDFile (Alias of GeneralMapMol{SDFileAtom, SDFileBond})
      - SMILES (Alias of GeneralMapMol{SmilesAtom, SmilesBond})
  - QueryMolGraph
    - ConnectedQueryMol{A<:QueryAtom, B<:QueryBond}
      - ConnectedSMARTS (Alias of ConnectedQueryMol{SmartsAtom, SmartsBond})
    - DisconnectedQueryMol{A<:QueryAtom, B<:QueryBond}
      - SMARTS (Alias of DisconnectedQueryMol{SmartsAtom, SmartsBond})
  - VectorMolGraph
    - GeneralVectorMol{A<:Atom, B<:Bond}
  - MapMolView
  - QueryMolView
  - VectorMolView



### AbstractMol

`MolGraph` is the base type of all molecular objects. Fundamental graph element accessor functions and counting functions are available in all molecular model objects.


#### Methods

- getatom
- getbond
- neighbors
- neighborcount (or degree)
- atomcount
- bondcount



### MapMol

`MapMol` is used as a molecular model builder for general purpose.

This type inherits `AbstractMapMol`, a molecular graph model which have map(dict)-based structure. The map-based molecular graph can insert and delete elements (atoms and bonds). This can be easily converted to `VectorMol` object by using `vectormol` method


#### Methods

- updateatom!
- updatebond!
- unlinkatom!
- unlinkbond!
- vectormol


### QueryMol

`QueryMol` consists of `QueryAtom`s and `QueryBond`s that represents molecular query (ex. atom symbol is 'O' and charge is -1, bond order is 1 and not in rings, ...). This type of objects typically built from SMARTS query.

#### Methods

  - updateatom!
  - updatebond!
  - unlinkatom!
  - unlinkbond!



### VectorMol

`VectorMol` is vector(array)-based molecular model which is specialized for element-wise fast computation of molecular properties.

`VectorMol` can iterate over atom properties faster than `MapMol`and can store calculated molecular properties and annotation arrays which are suitable for vector computation. On the other hand, `VectorMol` does not have abilities to modify its graph structure (adding or removing elements). `VectorMol` can be converted to `MapMol` but the calculated properties and annotations will be lost.


 #### Common methods

  - mapmol



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
