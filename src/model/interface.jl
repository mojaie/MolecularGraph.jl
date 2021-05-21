#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Atom, QueryAtom,
    Bond, QueryBond


abstract type Atom <: AbstractNode end
abstract type QueryAtom <: AbstractNode end

abstract type Bond <: UndirectedEdge end
abstract type QueryBond <: UndirectedEdge end


# Aliases

getatom(mol, i) = nodeattr(mol, i)
getbond(mol, i) = edgeattr(mol, i)
getbond(mol, u, v) = edgeattr(mol, u, v)
hasbond(mol, u, v) = hasedge(mol, u, v)

setatom!(mol, i, attr) = setnodeattr!(mol, i, attr)
setbond!(mol, i, attr) = setedgeattr!(mol, i, attr)
addatom!(mol, attr) = addnode!(mol, attr)
addbond!(mol, u, v, attr) = addedge!(mol, u, v, attr)

atomcount(mol) = nodecount(mol)
bondcount(mol) = edgecount(mol)