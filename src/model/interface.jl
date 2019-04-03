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
