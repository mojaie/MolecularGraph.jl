#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolGraph, GeneralMol,
    Atom, QueryAtom,
    Bond, QueryBond


abstract type MolGraph <: GraphView end
abstract type GeneralMol <: MolGraph end

# Union types
# TODO: use traits
# https://github.com/JuliaLang/julia/issues/2345


abstract type Atom <: AbstractNode end
abstract type QueryAtom <: AbstractNode end

abstract type Bond <: UndirectedEdge end
abstract type QueryBond <: UndirectedEdge end
