#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolGraph, GeneralMol, GeneralMolView,
    Atom, QueryAtom,
    Bond, QueryBond


abstract type GraphMol <: GraphView end
abstract type GeneralMol <: GraphMol end
abstract type GeneralMolView <: GeneralMol end

# Union types
# TODO: use traits
# https://github.com/JuliaLang/julia/issues/2345


abstract type Atom <: AbstractNode end
abstract type QueryAtom <: AbstractNode end

abstract type Bond <: UndirectedEdge end
abstract type QueryBond <: UndirectedEdge end
