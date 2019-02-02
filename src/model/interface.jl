#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolGraph, VectorMolGraph, MapMolGraph, QueryMolGraph,
    VectorMolView, MapMolView, QueryMolView,
    VectorMol, MapMol, QueryMol,
    MolView, MutableMol,
    AbstractAtom, Atom, QueryAtom,
    AbstractBond, Bond, QueryBond,
    Annotation


abstract type MolGraph <: UndirectedGraphView end
abstract type VectorMolGraph <: MolGraph end
abstract type MapMolGraph <: MolGraph end
abstract type QueryMolGraph <: MolGraph end
abstract type VectorMolView <: MolGraph end
abstract type MapMolView <: MolGraph end
abstract type QueryMolView <: MolGraph end

# Union types
# TODO: use traits
# https://github.com/JuliaLang/julia/issues/2345

VectorMol = Union{VectorMolGraph,VectorMolView}
QueryMol = Union{QueryMolGraph,QueryMolView}
MapMol = Union{MapMolGraph,MapMolView,QueryMol}
MolView = Union{VectorMolView,MapMolView,QueryMolView}
MutableMol = Union{MapMolGraph,QueryMolGraph}


abstract type AbstractAtom <: AbstractNode end
abstract type Atom <: AbstractAtom end
abstract type QueryAtom <: AbstractAtom end

abstract type AbstractBond <: UndirectedEdge end
abstract type Bond <: AbstractBond end
abstract type QueryBond <: AbstractBond end


abstract type Annotation end
