#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
module UtilTest
    using Test
    using MolecularGraph.Util
    include("./util/iterator.jl")
end
"""

module GeometryTest
    using Test
    using LinearAlgebra
    using StaticArrays
    using Formatting
    using MolecularGraph.Geometry

    include("./geometry/coords2d.jl")
    # include("./geometry/embed2d.jl")
end

"""
module GraphTest
    using Test
    using MolecularGraph.Graph

    include("./graph/ugraph.jl")
    include("./graph/dgraph.jl")

    include("./graph/view/inducedsubgraph.jl")
    include("./graph/view/reversegraph.jl")

    include("./graph/generator.jl")
    include("./graph/merge.jl")
    include("./graph/linegraph.jl")
    include("./graph/dag.jl")

    include("./graph/shortestpath.jl")
    include("./graph/triangle.jl")
    include("./graph/clique.jl")
    include("./graph/bipartite.jl")
    include("./graph/bridge.jl")
    include("./graph/component.jl")
    include("./graph/cycle.jl")
    include("./graph/modularproduct.jl")
    include("./graph/isomorphism/vf2.jl")
    include("./graph/isomorphism/cliquebased.jl")
end


module MolecularGraphTest
    using Test
    using StaticArrays
    using MolecularGraph
    using MolecularGraph.Geometry
    using MolecularGraph.Graph

    include("./model/atom.jl")
    include("./model/molgraph.jl")

    include("sdfilereader.jl")
    include("sdfilewriter.jl")

    include("./smarts/base.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/molecule.jl")
    include("./smarts/smiles.jl")
    include("./smarts/smarts.jl")

    include("./annotation/topology.jl")
    include("./annotation/elemental.jl")
    include("./annotation/rotatable.jl")
    include("./annotation/aromatic.jl")
    include("./annotation/wclogp.jl")
    include("./annotation/funcgroup.jl")

    include("preprocess.jl")
    include("substructure.jl")
    include("mcs.jl")

    include("./draw/base.jl")
    include("./draw/svg.jl")
end
"""
