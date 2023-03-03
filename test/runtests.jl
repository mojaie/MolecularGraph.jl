#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraphTest
    using Test
    using LinearAlgebra
    using Logging
    using Graphs
    using MolecularGraph

    include("./util/iterator.jl")
    include("./util/math.jl")

    include("./geometry/cartesian.jl")
    include("./geometry/internal.jl")

    include("./graph/cycle.jl")
    include("./graph/matching.jl")

    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/molgraph.jl")
    include("./model/query.jl")

    include("sdfilereader.jl")
    include("sdfilewriter.jl")

    include("./smarts/base.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/molecule.jl")
    include("./smarts/smiles.jl")
    include("./smarts/smarts.jl")

    include("stereo.jl")
    include("preprocess.jl")
    include("coords.jl")

    include("properties.jl")
    include("mass.jl")
    include("wclogp.jl")
    include("inchi.jl")

    include("./draw/base.jl")
    include("./draw/svg.jl")
    include("./draw/3d.jl")

    """
    include("./graph/dfs.jl")
    include("./graph/shortestpath.jl")
    include("./graph/generator.jl")

    include("./graph/plaingraph.jl")
    include("./graph/plaindigraph.jl")
    include("./graph/plainhypergraph.jl")

    include("./graph/multigraph.jl")
    include("./graph/dag.jl")
    include("./graph/connectivity.jl")
    include("./graph/subgraphview.jl")
    include("./graph/disjointunion.jl")
    include("./graph/linegraph.jl")

    include("./graph/triangle.jl")
    include("./graph/clique.jl")
    include("./graph/bipartite.jl")
    include("./graph/planarity.jl")
    include("./graph/product.jl")
    include("./graph/isomorphism/vf2.jl")
    include("./graph/isomorphism/cliquebased.jl")

    include("structurematch.jl"

    """
end
