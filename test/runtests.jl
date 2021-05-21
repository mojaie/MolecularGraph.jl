#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


module MolecularGraphUtilTest
    using Test
    using MolecularGraph.Util
    include("./util/iterator.jl")
    include("./util/math.jl")
end


module MolecularGraphGraphTest
    using Test
    using MolecularGraph.Graph

    include("./graph/dfs.jl")
    include("./graph/shortestpath.jl")
    include("./graph/generator.jl")

    include("./graph/plaingraph.jl")
    include("./graph/plaindigraph.jl")
    include("./graph/plainhypergraph.jl")

    include("./graph/multigraph.jl")
    include("./graph/dag.jl")
    include("./graph/cycle.jl")
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
end


module MolecularGraphGeometryTest
    using Test
    using LinearAlgebra
    using MolecularGraph.Geometry

    include("./geometry/cartesian.jl")
    include("./geometry/internal.jl")
end


module MolecularGraphTest
    using Test
    using MolecularGraph
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

    include("properties.jl")
    include("preprocess.jl")
    include("stereo.jl")
    include("mass.jl")
    include("wclogp.jl")
    include("structurematch.jl")
    include("funcgroup.jl")
    include("mcs.jl")
    include("inchi.jl")

    include("./draw/base.jl")
    include("./draw/svg.jl")
    include("./draw/3d.jl")
end
