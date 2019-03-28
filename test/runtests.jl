#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


module MolecularGraphUtilTest
    using Test
    using MolecularGraph.MolecularGraphUtil
    include("./util/iterator.jl")
end


module MolecularGraphModelTest
    using Test
    using MolecularGraph.MolecularGraphModel

    include("./graph/generator.jl")
    include("./graph/ugraph.jl")
    include("./graph/dgraph.jl")
    include("./graph/multigraph.jl")

    include("./graph/view/inducedsubgraph.jl")
    include("./graph/view/reversegraph.jl")

    include("./graph/merge.jl")
    include("./graph/linegraph.jl")
    include("./graph/dag.jl")

    include("./graph/shortestpath.jl")
    include("./graph/triangle.jl")
    include("./graph/clique.jl")
    include("./graph/bipartite.jl")
    include("./graph/connectivity.jl")
    include("./graph/cycle.jl")
    include("./graph/planarity.jl")
    include("./graph/modularproduct.jl")
    include("./graph/isomorphism/vf2.jl")
    include("./graph/isomorphism/cliquebased.jl")
end


module MolecularGraphGeometryTest
    using Test
    using LinearAlgebra
    using Formatting
    using MolecularGraph.MolecularGraphGeometry

    include("./geometry/coords2d.jl")
    include("./geometry/coords3d.jl")
    include("./geometry/coordsinternal.jl")
end


module MolecularGraphTest
    using Test
    using MolecularGraph
    using MolecularGraph.MolecularGraphGeometry
    using MolecularGraph.MolecularGraphModel

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

    include("preprocess.jl")
    include("substructure.jl")
    include("mcs.jl")
    include("properties.jl")
    include("wclogp.jl")
    include("funcgroup.jl")

    include("./draw/base.jl")
    include("./draw/svg.jl")
end
