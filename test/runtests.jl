#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GeometryTest
    using Test
    using LinearAlgebra
    using StaticArrays
    using MolecularGraph.Geometry
    using MolecularGraph.Error
    include("geometry.jl")
end


module GraphTest
    using Test
    using MolecularGraph.Graph
    using MolecularGraph.Error
    include("./graph/ugraph.jl")
    include("./graph/dgraph.jl")
    include("./graph/graphview.jl")

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
    include("./graph/vf2.jl")
    include("./graph/vf2edge.jl")
end

module MolecularGraphTest
    using Test
    using StaticArrays
    using MolecularGraph
    using MolecularGraph.Error

    include("./model/atom.jl")
    include("./model/molgraph.jl")

    include("./smarts/base.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/molecule.jl")
    include("./smarts/smiles.jl")
    include("./smarts/smarts.jl")

    include("substructure.jl")

    include("./annotation/topology.jl")
    include("./annotation/elemental.jl")
    include("./annotation/rotatable.jl")
    include("./annotation/aromatic.jl")
    include("./annotation/funcgroup.jl")

    include("remover.jl")

    include("./draw/base.jl")
    # include("./draw/coords2d.jl")
    include("./draw/svg.jl")

    include("sdfilereader.jl")
end
