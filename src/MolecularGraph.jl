#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraph

    export
        Util,
        Geometry,
        Graph

    module Util
        include("./util/iterator.jl")
    end

    module Geometry
        using Formatting
        using LinearAlgebra
        using StaticArrays

        include("./geometry/interface.jl")

        include("./geometry/coords2d.jl")
        # include("./geometry/coords3d.jl")
        # include("./geometry/coordsinternal.jl")
        # include("./geometry/embed2d.jl")
        # include("./geometry/embed3d.jl")
    end

    module Graph
        using StaticArrays
        using ..Util

        include("./graph/interface.jl")

        include("./graph/ugraph.jl")
        include("./graph/dgraph.jl")
        include("./graph/graphutil.jl")

        include("./graph/view/base.jl")
        include("./graph/view/inducedsubgraph.jl")
        include("./graph/view/complementgraph.jl")
        include("./graph/view/reversegraph.jl")

        include("./graph/generator.jl")
        include("./graph/merge.jl")
        include("./graph/linegraph.jl")
        include("./graph/dag.jl")
        include("./graph/modularproduct.jl")

        include("./graph/shortestpath.jl")
        include("./graph/bipartite.jl")
        include("./graph/triangle.jl")
        include("./graph/clique.jl")
        include("./graph/bridge.jl")
        include("./graph/component.jl")
        include("./graph/cycle.jl")

        include("./graph/isomorphism/base.jl")
        include("./graph/isomorphism/cliquebased.jl")
        include("./graph/isomorphism/vf2.jl")
    end

    using LinearAlgebra
    using Printf
    using Formatting
    using StaticArrays
    using Statistics
    using YAML
    using ..Util
    using ..Geometry
    using ..Graph

    include("./annotation/interface.jl")
    include("./model/interface.jl")
    include("./draw/interface.jl")

    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/molgraph.jl")
    include("./model/molgraphview.jl")
    include("./model/molgraphutil.jl")

    include("./annotation/base.jl")
    include("./annotation/topology.jl")
    include("./annotation/elemental.jl")
    include("./annotation/rotatable.jl")
    include("./annotation/aromatic.jl")
    include("./annotation/wclogp.jl")
    include("./annotation/funcgroup.jl")

    include("properties.jl")
    include("substructure.jl")
    include("mcs.jl")
    include("preprocess.jl")

    include("./smarts/base.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/molecule.jl")

    include("./draw/color.jl")
    include("./draw/base.jl")
    include("./draw/draw2d.jl")
    include("./draw/svg.jl")

    include("download.jl")
    include("sdfilereader.jl")

end
