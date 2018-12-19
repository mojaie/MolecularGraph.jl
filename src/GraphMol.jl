#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GraphMol

    export
        Error,
        Geometry,
        Graph

    module Error
        include("exception.jl")
    end

    module Geometry
        using LinearAlgebra
        using StaticArrays
        include("geometry.jl")
    end

    module Graph
        using StaticArrays
        using ..Error
        include("./graph/interface.jl")
        include("./graph/ugraph.jl")
        include("./graph/ugraphview.jl")
        include("./graph/dgraph.jl")
        include("./graph/dgraphview.jl")

        include("./graph/merge.jl")
        include("./graph/linegraph.jl")
        include("./graph/dag.jl")

        include("./graph/shortestpath.jl")
        include("./graph/bipartite.jl")
        include("./graph/triangle.jl")
        include("./graph/bridge.jl")
        include("./graph/component.jl")
        include("./graph/cycle.jl")
        include("./graph/vf2.jl")
        include("./graph/vf2edge.jl")
    end

    using LinearAlgebra
    using Printf
    using StaticArrays
    using Statistics
    using YAML
    using ..Error
    using ..Geometry
    using ..Graph

    include("./model/interface.jl")
    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/molgraph.jl")

    include("./annotation/base.jl")
    include("./annotation/topology.jl")
    include("./annotation/elemental.jl")
    include("./annotation/rotatable.jl")
    include("./annotation/aromatic.jl")
    include("./annotation/funcgroup.jl")

    include("properties.jl")
    include("substructure.jl")
    include("remover.jl")

    include("./smarts/base.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/molecule.jl")

    include("./draw/base.jl")
    include("./draw/coords2d.jl")
    include("./draw/draw2d.jl")
    include("./draw/svg.jl")

    include("download.jl")
    include("sdfilereader.jl")

end
