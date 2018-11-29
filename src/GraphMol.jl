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

        include("./graph/bridge.jl")
        include("./graph/component.jl")
        include("./graph/cycle.jl")
        include("./graph/isomorphism.jl")
        include("./graph/linegraph.jl")
        include("./graph/merge.jl")
        include("./graph/shortestpath.jl")
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

    include("topology.jl")
    include("annotation.jl")
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
