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
        include("./graph/udgraph.jl")
        include("./graph/isomorphism.jl")
        include("./graph/translate.jl")
    end

    using LinearAlgebra
    using Printf
    using StaticArrays
    using Statistics
    using YAML
    using ..Error
    using ..Geometry
    using ..Graph

    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/molgraph.jl")

    include("topology.jl")
    include("annotation.jl")
    include("substructure.jl")
    include("remover.jl")

    include("./draw/base.jl")
    include("./draw/coords2d.jl")
    include("./draw/draw2d.jl")
    include("./draw/svg.jl")

    include("download.jl")
    include("sdfilereader.jl")
    include("smilesreader.jl")

end
