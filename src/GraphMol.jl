#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GraphMol

    export
        Error,
        Geometry,
        GraphModel

    module Error
        include("exception.jl")
    end

    module Geometry
        using LinearAlgebra
        using StaticArrays
        include("geometry.jl")
    end

    module GraphModel
        using StaticArrays
        using ..Error
        include("./model/udgraph.jl")
    end

    using LinearAlgebra
    using Printf
    using StaticArrays
    using Statistics
    using YAML
    using ..Error
    using ..GraphModel
    using ..Geometry

    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/molgraph.jl")

    include("topology.jl")
    include("annotation.jl")

    include("./draw/base.jl")
    include("./draw/coords2d.jl")
    include("./draw/draw2d.jl")
    include("./draw/svg.jl")

    include("download.jl")
    include("sdfilereader.jl")
    include("smilesreader.jl")

end
