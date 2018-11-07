#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GeometryTest
    using Test
    using LinearAlgebra
    using StaticArrays
    using GraphMol.Geometry
    using GraphMol.Error
    include("geometry.jl")
end


module GraphModelTest
    using Test
    using GraphMol.GraphModel
    using GraphMol.Error
    include("./model/udgraph.jl")
end


module GraphMolTest
    using Test
    using StaticArrays
    using GraphMol
    using GraphMol.Error
    using GraphMol: resolve_inclusion, canonicalize_cycle

    include("./model/atom.jl")
    include("./model/molgraph.jl")

    include("topology.jl")
    include("annotation.jl")

    include("./draw/base.jl")
    # include("./draw/coords2d.jl")
    # include("./draw/svg.jl")

    include("sdfilereader.jl")
    include("smilesreader.jl")

end
