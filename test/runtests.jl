#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GeometryTest

    using Test
    using GraphMol.Geometry

    include("geometry.jl")

end


module GraphModelTest

    using Test
    using GraphMol.GraphModel

    include("./model/undirectedgraph.jl")

end


module GraphMolTest

    using Test
    import GraphMol

    include("./model/atom.jl")
    include("./model/moleculargraph.jl")
    include("topology.jl")
    include("sdfilereader.jl")

end
