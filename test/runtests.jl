#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
module GeometryTest
    using Test
    using GraphMol.GraphMolError
    using GraphMol.Geometry
    include("geometry.jl")
end


module GraphModelTest
    using Test
    using GraphMol.GraphMolError
    using GraphMol.GraphModel
    include("./model/undirectedgraph.jl")
end


module MolecularModelTest
    using Test
    using GraphMol.GraphMolError
    using GraphMol.MolecularModel
    include("./model/atom.jl")
    include("./model/moleculargraph.jl")
end


module DescriptorTest
    using Test
    using GraphMol.GraphMolError
    using GraphMol.Descriptor
    using GraphMol.GraphMolIO
    using GraphMol.Descriptor: resolve_inclusion
    include("topology.jl")
end


module GraphMolIOTest
    using Test
    using GraphMol.GraphMolError
    using GraphMol.GraphMolIO
    using GraphMol.GraphMolIO: parseatoms, parsebonds, parsemol
    include("sdfilereader.jl")
end
"""

module DrawingTest
    using Test
    using GraphMol.GraphMolIO
    using GraphMol.Drawing
    include("./draw/svg.jl")
end
