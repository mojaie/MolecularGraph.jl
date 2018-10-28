#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GeometryTest
    using Test
    using LinearAlgebra
    using StaticArrays
    using GraphMol.Geometry
    using GraphMol.GraphMolError
    include("geometry.jl")
end


module GraphModelTest
    using Test
    using GraphMol.GraphModel
    using GraphMol.GraphMolError
    include("./model/udgraph.jl")
end


module MolecularModelTest
    using Test
    using GraphMol.GraphMolError
    using GraphMol.MolecularModel
    # include("./model/atom.jl")
    include("./model/molgraph.jl")
end


module GraphMolIOTest
    using Test
    using GraphMol.GraphMolError
    using GraphMol.GraphMolIO
    using GraphMol.GraphMolIO: parsesdfatom, parsesdfbond, parsesdfmol
    using GraphMol.GraphMolIO: tokenize, parsesmiles!, parsesmilestoken,
                               parsesmilesatom, parsesmilesbond
    using GraphMol.MolecularModel
    include("sdfilereader.jl")
    include("smilesreader.jl")
end


module BaseAnnotationTest
    using Test
    using GraphMol.BaseAnnotation
    using GraphMol.BaseAnnotation: resolve_inclusion, canonicalize_cycle
    using GraphMol.GraphMolError
    using GraphMol.GraphMolIO
    using GraphMol.MolecularModel
    include("topology.jl")
    include("annotation.jl")
end


module DrawingTest
    using Test
    using GraphMol.GraphMolIO
    using GraphMol.MolecularModel
    using GraphMol.Drawing
    include("./draw/base.jl")
    # include("./draw/coords2d.jl")
    # include("./draw/svg.jl")
end
