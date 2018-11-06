#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GraphMol

    export
        Error,
        Geometry,
        GraphModel,
        MolecularModel,
        Descriptor,
        GraphMolIO

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

    module MolecularModel
        using StaticArrays
        using YAML
        using ..GraphModel
        include("./model/atom.jl")
        include("./model/bond.jl")
        include("./model/molgraph.jl")
    end

    module Drawing
        using LinearAlgebra
        using Printf
        using StaticArrays
        using Statistics
        using ..MolecularModel
        using ..Geometry
        include("./draw/base.jl")
        include("./draw/coords2d.jl")
        include("./draw/draw2d.jl")
        include("./draw/svg.jl")
    end

    module BaseAnnotation
        using ...MolecularModel
        include("topology.jl")
        include("annotation.jl")
    end

    module GraphMolIO
        using StaticArrays
        using ..Error
        using ..MolecularModel
        using ..BaseAnnotation
        include("download.jl")
        include("sdfilereader.jl")
        include("smilesreader.jl")
    end

end
