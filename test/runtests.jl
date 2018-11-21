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


module GraphTest
    using Test
    using GraphMol.Graph
    using GraphMol.Error
    include("./graph/udgraph.jl")
    include("./graph/isomorphism.jl")
    include("./graph/translate.jl")
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
    include("substructure.jl")
    include("remover.jl")

    include("./smarts/base.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/molecule.jl")
    include("./smarts/smiles.jl")
    include("./smarts/smarts.jl")
    include("./smarts/smilesreader.jl")

    include("./draw/base.jl")
    # include("./draw/coords2d.jl")
    # include("./draw/svg.jl")

    include("sdfilereader.jl")

end
