#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraphTest
    using Test
    using LinearAlgebra
    using Logging
    using Graphs
    using MolecularGraph

    """
    include("./util/iterator.jl")
    include("./util/math.jl")

    include("./geometry/cartesian.jl")
    include("./geometry/internal.jl")

    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/molgraph.jl")
    include("./model/query.jl")

    include("./graph/operators.jl")
    include("./graph/cycle.jl")
    include("./graph/bipartite.jl")
    include("./graph/matching.jl")
    include("./graph/planarity.jl")
    include("./graph/clique.jl")
    include("./graph/isomorphism_vf2.jl")
    include("./graph/isomorphism_clique.jl")

    include("sdfilereader.jl")
    include("sdfilewriter.jl")

    include("./smarts/base.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/molecule.jl")
    include("./smarts/smiles.jl")
    include("./smarts/smarts.jl")
    """

    include("stereo.jl")
    include("preprocess.jl")
    include("coords.jl")

    include("properties.jl")
    include("mass.jl")
    include("wclogp.jl")
    include("inchi.jl")

    include("structurematch.jl")
    include("querycontainment.jl")

    include("./draw/base.jl")
    include("./draw/svg.jl")
    include("./draw/3d.jl")

end
