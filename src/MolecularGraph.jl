#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraph

    export
        MolecularGraphUtil,
        MolecularGraphGeometry,
        MolecularGraphModel

    module MolecularGraphUtil
        include("./util/meta.jl")
        include("./util/iterator.jl")
    end

    module MolecularGraphModel
        using ..MolecularGraphUtil

        include("./graph/interface.jl")

        include("./graph/ugraph.jl")
        include("./graph/dgraph.jl")
        include("./graph/multigraph.jl")
        include("./graph/graphutil.jl")

        include("./graph/view/base.jl")
        include("./graph/view/inducedsubgraph.jl")
        include("./graph/view/complementgraph.jl")
        include("./graph/view/reversegraph.jl")

        include("./graph/generator.jl")
        include("./graph/merge.jl")
        include("./graph/linegraph.jl")
        include("./graph/dag.jl")
        include("./graph/modularproduct.jl")

        include("./graph/shortestpath.jl")
        include("./graph/bipartite.jl")
        include("./graph/triangle.jl")
        include("./graph/clique.jl")
        include("./graph/connectivity.jl")
        include("./graph/cycle.jl")
        include("./graph/planarity.jl")

        include("./graph/isomorphism/base.jl")
        include("./graph/isomorphism/cliquebased.jl")
        include("./graph/isomorphism/vf2.jl")
    end

    module MolecularGraphGeometry
        using Formatting
        using LinearAlgebra
        using ..MolecularGraphUtil
        using ..MolecularGraphModel

        include("./geometry/interface.jl")

        include("./geometry/coords2d.jl")
        include("./geometry/coords3d.jl")
        include("./geometry/coordsinternal.jl")
    end

    using LinearAlgebra
    using Printf
    using Formatting
    using Statistics
    using YAML
    using ..MolecularGraphUtil
    using ..MolecularGraphGeometry
    using ..MolecularGraphModel

    include("./model/interface.jl")
    include("./draw/interface.jl")

    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/molgraph.jl")
    include("./model/molgraphview.jl")
    include("./model/molgraphutil.jl")

    include("substructure.jl")
    include("mcs.jl")
    include("preprocess.jl")
    include("properties.jl")
    include("wclogp.jl")
    include("funcgroup.jl")

    include("./smarts/base.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/molecule.jl")

    include("./draw/draw2d.jl")
    include("./draw/svg.jl")

    include("download.jl")
    include("sdfilereader.jl")
    include("sdfilewriter.jl")

end
