#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraph

    module Util
        include("./util/meta.jl")
        include("./util/iterator.jl")
        include("./util/math.jl")
    end

    module Graph
        using ..Util

        include("./graph/interface.jl")

        include("./graph/plaingraph.jl")
        include("./graph/plaindigraph.jl")
        include("./graph/plainhypergraph.jl")

        include("./graph/generator.jl")
        include("./graph/multigraph.jl")
        include("./graph/dag.jl")
        include("./graph/inducedsubgraph.jl")
        include("./graph/product.jl")
        include("./graph/disjointunion.jl")
        include("./graph/linegraph.jl")

        include("./graph/dfs.jl")
        include("./graph/shortestpath.jl")
        include("./graph/bipartite.jl")
        include("./graph/triangle.jl")
        include("./graph/clique.jl")
        include("./graph/connectivity.jl")
        include("./graph/cycle.jl")
        include("./graph/planarity.jl")

        include("./graph/isomorphism/edgeinduced.jl")
        include("./graph/isomorphism/vf2.jl")
        include("./graph/isomorphism/cliquemcs.jl")
    end

    module Geometry
        using LinearAlgebra
        using ..Util

        include("./geometry/interface.jl")

        include("./geometry/cartesian.jl")
        include("./geometry/internal.jl")
    end

    using LinearAlgebra
    using Statistics
    using YAML
    using DelimitedFiles
    using .Util
    using .Graph
    using .Geometry

    export Util, Graph, Geometry

    include("./model/interface.jl")
    include("./draw/interface.jl")

    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/molgraph.jl")

    include("substructure.jl")
    include("mcs.jl")
    include("stereo.jl")
    include("preprocess.jl")
    include("properties.jl")
    include("mass.jl")
    include("wclogp.jl")
    include("funcgroup.jl")
    include("inchi.jl")

    include("./smarts/base.jl")
    include("./smarts/atom.jl")
    include("./smarts/bond.jl")
    include("./smarts/logicaloperator.jl")
    include("./smarts/molecule.jl")

    include("coordgen.jl")
    include("./draw/draw2d.jl")
    include("./draw/svg.jl")

    include("download.jl")
    include("sdfilereader.jl")
    include("sdfilewriter.jl")

    using Requires

    function __init__()
        @require AbstractPlotting="537997a7-5e4e-5d89-9595-2241ea00577e" include("draw/draw3d.jl")
    end
end
