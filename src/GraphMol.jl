#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GraphMol

    export
        Geometry,
        GraphModel,
        loadsdfiter,
        loadsdfmol

    module Geometry
        include("geometry.jl")
    end

    module GraphModel
        include("./model/undirectedgraph.jl")
    end

    using GraphMol.GraphModel

    include("./model/atom.jl")
    include("./model/bond.jl")
    include("./model/moleculargraph.jl")
    include("sdfilereader.jl")

end
