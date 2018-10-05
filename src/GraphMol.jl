#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GraphMol

    export
        Geometry,
        GraphModel,
        Atom,
        MolecularGraph,
        weight,
        number,
        name,
        color,
        addhydrogen!,
        getatom,
        getbond,
        updateatom!,
        updatebond!,
        required_descriptor,
        loadsdfiter,
        loadsdfmol

    module Geometry
        include("geometry.jl")
    end

    module GraphModel
        include("./model/undirectedgraph.jl")
    end

    include("./model/atom.jl")
    include("./model/moleculargraph.jl")
    include("sdfilereader.jl")

end
