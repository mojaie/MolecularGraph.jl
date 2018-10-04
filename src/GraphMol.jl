#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module GraphMol

    export
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


    module GraphModel

        export
            UndirectedGraph,
            getnode,
            getedge,
            updatenode!,
            updateedge!

        include("./model/undirectedgraph.jl")

    end

    include("./model/atom.jl")
    include("./model/moleculargraph.jl")
    include("sdfilereader.jl")

end
