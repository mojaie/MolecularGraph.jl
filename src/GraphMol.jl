
module GraphMol

    export
        GraphModel,
        MolecularModel,
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


    module MolecularModel

        export
            Atom,
            MolecularGraph,
            weight,
            number,
            name,
            color,
            addH!,
            getatom,
            getbond,
            updateatom!,
            updatebond!

        include("./model/atom.jl")
        include("./model/moleculargraph.jl")

    end


    include("sdfilereader.jl")

end
