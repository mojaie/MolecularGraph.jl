
module GraphMol

    export
        GraphModel,
        MolecularModel,
        v2000reader


    module GraphModel

        export
            UndirectedGraph,
            node,
            edge,
            addnode!,
            addedge!

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
            atom,
            bond,
            addatom!,
            addbond!

        include("./model/atom.jl")
        include("./model/moleculargraph.jl")

    end


    include("v2000reader.jl")

end
