
module GraphMol
    export GraphModel, MolecularModel

    module GraphModel

        export
            UndirectedGraph

        include("undirectedgraph.jl")

    end

    module MolecularModel

        export
            Atom,
            add_atom,
            MolecularGraph

        include("atom.jl")
        include("moleculargraph.jl")

    end




end
