
using Documenter, MolecularGraph


makedocs(
    sitename="MolecularGraph.jl",
    pages = [
        "Home" => "index.md",
        "MolecularGraph" => [
            "I/O" => "moleculargraph/io.md",
            "Structure drawing" => "moleculargraph/draw.md",
            "Chemical properties" => "moleculargraph/properties.md",
            "Preprocessing" => "moleculargraph/preprocess.md",
            "Substructure match" => "moleculargraph/substructure.md",
            "MCS" => "moleculargraph/mcs.md",
            "Functional group detection" => "moleculargraph/funcgroup.md"
        ],
        "GraphModel" => [
            "Interface" => "graph/interface.md",
            "Generator" => "graph/generator.md",
            "Traversal" => "graph/traversal.md",
            "Topology" => "graph/topology.md",
            "Operations" => "graph/operation.md",
            "Clique" => "graph/clique.md",
            "Isomorphism" => "graph/isomorphism.md"
        ],
        "Python interface" => "python.md",
        "Design of molecular graph models" => "design.md"
    ]
)

deploydocs(repo="github.com/mojaie/MolecularGraph.jl.git")
