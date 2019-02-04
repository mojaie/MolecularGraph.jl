
using Documenter, MolecularGraph


makedocs(
    sitename="MolecularGraph.jl",
    pages = [
        "Home" => "index.md",
        "MolecularGraph" => [
            "Molecule I/O" => "moleculargraph/io.md",
            "Basic chemical properties" => "moleculargraph/properties.md",
            "Preprocessing" => "moleculargraph/preprocess.md",
            "Structure match" => "moleculargraph/structure.md"
        ],
        "Graph" => [
            "Clique" => "graph/clique.md"
        ],
        "Python interface" => "python.md",
        "Design of molecular graph models" => "design.md"
    ]
)

deploydocs(repo="github.com/mojaie/MolecularGraph.jl.git")
