
using Documenter, MolecularGraph


makedocs(
    sitename="MolecularGraph.jl",
    pages = [
        "Home" => "index.md",
        "API" => [
            "Molecular properties" => "api/properties.md",
            "Graph" => "api/graph.md"
        ],
        "Further readings" => [
            "Python interface" => "python.md",
            "Design of molecular graph models" => "moleculargraph.md"
        ]
    ]
)

deploydocs(repo="github.com/mojaie/MolecularGraph.jl.git")
