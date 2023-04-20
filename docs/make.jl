
using Pkg
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using Documenter, MolecularGraph


makedocs(
    sitename="MolecularGraph.jl",
    pages = [
        "Home" => "index.md",
        "MolecularGraph" => [
            "Molecular graph models" => "moleculargraph/model.md",
            "I/O" => "moleculargraph/io.md",
            "Structure drawing" => "moleculargraph/draw.md",
            "Coordinates" => "moleculargraph/coordinates.md",
            "Properties and descriptors" => "moleculargraph/properties.md",
            "Structure match" => "moleculargraph/structurematch.md",
            "Molecular queries" => "moleculargraph/query.md",
            "InChI" => "moleculargraph/inchi.md",
            "Graph algorithims" => "moleculargraph/graph.md",
            "Preprocessing" => "moleculargraph/preprocess.md",
            "Molecular mass/weight" => "moleculargraph/mass.md",
            "Stereochemistry" => "moleculargraph/stereo.md"
        ],
        "Implementation notes" => [
            "Concept of molecular graph models" => "notes/concept.md",
            "Miscellaneous" => "notes/misc.md"
        ]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo="github.com/mojaie/MolecularGraph.jl.git",
    push_preview=true
)
