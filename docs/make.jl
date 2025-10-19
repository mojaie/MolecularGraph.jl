
using Pkg
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using Documenter, MolecularGraph


makedocs(
    sitename="MolecularGraph.jl",
    pages = [
        "Home" => "index.md",
        "Getting started" => "gettingstarted.md"
        "Basics" => [
            "Concepts" => "basics/concepts.md",
            "Properties and descriptors" => "basics/properties.md",
            "Drawing" => "basics/drawing.md"
        ],
        "Advanced topics" => [
            "Interfaces" => "basics/interfaces.md",
            "Editing molecules" => "basics/editmol.md",
            "Customize atom properties" => "basics/propinterface.md"
        ],
        "API References" => [
            "Molecular graph models" => "api/model.md",
            "I/O" => "api/io.md",
            "Structure drawing" => "api/draw.md",
            "Coordinates" => "api/coordinates.md",
            "Properties and descriptors" => "api/properties.md",
            "Structure match" => "api/structurematch.md",
            "Molecular queries" => "api/query.md",
            "InChI" => "api/inchi.md",
            "Graph algorithims" => "api/graph.md",
            "Preprocessing" => "api/preprocess.md",
            "Molecular mass/weight" => "api/mass.md",
            "Stereochemistry" => "api/stereo.md"
        ],
        "Miscellaneous" => "misc.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo="github.com/mojaie/MolecularGraph.jl.git",
    push_preview=true
)
