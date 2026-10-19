using Documenter, MolecularGraph

DocMeta.setdocmeta!(
    MolecularGraph, :DocTestSetup, :(using MolecularGraph); recursive=true)

makedocs(
    sitename="MolecularGraph.jl",
    modules=[MolecularGraph],
    checkdocs=:exports,
    pages = [
        "Home" => "index.md",
        "Getting started" => "gettingstarted.md",
        "Basics" => [
            "Concepts" => "basics/concepts.md",
            "Core functions" => "basics/core.md",
            "Properties and descriptors" => "basics/properties.md",
            "Drawing molecules" => "basics/drawing.md",
            "Preprocessing" => "basics/preprocessing.md",
            "Stereochemistry" => "basics/stereo.md",
            "Mass and isotopes" => "basics/mass.md",
            "Substructure and query" => "basics/substructure.md",
            "Maximum common substructure" => "basics/mcs.md"
        ],
        "Advanced topics" => [
            "Interfaces" => "advanced/interface.md",
            "Editing molecules" => "advanced/editmol.md",
            "Customize atom properties" => "advanced/propinterface.md"
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
            "RDKit interface" => "api/rdkit.md",
            "Graph algorithims" => "api/graph.md",
            "Preprocessing" => "api/preprocess.md",
            "Molecular mass/weight" => "api/mass.md",
            "Stereochemistry" => "api/stereo.md"
        ],
        "Miscellaneous" => "misc.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        # tutorial pages embed structure images as inline SVG
        size_threshold = 500 * 2^10,
        size_threshold_warn = 250 * 2^10,
        example_size_threshold = 64 * 2^10
    )
)

deploydocs(
    repo="github.com/mojaie/MolecularGraph.jl.git",
    devbranch="master",
    push_preview=true
)
