
using Pkg
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using Documenter, MolecularGraph


makedocs(
    sitename="MolecularGraph.jl",
    pages = [
        "Home" => "index.md",
        "MolecularGraph" => [
            "Molecule" => "moleculargraph/molecule.md",
            "I/O" => "moleculargraph/io.md",
            "Structure drawing" => "moleculargraph/draw.md",
            "Chemical properties" => "moleculargraph/properties.md",
            "Preprocessing" => "moleculargraph/preprocess.md",
            "Molecular mass/weight" => "moleculargraph/mass.md",
            "Structure match" => "moleculargraph/structurematch.md",
            "Functional group detection" => "moleculargraph/funcgroup.md"
        ],
        "MolecularGraph.Graph" => "graph.md",
        "MolecularGraph.Geometry" => "geometry.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo="github.com/mojaie/MolecularGraph.jl.git",
    push_preview=true
)
