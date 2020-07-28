
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
            "Substructure match" => "moleculargraph/substructure.md",
            "MCS" => "moleculargraph/mcs.md",
            "Functional group detection" => "moleculargraph/funcgroup.md"
        ],
        "MolecularGraph.Graph" => [
            "Interface" => "graph/interface.md",
            "Graph models" => "graph/model.md",
            "Generator" => "graph/generator.md",
            "Traversal" => "graph/traversal.md",
            "Topology" => "graph/topology.md",
            "Operations" => "graph/operation.md",
            "Clique" => "graph/clique.md",
            "Isomorphism" => "graph/isomorphism.md"
        ],
        "MolecularGraph.Geometry" => [
            "Cartesian" => "geometry/cartesian.md"
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
