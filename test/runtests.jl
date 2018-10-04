
using Test
using GraphMol.GraphModel
using GraphMol.MolecularModel
import GraphMol

const tests = [
    "sdfilereader",
    "./model/atom", "./model/moleculargraph", "./model/undirectedgraph"
]

for test in tests
    include("$(test).jl")
end
