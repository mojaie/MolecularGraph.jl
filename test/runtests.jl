
using Test
using GraphMol.GraphModel
using GraphMol.MolecularModel

const tests = [
    "atom", "moleculegraph", "undirectedgraph"
]

for test in tests
    include("$(test).jl")
end
