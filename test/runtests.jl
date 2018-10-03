
using Test
using GraphMol.GraphModel
using GraphMol.MolecularModel

const tests = [
    "./model/atom", "./model/moleculargraph", "./model/undirectedgraph"
]

for test in tests
    include("$(test).jl")
end
