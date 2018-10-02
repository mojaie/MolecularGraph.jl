
using Test

const tests = [
    "atom", "moleculegraph", "undirectedgraph"
]

for test in tests
    include("$(test).jl")
end
