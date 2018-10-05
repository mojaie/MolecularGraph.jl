#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using Test
using GraphMol
using GraphMol.Geometry
using GraphMol.GraphModel
using GraphMol.MolecularModel

const tests = [
    "sdfilereader",
    "./model/atom", "./model/moleculargraph", "./model/undirectedgraph"
]

for test in tests
    include("$(test).jl")
end
