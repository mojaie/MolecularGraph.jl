#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const VERSION = begin
    io = open(joinpath(dirname(@__FILE__), "..", "..", "Project.toml"))
    readuntil(io, "version = \"")
    VersionNumber(readuntil(io, "\""))
end


# for test
testfile(file::String) = joinpath(dirname(@__FILE__), "..", "..", "assets", "test", file)