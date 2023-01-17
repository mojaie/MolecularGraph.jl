#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    VERSION, MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION

const VERSION = begin
    io = open(joinpath(dirname(@__FILE__), "..", "..", "Project.toml"))
    readuntil(io, "version = \"")
    readuntil(io, "\"\n")
end

const MAJOR_VERSION = parse(Int, split(VERSION, ".")[1])
const MINOR_VERSION = parse(Int, split(VERSION, ".")[2])
const PATCH_VERSION = parse(Int, split(VERSION, ".")[3])
