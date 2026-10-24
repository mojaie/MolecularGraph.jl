#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# Package version read from Project.toml at precompile time.
# Project.toml is registered as a dependency so that the precompile cache is invalidated
# when the version is changed.
const VERSION = let
    projectfile = joinpath(dirname(@__FILE__), "..", "..", "Project.toml")
    include_dependency(projectfile)
    open(projectfile) do io
        readuntil(io, "version = \"")
        VersionNumber(readuntil(io, "\""))
    end
end


# for test
testfile(file::String) = joinpath(dirname(@__FILE__), "..", "..", "assets", "test", file)