#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    VERSION

const VERSION = begin
    project = joinpath(dirname(@__FILE__), "..", "..", "Project.toml")
    v = []
    for line in eachline(project)
        if startswith(line, "version")
            push!(v, strip(strip(split(line, " = ")[2]), ['"']))
            break
        end
    end
    v[1]
end
