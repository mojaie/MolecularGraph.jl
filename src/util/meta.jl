#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    version


function version()
    project = joinpath(dirname(@__FILE__), "..", "..", "Project.toml")
    for line in eachline(project)
        if startswith(line, "version")
            return strip(strip(split(line, " = ")[2]), ['"'])
        end
    end
end
