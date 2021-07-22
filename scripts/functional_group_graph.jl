
const PROJECT_DIR = joinpath(@__DIR__, "..")

using Pkg
Pkg.activate(PROJECT_DIR)

using YAML
using MolecularGraph
using MolecularGraph.Graph

const INPUT_DIR = joinpath(PROJECT_DIR, "./assets/raw/functionalgroup")
const OUTPUT_FILE = joinpath(PROJECT_DIR, "./assets/const/functionalgroup.yaml")


function isaatommatch(q1::QueryMol, q2::QueryMol)
    return function (qa1, qa2)
        return query_contains(nodeattr(q2, qa2).query, nodeattr(q1, qa1).query)
    end
end

function isabondmatch(q1::QueryMol, q2::QueryMol)
    return function (qb1, qb2)
        return query_contains(edgeattr(q2, qb2).query, edgeattr(q1, qb1).query)
    end
end


function hasatommatch(q1::QueryMol, q2::QueryMol)
    return function (qa1, qa2)
        return isequivalent(nodeattr(q1, qa1).query, nodeattr(q2, qa2).query)
    end
end

function hasbondmatch(q1::QueryMol, q2::QueryMol)
    return function (qb1, qb2)
        return isequivalent(edgeattr(q1, qb1).query, edgeattr(q2, qb2).query)
    end
end



function run()
    # read list of functional groups from .yaml files
    fgrecords = []
    for f in readdir(INPUT_DIR)
        append!(fgrecords, YAML.load(open(joinpath(INPUT_DIR, f))))
    end
    for r in fgrecords
        delete!(r, "isa")
        delete!(r, "has")
    end
    isaedges = Set{Tuple{Int,Int}}()
    hasedges = Set{Tuple{Int,Int}}()

    # Detect query relationship
    for (u, rcd1) in enumerate(fgrecords)
        println(rcd1["key"])
        for (v, rcd2) in enumerate(fgrecords)
            rcd1["key"] == rcd2["key"] && continue
            q1s = smartstomol.(haskey(rcd1, "any") ? rcd1["any"] : [rcd1["query"]])
            q2s = smartstomol.(haskey(rcd2, "any") ? rcd2["any"] : [rcd2["query"]])
            # isa
            for q1 in q1s
                for q2 in q2s
                    if hasexactmatch(q1, q2, atommatcher=isaatommatch, bondmatcher=isabondmatch)
                        push!(isaedges, (u, v))
                        println(rcd1["key"], " isa ", rcd2["key"])
                        break
                    end
                end
                (u, v) in isaedges && break
            end
            # has
            for q1 in q1s
                for q2 in q2s
                    if hassubstructmatch(q1, q2, atommatcher=isaatommatch, bondmatcher=isabondmatch)
                        push!(hasedges, (u, v))
                        println(rcd1["key"], " has ", rcd2["key"])
                        break
                    end
                end
                (u, v) in hasedges && break
            end
        end
    end
    # Generate relationship graph
    isagraph = plaindigraph(length(fgrecords), isaedges)
    redisaedges = transitive_reduction(isagraph)
    hasgraph = plaindigraph(length(fgrecords), hasedges)
    redhasedges = transitive_reduction(hasgraph)
    # Add relationship to records
    for e in redisaedges
        (u, v) = getedge(isagraph, e)
        if !haskey(fgrecords[u], "isa")
            fgrecords[u]["isa"] = []
        end
        push!(fgrecords[u]["isa"], fgrecords[v]["key"])
    end
    for e in redhasedges
        (u, v) = getedge(hasgraph, e)
        (u, v) in isaedges && continue
        if !haskey(fgrecords[u], "has")
            fgrecords[u]["has"] = []
        end
        push!(fgrecords[u]["has"], fgrecords[v]["key"])
    end
    # Write records to the file
    YAML.write_file(OUTPUT_FILE, fgrecords)
end


run()