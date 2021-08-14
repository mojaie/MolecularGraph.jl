
const PROJECT_DIR = joinpath(@__DIR__, "..")

using Pkg
Pkg.activate(PROJECT_DIR)

using YAML
using MolecularGraph
using MolecularGraph.Graph
using MolecularGraph.Util

const INPUT_DIR = joinpath(PROJECT_DIR, "./assets/raw/functionalgroup")
const OUTPUT_FILE = joinpath(PROJECT_DIR, "./assets/const/functionalgroup.yaml")


function run()
    # read file directories
    paths = [joinpath(INPUT_DIR, f) for f in readdir(INPUT_DIR)]
    # read list of functional groups from .yaml files
    fgrecords = []
    for p in paths
        for rcd in YAML.load(open(p))
            qs = haskey(rcd, "any") ? rcd["any"] : [rcd["query"]]
            rcd["qmols"] = smartstomol.(qs)
            push!(fgrecords, rcd)
        end
    end
    # Detect query relationships
    isaedges = Set{Tuple{Int,Int}}()
    hasedges = Set{Tuple{Int,Int}}()
    for u in 1:length(fgrecords)
        println(fgrecords[u]["key"])
        for v in 1:length(fgrecords)
            u == v && continue
            # isa
            for qu in fgrecords[u]["qmols"]
                for qv in fgrecords[v]["qmols"]
                    if hasexactmatch(qu, qv)
                        push!(isaedges, (u, v))
                        println(fgrecords[u]["key"], " isa ", fgrecords[v]["key"])
                        break
                    end
                end
                (u, v) in isaedges && break
            end
            # has
            for qu in fgrecords[u]["qmols"]
                for qv in fgrecords[v]["qmols"]
                    if hassubstructmatch(qu, qv, atommatcher=recursiveatommatch)
                        push!(hasedges, (u, v))
                        println(fgrecords[u]["key"], " has ", fgrecords[v]["key"])
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
    for rcd in fgrecords
        delete!(rcd, "qmols")
    end
    # Write records to the file
    YAML.write_file(OUTPUT_FILE, fgrecords)
end


run()