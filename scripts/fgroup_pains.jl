
const PROJECT_DIR = joinpath(@__DIR__, "..")

using Pkg
Pkg.activate(PROJECT_DIR)

using YAML
using MolecularGraph
using MolecularGraph.Graph
using MolecularGraph.Util

const INPUT_DIR = joinpath(PROJECT_DIR, "./assets/raw/functionalgroup")
const INPUT_SA_DIR = joinpath(PROJECT_DIR, "./assets/raw/structuralalerts")
const OUTPUT_FILE = joinpath(PROJECT_DIR, "./assets/const/fgroup_pains.yaml")


function run()
    # read file directories
    paths = [joinpath(INPUT_DIR, f) for f in readdir(INPUT_DIR)]
    # read list of functional groups from .yaml files
    fgrecords = []
    for p in paths
        for rcd in YAML.load(open(p))
            rcd["qmol"] = smartstomol(rcd["query"])
            push!(fgrecords, rcd)
        end
    end
    # PAINS
    for rcd in YAML.load(open(joinpath(INPUT_SA_DIR, "PAINS.yaml")))
        m = smartstomol(rcd["query"])
        m = removehydrogens(m)
        inferatomaromaticity!(m)
        rcd["qmol"] = m
        push!(fgrecords, rcd)
    end
    
    # Detect query relationships
    isaedges = Set{Tuple{Int,Int}}()
    hasedges = Set{Tuple{Int,Int}}()
    for u in 1:length(fgrecords)
        println(fgrecords[u]["key"])
        for v in 1:length(fgrecords)
            u == v && continue
            # isa
            if hasexactmatch(fgrecords[u]["qmol"], fgrecords[v]["qmol"])
                if (v, u) in isaedges
                    pop!(isaedges, (v, u))  # duplicate
                    if !haskey(fgrecords[u], "aliases")
                        fgrecords[u]["aliases"] = []
                    end
                    if !haskey(fgrecords[v], "aliases")
                        fgrecords[v]["aliases"] = []
                    end
                    push!(fgrecords[u]["aliases"], fgrecords[v]["key"])
                    push!(fgrecords[v]["aliases"], fgrecords[u]["key"])
                    println("dupes: ", fgrecords[u]["key"], " ", fgrecords[v]["key"])
                else
                    push!(isaedges, (u, v))
                    println(fgrecords[u]["key"], " isa ", fgrecords[v]["key"])
                end
            end
            # has
            if hassubstructmatch(
                    fgrecords[u]["qmol"], fgrecords[v]["qmol"], atommatcher=recursiveatommatch)
                if (v, u) in hasedges
                    pop!(hasedges, (v, u))  # duplicate
                else
                    push!(hasedges, (u, v))
                    println(fgrecords[u]["key"], " has ", fgrecords[v]["key"])
                end
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
        delete!(rcd, "qmol")
    end
    # Write records to the file
    YAML.write_file(OUTPUT_FILE, fgrecords)
end


run()