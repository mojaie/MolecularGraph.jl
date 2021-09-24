
const PROJECT_DIR = joinpath(@__DIR__, "..")

using Pkg
Pkg.activate(PROJECT_DIR)

using YAML
using MolecularGraph
using MolecularGraph.Graph
using MolecularGraph.Util

const INPUT_DIR = joinpath(PROJECT_DIR, "./assets/raw/functionalgroup")
const INPUT_SA_DIR = joinpath(PROJECT_DIR, "./assets/raw/structuralalerts")
const OUTPUT_FILE = joinpath(PROJECT_DIR, "./assets/const/default_query_relations_new.yaml")


function run()
    # read file directories
    paths = [joinpath(INPUT_DIR, f) for f in readdir(INPUT_DIR)]
    append!(paths, [joinpath(INPUT_SA_DIR, f) for f in readdir(INPUT_SA_DIR)])

    # read list of functional groups from .yaml files
    fgrecords = []
    for p in paths
        basename(p) == ".DS_Store" && continue
        for rcd in YAML.load(open(p))
            m = smartstomol(rcd["query"])
            if basename(p) == "PAINS.yaml"
                m = removehydrogens(m)
                inferatomaromaticity!(m)
            end
            convertnotquery!(m)
            rcd["qmol"] = m
            rcd["source"] = split(basename(p), ".")[1]
            push!(fgrecords, rcd)
        end
    end
    
    # Detect query relationships
    isaedges = Set{Tuple{Int,Int}}()
    hasedges = Set{Tuple{Int,Int}}()
    for (u, v) in combinations(length(fgrecords))
        if u + 1 == v
            @debug "---" fgrecords[u]["key"]
            u % 50 == 0 && println(u, " records processed...")
        end
        umol = fgrecords[u]["qmol"]
        vmol = fgrecords[v]["qmol"]
        ukey = fgrecords[u]["key"]
        vkey = fgrecords[v]["key"]
        # isa
        uv = hasexactmatch(umol, vmol)
        vu = hasexactmatch(vmol, umol)
        if uv && vu  # duplicate
            if !haskey(fgrecords[u], "aliases")
                fgrecords[u]["aliases"] = []
            end
            if !haskey(fgrecords[v], "aliases")
                fgrecords[v]["aliases"] = []
            end
            push!(fgrecords[u]["aliases"], vkey)
            push!(fgrecords[v]["aliases"], ukey)
            @debug "dupes" ukey vkey
        elseif uv
            push!(isaedges, (u, v))
            @debug "isa" ukey vkey
        elseif vu
            push!(isaedges, (v, u))
            @debug "isa" vkey ukey
        else
            # has
            uvsub = hassubstructmatch(umol, vmol)
            vusub = hassubstructmatch(vmol, umol)
            @assert !(uvsub && vusub)
            if uvsub
                push!(hasedges, (u, v))
                @debug "has" ukey vkey
            elseif vusub
                push!(hasedges, (v, u))
                @debug "has" vkey ukey
            end
        end
    end
    # Generate relationship graph
    @assert isempty(intersect(isaedges, hasedges))
    g = plaindigraph(length(fgrecords), union(isaedges, hasedges))
    # Add relationship to records
    for e in transitive_reduction(g)
        (u, v) = getedge(g, e)
        if (u, v) in isaedges
            if !haskey(fgrecords[u], "isa")
                fgrecords[u]["isa"] = []
            end
            push!(fgrecords[u]["isa"], fgrecords[v]["key"])
        elseif (u, v) in hasedges
            if !haskey(fgrecords[u], "has")
                fgrecords[u]["has"] = []
            end
            push!(fgrecords[u]["has"], fgrecords[v]["key"])
        end
    end
    for rcd in fgrecords
        delete!(rcd, "qmol")
    end
    # Write records to the file
    YAML.write_file(OUTPUT_FILE, fgrecords)
end


run()