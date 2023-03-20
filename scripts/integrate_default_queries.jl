
const PROJECT_DIR = joinpath(@__DIR__, "..")

using Pkg
Pkg.activate(PROJECT_DIR)

using YAML
using Graphs
using MolecularGraph

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
            rcd["qmol"] = smartstomol(rcd["query"])
            rcd["source"] = split(basename(p), ".")[1]
            if !haskey(rcd, "name")
                rcd["name"] = rcd["key"]
            end
            push!(fgrecords, rcd)
        end
    end
    
    # Detect query relationships
    dupes = Set{Int}()
    isaedges = Set{Edge{Int}}()
    hasedges = Set{Edge{Int}}()
    for (u, v) in combinations(length(fgrecords))
        if u + 1 == v
            @debug "---" fgrecords[u]["key"]
            u % 50 == 0 && println(u, " records processed...")
        end
        u in dupes && continue
        umol = fgrecords[u]["qmol"]
        vmol = fgrecords[v]["qmol"]
        ukey = fgrecords[u]["key"]
        vkey = fgrecords[v]["key"]
        uvmap = collect(substruct_matches(umol, vmol))
        vumap = collect(substruct_matches(vmol, umol))
        uvsub = !isempty(uvmap)
        vusub = !isempty(vumap)
        uv = uvsub && length(uvmap[1]) == nv(umol)
        vu = vusub && length(vumap[1]) == nv(vmol)
        if uv && vu  # equivalent query
            if !haskey(fgrecords[u], "aliases")
                fgrecords[u]["aliases"] = []
            end
            if !haskey(fgrecords[v], "aliases")
                fgrecords[v]["aliases"] = []
            end
            aliastext = join([fgrecords[v]["source"], fgrecords[v]["name"], fgrecords[v]["query"]], ": ")
            push!(fgrecords[u]["aliases"], aliastext)
            push!(dupes, v)
            @debug "dupes" ukey vkey
        elseif uv  # has_exact_match(umol, vmol)
            push!(isaedges, Edge(u, v))
            @debug "isa" ukey vkey
        elseif vu  # has_exact_match(vmol, umol)
            push!(isaedges, Edge(v, u))
            @debug "isa" vkey ukey
        else
            if uvsub
                push!(hasedges, Edge(u, v))
                @debug "has" ukey vkey
            elseif vusub
                push!(hasedges, Edge(v, u))
                @debug "has" vkey ukey
            end
        end
    end
    # Generate relationship graph
    @assert isempty(intersect(isaedges, hasedges))
    g = SimpleDiGraph(collect(union(isaedges, hasedges)))
    g = transitivereduction(g)
    # Add relationship to records
    for e in edges(g)
        (src(e) in dupes || dst(e) in dupes) && continue
        if e in isaedges
            if !haskey(fgrecords[src(e)], "isa")
                fgrecords[src(e)]["isa"] = []
            end
            push!(fgrecords[src(e)]["isa"], fgrecords[dst(e)]["key"])
        elseif e in hasedges
            if !haskey(fgrecords[src(e)], "has")
                fgrecords[src(e)]["has"] = []
            end
            push!(fgrecords[src(e)]["has"], fgrecords[dst(e)]["key"])
        end
    end
    filtered = []
    for (i, rcd) in enumerate(fgrecords)
        i in dupes && continue
        delete!(rcd, "qmol")
        # sort terms for consistency in diff check
        haskey(rcd, "isa") && (rcd["isa"] = sort(rcd["isa"]))
        haskey(rcd, "has") && (rcd["has"] = sort(rcd["has"]))
        haskey(rcd, "aliases") && (rcd["aliases"] = sort(rcd["aliases"]))
        push!(filtered, rcd)
    end
    # Write records to the file
    YAML.write_file(OUTPUT_FILE, filtered)
end


run()