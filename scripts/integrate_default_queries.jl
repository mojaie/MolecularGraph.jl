
const PROJECT_DIR = joinpath(@__DIR__, "..")

using Pkg
Pkg.activate(PROJECT_DIR)

using YAML
using Graphs
using MolecularGraph

const INPUT_DIR = joinpath(PROJECT_DIR, "./assets/raw_queries/default")
const INPUT_SA_DIR = joinpath(PROJECT_DIR, "./assets/raw_queries/ChEMBL")
const OUTPUT_FILE = joinpath(PROJECT_DIR, "./assets/const/default_queries_new.yaml")


function run()
    # read file directories
    paths = [joinpath(INPUT_DIR, f) for f in readdir(INPUT_DIR)]
    append!(paths, [joinpath(INPUT_SA_DIR, f) for f in readdir(INPUT_SA_DIR)])

    # read list of functional groups from .yaml files
    vprops = Dict{String,Any}[]
    for p in paths
        basename(p) == ".DS_Store" && continue
        for rcd in YAML.load(open(p))
            rcd["qmol"] = smartstomol(rcd["query"])
            rcd["source"] = split(basename(p), ".")[1]
            haskey(rcd, "name") || (rcd["name"] = rcd["key"])
            push!(vprops, rcd)
        end
    end
    
    # Detect query relationships
    isaedges = Set{Edge{Int}}()
    hasedges = Set{Edge{Int}}()
    skip = [  # has atoms with too complicated queries
        "str_alert_221", "str_alert_222", "str_alert_247", "str_alert_1003",
        "str_alert_246", "str_alert_381"
    ]
    for (u, v) in combinations(length(vprops))
        ukey = vprops[u]["key"]
        vkey = vprops[v]["key"]
        (ukey in skip || vkey in skip) && continue
        if u + 1 == v
            @info "--- $(ukey)"
            u % 50 == 0 && @info "===== $(u) records processed..."
        end
        # @info fgrecords[u]["key"] fgrecords[v]["key"]
        umol = vprops[u]["qmol"]
        vmol = vprops[v]["qmol"]
        uvmap = collect(substruct_matches(umol, vmol))
        vumap = collect(substruct_matches(vmol, umol))
        uvsub = !isempty(uvmap)
        vusub = !isempty(vumap)
        uv = uvsub && length(uvmap[1]) == nv(umol)
        vu = vusub && length(vumap[1]) == nv(vmol)
        if uv && vu  # equivalent query
            if !haskey(vprops[u], "aliases")
                vprops[u]["aliases"] = []
            end
            if !haskey(vprops[v], "aliases")
                vprops[v]["aliases"] = []
            end
            push!(vprops[u]["aliases"], vkey)
            push!(vprops[v]["aliases"], ukey)
            # @info "dupes" ukey vkey
        end
        if uv  # has_exact_match(umol, vmol)
            haskey(vprops[u], "isa") || (vprops[u]["isa"] = [])
            push!(vprops[u]["isa"], vkey)
            # @info "isa" ukey vkey
        elseif vu  # has_exact_match(vmol, umol)
            haskey(vprops[v], "isa") || (vprops[v]["isa"] = [])
            push!(vprops[v]["isa"], ukey)
            # @info "isa" vkey ukey
        elseif uvsub
            haskey(vprops[u], "has") || (vprops[u]["has"] = [])
            push!(vprops[u]["has"], vkey)
            # @info "has" ukey vkey
        elseif vusub
            haskey(vprops[v], "has") || (vprops[v]["has"] = [])
            push!(vprops[v]["has"], ukey)
            # @info "has" vkey ukey
        end
    end
    # cleanup
    for (i, rcd) in enumerate(vprops)
        delete!(rcd, "qmol")
        # sort terms for consistency in diff check
        haskey(rcd, "isa") && (vprops[i]["isa"] = sort(rcd["isa"]))
        haskey(rcd, "has") && (vprops[i]["has"] = sort(rcd["has"]))
        haskey(rcd, "aliases") && (vprops[i]["aliases"] = sort(rcd["aliases"]))
    end
    # Write records to the file
    YAML.write_file(OUTPUT_FILE, vprops)
end


run()