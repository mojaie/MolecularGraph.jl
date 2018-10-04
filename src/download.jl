
import HTTP

const RESOURCE_DIR = joinpath(dirname(@__FILE__), "..", "_resources")


function initialize()
    if not isdir(RESOURCE_DIR)
        mkdir(RESOURCE_DIR)
        print("Resource directory created: $(RESOURCE_DIR)")
end


function clear()
    initialize()
    for p in readdir(RESOURCE_DIR)
        rm(joinpath(RESOURCE_DIR, p))
    print("Resource directory is now empty: $(RESOURCE_DIR)")
end


function download(filename::AbstractString, url::AbstractString,
                  decode="utf-8")
    initialize()
    chunksize = 1024 * 1024
    print("Download started: $(url)")
    r = HTTP.request("GET", url)
    println(r.status)
    println(String(r.body))
end


function pubchemsdf(cid::AbstractString)
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
    download(dest)
