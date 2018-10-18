
export
    clear,
    pubchemsdf,
    molsfortest

import HTTP

const RESOURCE_DIR = joinpath(dirname(@__FILE__), "..", "_resources")


function initialize()
    if !isdir(RESOURCE_DIR)
        mkdir(RESOURCE_DIR)
        println("Resource directory created: $(RESOURCE_DIR)")
    end
end


function clear()
    initialize()
    for p in readdir(RESOURCE_DIR)
        rm(p, recursive=true)
    end
    println("Resource directory is now empty: $(RESOURCE_DIR)")
end


function pubchemsdf(cid::AbstractString, name::AbstractString)
    pubchemdir = joinpath(RESOURCE_DIR, "PubChem")
    if !isdir(pubchemdir)
        initialize()
        mkdir(pubchemdir)
        println("PubChem data directory created: $(pubchemdir)")
    end
    dest = joinpath(RESOURCE_DIR, "PubChem", "$(name).mol")
    if isfile(dest)
        println("file: $(name).mol already exists")
        return
    end
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(cid)/SDF"
    println("Downloading: $(url)")
    download(url, dest)
end


function molsfortest()
    pubchemfilelist = [
        ("Buckminsterfullerene", 123591)
    ]
    testdir = joinpath(RESOURCE_DIR, "testmols")
    if !isdir(testdir)
        initialize()
        mkdir(testdir)
        println("Test data directory created: $(testdir)")
    end
    for (name, id) in pubchemfilelist
        dest = joinpath(testdir, "$(name).sdf")
        if isfile(dest)
            println("file: $(name).mol already exists")
            return
        end
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$(id)/SDF"
        println("Downloading: $(url)")
        download(url, dest)
    end
    println("Finished: $(testdir)")
end
