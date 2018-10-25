#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    clear,
    pubchemsdf,
    molsfortest,
    resource_dir,
    pubchem_dir


resource_dir = joinpath(dirname(@__FILE__), "..", "_resources")
pubchem_dir = joinpath(resource_dir, "PubChem")

function initialize()
    if !isdir(resource_dir)
        mkdir(resource_dir)
        println("Resource directory created: $(resource_dir)")
    end
end


function clear()
    initialize()
    for p in readdir(resource_dir)
        rm(p, recursive=true)
    end
    println("Resource directory is now empty: $(resource_dir)")
end


function pubchemsdf(cid::AbstractString, name::AbstractString)
    if !isdir(pubchem_dir)
        initialize()
        mkdir(pubchem_dir)
        println("PubChem data directory created: $(pubchem_dir)")
    end
    dest = joinpath(resource_dir, "PubChem", "$(name).mol")
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
    testdir = joinpath(resource_dir, "testmols")
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
