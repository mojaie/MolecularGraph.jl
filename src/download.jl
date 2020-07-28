#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# TODO: deprecated

export
    clear,
    fetchfile,
    fetchresource,
    pubchemsdf,
    molsfortest,
    resource_dir,
    pubchem_dir


resource_dir = joinpath(dirname(@__FILE__), "..", "_resources")
pubchem_dir = joinpath(resource_dir, "PubChem")

function initialize_resource_dir()
    if !isdir(resource_dir)
        mkdir(resource_dir)
        println("Resource directory created: $(resource_dir)")
    end
end


function clear()
    initialize_resource_dir()
    for p in readdir(resource_dir)
        rm(p, recursive=true)
    end
    println("Resource directory is now empty: $(resource_dir)")
end


function fetchfile(url, dest)
    if isfile(dest)
        println("file: $(basename(dest)) already exists")
        return
    end
    println("Downloading: $(url)")
    download(url, dest)
end


function fetchresource(url, filename)
    initialize_resource_dir()
    dest = joinpath(resource_dir, filename)
    fetchfile(url, dest)
end


function pubchemsdf(cid, name)
    if !isdir(pubchem_dir)
        initialize_resource_dir()
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
        initialize_resource_dir()
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
