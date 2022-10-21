
using JSON
export smilesmoldata, sdfmoldata, structmatch


record(mol::GraphMol) = Dict(
    "inchikey" => inchikey(mol),
    "structure" => todict(mol),
    "svg" => drawsvg(mol, 200, 200),
    "formula" => molecularformula(mol),
    "heavyatom" => heavyatomcount(mol),
    "mw" => standardweight(Float64, mol),
    "logp" => wclogp(mol),
    "donor" => hdonorcount(mol),
    "acceptor" => hacceptorcount(mol),
    "rotatable" => rotatablecount(mol)
)

Base.@ccallable function smilesmoldata(smiles::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = smilestomol(unsafe_string(smiles))
        length(Graph.connectedcomponents(mol)) == 1 || return
        precalculate!(mol)
        buf = IOBuffer(write=true)
        JSON.print(buf, record(mol))
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function sdfmoldata(sdf::Ptr{UInt8})::Ptr{UInt8}
    return try
        mol = sdftomol(split(unsafe_string(sdf), "\n"))
        precalculate!(mol)
        buf = IOBuffer(write=true)
        JSON.print(buf, record(mol))
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

Base.@ccallable function structmatch(jsondump::Ptr{UInt8})::Ptr{UInt8}
    return try
        query = JSON.parse(unsafe_string(jsondump))
        fr = graphmol(query["first"]["structure"])
        sc = graphmol(query["second"]["structure"])
        matchtype = get("matchtype", symbol(query["params"]), :substruct)
        ignoreh = get("ignoreHydrogen", query["params"], true)
        superstr = get("superstructure", query["params"], false)
        if ignoreh
            fr = graphmol(removehydrogen(fr, all=true))
            sc = graphmol(removehydrogen(sc, all=true))
        end
        if superstr
            issub = isstructmatch(sc, fr, matchtype)
        else
            issub = isstructmatch(fr, sc, matchtype)
        end
        res = Dict(
            "first" => query["first"]["id"],
            "second" => query["second"]["id"],
            "result" => issub
        )
        buf = IOBuffer(write=true)
        JSON.print(buf, res)
        return pointer(buf.data)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    smiles = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    println(unsafe_string(smilesmoldata(pointer(smiles))))
    fpath = joinpath(dirname(@__FILE__), "../assets/test/demo.mol")
    sdf = read(open(fpath), String)
    println(unsafe_string(sdfmoldata(pointer(sdf))))
end

