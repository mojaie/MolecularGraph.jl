

# For testing, execute following command at the project dir.
# JULIA_DEBUG=Script julia --project=. -- ./build/generate_precompile.jl


module Script

using Base: unsafe_convert
using Base64
using JSON
using MolecularGraph

function run()
    op = JSON.json(Dict())
    paclitaxel_str = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    paclitaxel = unsafe_string(smilestomol(
        unsafe_convert(Cstring, paclitaxel_str), unsafe_convert(Cstring, op)))
    @debug "nv(paclitaxel)" vertex_count(unsafe_convert(Cstring, paclitaxel))
    @debug "inchikey(paclitaxel)" unsafe_string(inchikey(unsafe_convert(Cstring, paclitaxel)))

    tc99m_str = "[99Tc+4].[O-]P([O-])(=O)OP([O-])([O-])=O"
    tc99m = unsafe_string(smilestomol(
        unsafe_convert(Cstring, tc99m_str), unsafe_convert(Cstring, op)))
    @debug "ne(tc99m)" edge_count(unsafe_convert(Cstring, tc99m))
    @debug "inchikey(tc99m)" unsafe_string(inchikey(unsafe_convert(Cstring, tc99m)))

    acalabrutinib_str = "CC#CC(=O)N1CCC[C@H]1C1=NC(=C2N1C=CN=C2N)C1=CC=C(C=C1)C(=O)NC1=CC=CC=N1"
    acalabrutinib = unsafe_string(smilestomol(
        unsafe_convert(Cstring, acalabrutinib_str), unsafe_convert(Cstring, op)))
    @debug "standard_weight(acalabrutinib)" standard_weight(unsafe_convert(Cstring, acalabrutinib))
    @debug "inchikey(acalabrutinib)" unsafe_string(inchikey(unsafe_convert(Cstring, acalabrutinib)))

    nullsmiles_str = ""
    nullsmiles = unsafe_string(smilestomol(
        unsafe_convert(Cstring, nullsmiles_str), unsafe_convert(Cstring, op)))
    @debug "inchikey(nullsmiles)" unsafe_string(inchikey(unsafe_convert(Cstring, nullsmiles)))

    demomol_sdf = read(open(joinpath(dirname(@__FILE__), "../assets/test/demo.mol")), String)
    demomol = unsafe_string(sdftomol(
        unsafe_convert(Cstring, demomol_sdf), unsafe_convert(Cstring, op)))
    @debug "standard_weight(demomol)" standard_weight(unsafe_convert(Cstring, demomol))
    @debug "length(drawsvg(demomol))" length(unsafe_string(drawsvg(unsafe_convert(Cstring, demomol))))
    str = drawpng(unsafe_convert(Cstring, demomol), UInt32(1000), UInt32(1000))
    img = base64decode(unsafe_string(str))
    # f = open("test.png", "w")
    # write(f, img)
    # close(f)
    @debug "length(drawpng(io, demomol, 1000, 1000))" length(img)
    @debug "length(molblock(demomol))" length(unsafe_string(molblock(unsafe_convert(Cstring, demomol))))
    @debug "length(sdfmolblock(demomol))" length(unsafe_string(sdfmolblock(unsafe_convert(Cstring, demomol))))

    nullmol_sdf = read(open(joinpath(dirname(@__FILE__), "../assets/test/null.mol")), String)
    nullmol = unsafe_string(sdftomol(
        unsafe_convert(Cstring, nullmol_sdf), unsafe_convert(Cstring, op)))
    @debug "standard_weight(nullsdf)" standard_weight(unsafe_convert(Cstring, nullmol))
    errormol_sdf = read(open(joinpath(dirname(@__FILE__), "../assets/test/error.mol")), String)
    errormol = unsafe_string(sdftomol(
        unsafe_convert(Cstring, errormol_sdf), unsafe_convert(Cstring, op)))
    @debug "standard_weight(errormol)" standard_weight(unsafe_convert(Cstring, errormol))
    @debug "errormol error msg" MolGraph(errormol).gprops[:error_sdfilereader]

    notamide = unsafe_string(smartstomol(
        unsafe_convert(Cstring, raw"[NX3;H2,H1;!$(NC=O)]")))
    notamide2 = unsafe_string(smartstomol(
        unsafe_convert(Cstring, raw"[NX3;H2,H1]")))
    rotatable = unsafe_string(smartstomol(
        unsafe_convert(Cstring, raw"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")))
    @debug "notamide rotatable" has_exact_match(
        unsafe_convert(Cstring, notamide), unsafe_convert(Cstring, rotatable),
        unsafe_convert(Cstring, op))
    @debug "notamide notamide2" has_substruct_match(
        unsafe_convert(Cstring, notamide), unsafe_convert(Cstring, notamide2),
        unsafe_convert(Cstring, op))

    @debug "nullsmiles nullmol" has_exact_match(
        unsafe_convert(Cstring, nullsmiles), unsafe_convert(Cstring, nullmol),
        unsafe_convert(Cstring, op))
    @debug "nullsmiles nullmol" has_substruct_match(
        unsafe_convert(Cstring, nullsmiles), unsafe_convert(Cstring, nullmol),
        unsafe_convert(Cstring, op))

    halide = unsafe_string(smartstomol(unsafe_convert(Cstring, "[#9,#17,#35]c1ccccc1")))
    halide2 = unsafe_string(smilestomol(
        unsafe_convert(Cstring, "[F]c1ccccc1"), unsafe_convert(Cstring, op)))
    @debug "halide2 halide" has_exact_match(
        unsafe_convert(Cstring, halide2), unsafe_convert(Cstring, halide),
        unsafe_convert(Cstring, op))

    tenatoprazole = unsafe_string(smilestomol(
        unsafe_convert(Cstring, "Cc1c(OC)c(C)cnc1CS(=O)c2[nH]c3ccc(OC)nc3n2"),
        unsafe_convert(Cstring, op)
    ))
    query1 = unsafe_string(smilestomol(
        unsafe_convert(Cstring, "c1ncccc1N"), unsafe_convert(Cstring, op)))
    query2 = unsafe_string(smartstomol(unsafe_convert(Cstring, "[r6R2][r6R2]")))
    @debug "tenatoprazole tenatoprazole" has_exact_match(
        unsafe_convert(Cstring, tenatoprazole), unsafe_convert(Cstring, tenatoprazole),
        unsafe_convert(Cstring, op))
    @debug "tenatoprazole query1" has_substruct_match(
        unsafe_convert(Cstring, tenatoprazole), unsafe_convert(Cstring, query1),
        unsafe_convert(Cstring, op))
    @debug "paclitaxel query2" has_substruct_match(
        unsafe_convert(Cstring, paclitaxel), unsafe_convert(Cstring, query2),
        unsafe_convert(Cstring, op))

    furosemide = unsafe_string(smilestomol(
        unsafe_convert(Cstring, "O=S(=O)(N)c1c(Cl)cc(c(C(=O)O)c1)NCc2occc2"),
        unsafe_convert(Cstring, op)
    ))
    acetazolamide = unsafe_string(smilestomol(
        unsafe_convert(Cstring, "O=S(=O)(c1nnc(s1)NC(=O)C)N"),
        unsafe_convert(Cstring, op)
    ))
    bortezomib = unsafe_string(smilestomol(
        unsafe_convert(Cstring, "O=C(N[C@H](C(=O)N[C@H](B(O)O)CC(C)C)Cc1ccccc1)c2nccnc2"),
        unsafe_convert(Cstring, op)
    ))
    tofacitinib = unsafe_string(smilestomol(
        unsafe_convert(Cstring, raw"CC2CCN(C(=O)CC#N)CC2N(C)c3ncnc1[nH]ccc13"),
        unsafe_convert(Cstring, op)
    ))
    @debug tdmcis_size(
        unsafe_convert(Cstring, nullmol), unsafe_convert(Cstring, nullsmiles),
        unsafe_convert(Cstring, JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_size(
        unsafe_convert(Cstring, furosemide), unsafe_convert(Cstring, acetazolamide),
        unsafe_convert(Cstring, JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_gls(
        unsafe_convert(Cstring, bortezomib), unsafe_convert(Cstring, tofacitinib),
        unsafe_convert(Cstring, JSON.json(Dict("diameter" => 8, "tolerance" => 1))))

    @debug tdmces_size(
        unsafe_convert(Cstring, nullsmiles), unsafe_convert(Cstring, nullmol),
        unsafe_convert(Cstring, JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmces_size(
        unsafe_convert(Cstring, acetazolamide), unsafe_convert(Cstring, furosemide),
        unsafe_convert(Cstring, JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmces_gls(
        unsafe_convert(Cstring, bortezomib), unsafe_convert(Cstring, tofacitinib),
        unsafe_convert(Cstring, JSON.json(Dict("diameter" => 8))))
end

end

Script.run()