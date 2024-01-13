

# For testing, execute following command at the project dir.
# JULIA_DEBUG=Script julia --project=. -- ./build/generate_precompile.jl


module Script

using MolecularGraph
using JSON

function run()
    op = JSON.json(Dict())
    paclitaxel = unsafe_string(smilestomol(pointer(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    ), pointer(op)))
    @debug "nv(paclitaxel)" vertex_count(pointer(paclitaxel))
    @debug "inchikey(paclitaxel)" unsafe_string(inchikey(pointer(paclitaxel)))

    tc99m = unsafe_string(smilestomol(pointer(
        "[99Tc+4].[O-]P([O-])(=O)OP([O-])([O-])=O"
    ), pointer(op)))
    @debug "ne(tc99m)" edge_count(pointer(tc99m))
    @debug "inchikey(tc99m)" unsafe_string(inchikey(pointer(tc99m)))

    acalabrutinib = unsafe_string(smilestomol(pointer(
        "CC#CC(=O)N1CCC[C@H]1C1=NC(=C2N1C=CN=C2N)C1=CC=C(C=C1)C(=O)NC1=CC=CC=N1"
    ), pointer(op)))
    @debug "standard_weight(acalabrutinib)" standard_weight(pointer(acalabrutinib))
    @debug "inchikey(acalabrutinib)" unsafe_string(inchikey(pointer(acalabrutinib)))

    nullsmiles = unsafe_string(smilestomol(pointer(""), pointer(op)))
    @debug "inchikey(nullsmiles)" unsafe_string(inchikey(pointer(nullsmiles)))

    halide = unsafe_string(smartstomol(pointer("[#9,#17,#35]c1ccccc1")))
    notamide = unsafe_string(smartstomol(pointer(raw"[NX3;H2,H1;!$(NC=O)]")))
    rotatable = unsafe_string(smartstomol(pointer(raw"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")))

    demomol_ = read(open(joinpath(dirname(@__FILE__), "../assets/test/demo.mol")), String)
    demomol = unsafe_string(sdftomol(pointer(demomol_), pointer(op)))
    dst = Array{UInt8}(undef, (200_000,))  # may be sufficient for 1000x1000 structure images
    # GC.@preserve demomol dst begin
    @debug "standard_weight(demomol)" standard_weight(pointer(demomol))
    @debug "length(drawsvg(demomol))" length(unsafe_string(drawsvg(pointer(demomol))))
    p = pointer(dst)
    imgsize = drawpng(p, pointer(demomol), UInt32(1000), UInt32(1000))
    o = IOBuffer(write=true)
    unsafe_write(o, p, imgsize)
    # f = open("test.png", "w")
    # write(f, o.data)
    # close(f)
    @debug "length(drawpng(io, demomol, 1000, 1000))" length(o.data)
    # end

    nullmol_ = read(open(joinpath(dirname(@__FILE__), "../assets/test/null.mol")), String)
    nullmol = unsafe_string(sdftomol(pointer(nullmol_), pointer(op)))
    @debug "standard_weight(nullsdf)" standard_weight(pointer(nullmol))
    errormol_ = read(open(joinpath(dirname(@__FILE__), "../assets/test/error.mol")), String)
    errormol = unsafe_string(sdftomol(pointer(errormol_), pointer(op)))
    @debug "standard_weight(errormol)" standard_weight(pointer(errormol))
    @debug "errormol error msg" MolGraph(errormol).gprops[:error_sdfilereader]

    tenatoprazole = unsafe_string(smilestomol(pointer(
        "Cc1c(OC)c(C)cnc1CS(=O)c2[nH]c3ccc(OC)nc3n2"
    ), pointer(op)))
    query1 = unsafe_string(smilestomol(pointer("c1ncccc1N"), pointer(op)))
    query2 = unsafe_string(smartstomol(pointer("[r6R2][r6R2]")))
    halide2 = unsafe_string(smilestomol(pointer("[F]c1ccccc1"), pointer(op)))

    @debug "nullsmiles nullmol" has_exact_match(
        pointer(nullsmiles), pointer(nullmol), pointer(op))
    @debug "halide2 halide" has_exact_match(
        pointer(halide2), pointer(halide), pointer(op))
    @debug "tenatoprazole tenatoprazole" has_exact_match(
        pointer(tenatoprazole), pointer(tenatoprazole), pointer(op))

    @debug "nullsmiles nullmol" has_substruct_match(
        pointer(nullsmiles), pointer(nullmol), pointer(op))
    @debug "tenatoprazole query1" has_substruct_match(
        pointer(tenatoprazole), pointer(query1), pointer(op))
    @debug "tenatoprazole query2" has_substruct_match(
        pointer(paclitaxel), pointer(query2), pointer(op))

    furosemide = unsafe_string(smilestomol(pointer(
        "O=S(=O)(N)c1c(Cl)cc(c(C(=O)O)c1)NCc2occc2"), pointer(op)))
    acetazolamide = unsafe_string(smilestomol(pointer(
        "O=S(=O)(c1nnc(s1)NC(=O)C)N"), pointer(op)))
    bortezomib = unsafe_string(smilestomol(pointer(
        "O=C(N[C@H](C(=O)N[C@H](B(O)O)CC(C)C)Cc1ccccc1)c2nccnc2"), pointer(op)))
    tofacitinib = unsafe_string(smilestomol(pointer(
        raw"CC2CCN(C(=O)CC#N)CC2N(C)c3ncnc1[nH]ccc13"), pointer(op)))
    @debug tdmcis_size(
        pointer(nullmol), pointer(nullsmiles),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_size(
        pointer(furosemide), pointer(acetazolamide),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_tanimoto(
        pointer(furosemide), pointer(acetazolamide),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_dist(
        pointer(bortezomib), pointer(tofacitinib),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_gls(
        pointer(bortezomib), pointer(tofacitinib),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))

    @debug tdmces_size(
        pointer(nullsmiles), pointer(nullmol),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmces_tanimoto(
        pointer(furosemide), pointer(acetazolamide),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmces_dist(
        pointer(furosemide), pointer(acetazolamide),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmces_gls(
        pointer(bortezomib), pointer(tofacitinib),
        pointer(JSON.json(Dict("diameter" => 8))))
end

end

Script.run()