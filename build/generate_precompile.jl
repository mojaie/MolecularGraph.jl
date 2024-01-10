

# For testing, execute following command at the project dir.
# JULIA_DEBUG=Script julia --project=. -- ./build/generate_precompile.jl


module Script

using MolecularGraph
using JSON

function run()
    op = JSON.json(Dict())
    paclitaxel = JSON.parse(unsafe_string(smilestomol(pointer(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    ), pointer(op))))
    ikey_paclitaxel = inchikey(pointer(JSON.json(paclitaxel)))
    nv_paclitaxel = vertex_count(pointer(JSON.json(paclitaxel)))
    @debug unsafe_string(ikey_paclitaxel)
    @debug nv_paclitaxel

    tc99m = JSON.parse(unsafe_string(smilestomol(pointer(
        "[99Tc+4].[O-]P([O-])(=O)OP([O-])([O-])=O"
    ), pointer(op))))
    ikey_tc99m = inchikey(pointer(JSON.json(tc99m)))
    ne_tc99m = edge_count(pointer(JSON.json(tc99m)))
    @debug unsafe_string(ikey_tc99m)
    @debug ne_tc99m

    acalabrutinib = JSON.parse(unsafe_string(smilestomol(pointer(
        "CC#CC(=O)N1CCC[C@H]1C1=NC(=C2N1C=CN=C2N)C1=CC=C(C=C1)C(=O)NC1=CC=CC=N1"
    ), pointer(op))))
    ikey_acalabrutinib = inchikey(pointer(JSON.json(acalabrutinib)))
    @debug unsafe_string(ikey_acalabrutinib)
    mw_acalabrutinib = standard_weight(pointer(JSON.json(acalabrutinib)))
    @debug mw_acalabrutinib

    nullsmiles = JSON.parse(unsafe_string(smilestomol(pointer(""), pointer(op))))
    ikey_nullsmiles = inchikey(pointer(JSON.json(nullsmiles)))
    @debug unsafe_string(ikey_nullsmiles)

    halide = JSON.parse(unsafe_string(smartstomol(pointer("[#9,#17,#35]c1ccccc1"))))
    notamide = JSON.parse(unsafe_string(smartstomol(pointer(raw"[NX3;H2,H1;!$(NC=O)]"))))
    rotatable = JSON.parse(unsafe_string(smartstomol(pointer(raw"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"))))

    demomol_ = read(open(joinpath(dirname(@__FILE__), "../assets/test/demo.mol")), String)
    demomol = JSON.parse(unsafe_string(sdftomol(pointer(demomol_), pointer(op))))
    mw_demomol = standard_weight(pointer(JSON.json(demomol)))
    @debug mw_demomol
    svg_demomol = drawsvg(pointer(JSON.json(demomol)))
    @debug unsafe_string(svg_demomol)
    @debug unsafe_string(sdftosvg(pointer(demomol_)))
    buf_size = Sys.isapple() ? 1_000_000 : 200_000  # MacOS seg fault workaround
    dst = Array{UInt8}(undef, (buf_size,))  # may be sufficient for 1000x1000 structure images
    p = pointer(dst)
    imgsize = drawpng(p, pointer(JSON.json(demomol)), UInt32(1000), UInt32(1000))
    o = IOBuffer(write=true)
    unsafe_write(o, p, imgsize)
    # f = open("test.png", "w")
    # write(f, o.data)
    # close(f)
    @debug length(o.data)

    nullmol_ = read(open(joinpath(dirname(@__FILE__), "../assets/test/null.mol")), String)
    nullmol = JSON.parse(unsafe_string(sdftomol(pointer(nullmol_), pointer(op))))
    mw_nullmol = standard_weight(pointer(JSON.json(nullmol)))
    @debug mw_nullmol

    tenatoprazole = JSON.parse(unsafe_string(smilestomol(pointer(
        "Cc1c(OC)c(C)cnc1CS(=O)c2[nH]c3ccc(OC)nc3n2"
    ), pointer(op))))
    @debug unsafe_string(smilestosvg(pointer("Cc1c(OC)c(C)cnc1CS(=O)c2[nH]c3ccc(OC)nc3n2")))
    query1 = JSON.parse(unsafe_string(smilestomol(pointer("c1ncccc1N"), pointer(op))))
    query2 = JSON.parse(unsafe_string(smartstomol(pointer("[r6R2][r6R2]"))))
    halide2 = JSON.parse(unsafe_string(smilestomol(pointer("[F]c1ccccc1"), pointer(op))))

    @debug has_exact_match(
        pointer(JSON.json(nullsmiles)),
        pointer(JSON.json(nullmol)),
        pointer(op))
    @debug has_exact_match(
        pointer(JSON.json(halide2)),
        pointer(JSON.json(halide)),
        pointer(op))
    @debug has_exact_match(
        pointer(JSON.json(tenatoprazole)),
        pointer(JSON.json(tenatoprazole)),
        pointer(op))

    @debug has_substruct_match(
        pointer(JSON.json(nullsmiles)),
        pointer(JSON.json(nullmol)),
        pointer(op))
    @debug has_substruct_match(
        pointer(JSON.json(tenatoprazole)),
        pointer(JSON.json(query1)),
        pointer(op))
    @debug has_substruct_match(
        pointer(JSON.json(paclitaxel)),
        pointer(JSON.json(query2)),
        pointer(op))

    furosemide = JSON.parse(unsafe_string(smilestomol(pointer(
        "O=S(=O)(N)c1c(Cl)cc(c(C(=O)O)c1)NCc2occc2"), pointer(op))))
    acetazolamide = JSON.parse(unsafe_string(smilestomol(pointer(
        "O=S(=O)(c1nnc(s1)NC(=O)C)N"), pointer(JSON.json(Dict())))))
    bortezomib = JSON.parse(unsafe_string(smilestomol(pointer(
        "O=C(N[C@H](C(=O)N[C@H](B(O)O)CC(C)C)Cc1ccccc1)c2nccnc2"), pointer(op))))
    tofacitinib = JSON.parse(unsafe_string(smilestomol(pointer(
        raw"CC2CCN(C(=O)CC#N)CC2N(C)c3ncnc1[nH]ccc13"), pointer(op))))
    @debug tdmcis_size(
        pointer(JSON.json(nullmol)),
        pointer(JSON.json(nullsmiles)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_size(
        pointer(JSON.json(furosemide)),
        pointer(JSON.json(acetazolamide)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_tanimoto(
        pointer(JSON.json(furosemide)),
        pointer(JSON.json(acetazolamide)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_dist(
        pointer(JSON.json(bortezomib)),
        pointer(JSON.json(tofacitinib)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmcis_gls(
        pointer(JSON.json(bortezomib)),
        pointer(JSON.json(tofacitinib)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))

    @debug tdmces_size(
        pointer(JSON.json(nullsmiles)),
        pointer(JSON.json(nullmol)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmces_tanimoto(
        pointer(JSON.json(furosemide)),
        pointer(JSON.json(acetazolamide)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmces_dist(
        pointer(JSON.json(furosemide)),
        pointer(JSON.json(acetazolamide)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tdmces_gls(
        pointer(JSON.json(bortezomib)),
        pointer(JSON.json(tofacitinib)),
        pointer(JSON.json(Dict("diameter" => 8))))

    @debug JSON.parse(unsafe_string(tdmces_gls_batch(
        pointer(JSON.json([
            [1,2,3],
            [2,3,4],
            [furosemide, bortezomib, tofacitinib],
            [bortezomib, tofacitinib, acetazolamide],
            Dict("diameter" => 8), 0.1
        ])))))
end

end

Script.run()