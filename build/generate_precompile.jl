

# For testing, execute following command at the project dir.
# JULIA_DEBUG=Script julia --project=. -- ./build/generate_precompile.jl


module Script

using MolecularGraph
using JSON

function run()
    paclitaxel = JSON.parse(unsafe_string(smilestomol(pointer(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    ))))
    ikey_paclitaxel = inchikey(pointer(JSON.json(paclitaxel)))
    @debug unsafe_string(ikey_paclitaxel)

    tc99m = JSON.parse(unsafe_string(smilestomol(pointer(
        "[99Tc+4].[O-]P([O-])(=O)OP([O-])([O-])=O"
    ))))
    ikey_tc99m = inchikey(pointer(JSON.json(tc99m)))
    @debug unsafe_string(ikey_tc99m)

    acalabrutinib = JSON.parse(unsafe_string(smilestomol(pointer(
        "CC#CC(=O)N1CCC[C@H]1C1=NC(=C2N1C=CN=C2N)C1=CC=C(C=C1)C(=O)NC1=CC=CC=N1"
    ))))
    ikey_acalabrutinib = inchikey(pointer(JSON.json(acalabrutinib)))
    @debug unsafe_string(ikey_acalabrutinib)
    mw_acalabrutinib = standard_weight(pointer(JSON.json(acalabrutinib)))
    @debug mw_acalabrutinib

    nullsmiles = JSON.parse(unsafe_string(smilestomol(pointer(""))))
    ikey_nullsmiles = inchikey(pointer(JSON.json(nullsmiles)))
    @debug unsafe_string(ikey_nullsmiles)

    halide = JSON.parse(unsafe_string(smartstomol(pointer(
        "[#9,#17,#35]c1ccccc1"
    ))))
    notamide = JSON.parse(unsafe_string(smartstomol(pointer(
        raw"[NX3;H2,H1;!$(NC=O)]"
    ))))
    rotatable = JSON.parse(unsafe_string(smartstomol(pointer(
        raw"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"
    ))))

    demomol_ = read(open(joinpath(dirname(@__FILE__), "../assets/test/demo.mol")), String)
    demomol = JSON.parse(unsafe_string(sdftomol(pointer(demomol_))))
    mw_demomol = standard_weight(pointer(JSON.json(demomol)))
    @debug mw_demomol

    nullmol_ = read(open(joinpath(dirname(@__FILE__), "../assets/test/null.mol")), String)
    nullmol = JSON.parse(unsafe_string(sdftomol(pointer(nullmol_))))
    mw_nullmol = standard_weight(pointer(JSON.json(nullmol)))
    @debug mw_nullmol

    tenatoprazole = JSON.parse(unsafe_string(smilestomol(pointer(
        "Cc1c(OC)c(C)cnc1CS(=O)c2[nH]c3ccc(OC)nc3n2"
    ))))
    query1 = JSON.parse(unsafe_string(smilestomol(pointer("c1ncccc1N"))))
    query2 = JSON.parse(unsafe_string(smartstomol(pointer("[r6R2][r6R2]"))))
    halide2 = JSON.parse(unsafe_string(smilestomol(pointer("[F]c1ccccc1"))))

    @debug has_exact_match(
        pointer(JSON.json(nullsmiles)),
        pointer(JSON.json(nullmol)),
        pointer(JSON.json(Dict())))
    @debug has_exact_match(
        pointer(JSON.json(halide2)),
        pointer(JSON.json(halide)),
        pointer(JSON.json(Dict())))
    @debug has_exact_match(
        pointer(JSON.json(tenatoprazole)),
        pointer(JSON.json(tenatoprazole)),
        pointer(JSON.json(Dict())))

    @debug has_substruct_match(
        pointer(JSON.json(nullsmiles)),
        pointer(JSON.json(nullmol)),
        pointer(JSON.json(Dict())))
    @debug has_substruct_match(
        pointer(JSON.json(tenatoprazole)),
        pointer(JSON.json(query1)),
        pointer(JSON.json(Dict())))
    @debug has_substruct_match(
        pointer(JSON.json(paclitaxel)),
        pointer(JSON.json(query2)),
        pointer(JSON.json(Dict())))

    furosemide = JSON.parse(unsafe_string(smilestomol(pointer("O=S(=O)(N)c1c(Cl)cc(c(C(=O)O)c1)NCc2occc2"))))
    acetazolamide = JSON.parse(unsafe_string(smilestomol(pointer("O=S(=O)(c1nnc(s1)NC(=O)C)N"))))
    bortezomib = JSON.parse(unsafe_string(smilestomol(pointer("O=C(N[C@H](C(=O)N[C@H](B(O)O)CC(C)C)Cc1ccccc1)c2nccnc2"))))
    tofacitinib = JSON.parse(unsafe_string(smilestomol(pointer(raw"CC2CCN(C(=O)CC#N)CC2N(C)c3ncnc1[nH]ccc13"))))
    @debug tcmcis(
        pointer(JSON.json(nullmol)),
        pointer(JSON.json(nullsmiles)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tcmcis(
        pointer(JSON.json(furosemide)),
        pointer(JSON.json(acetazolamide)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tcmcis(
        pointer(JSON.json(bortezomib)),
        pointer(JSON.json(tofacitinib)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tcmces(
        pointer(JSON.json(nullsmiles)),
        pointer(JSON.json(nullmol)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tcmces(
        pointer(JSON.json(furosemide)),
        pointer(JSON.json(acetazolamide)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
    @debug tcmces(
        pointer(JSON.json(bortezomib)),
        pointer(JSON.json(tofacitinib)),
        pointer(JSON.json(Dict("diameter" => 8, "tolerance" => 1))))
end

end

Script.run()