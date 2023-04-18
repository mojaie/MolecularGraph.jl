

# For testing, execute following command at the project dir.
# JULIA_DEBUG=Script julia --project=. -- ./build/generate_precompile.jl


module Script

using MolecularGraph
using JSON

function run()
    smiles = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    mol1 = JSON.parse(unsafe_string(smilestomol(pointer(smiles))))
    # mol1q = smartstomol(pointer(smiles)) TODO: serialize SMARTS query?
    ikey1 = inchikey(pointer(JSON.json(mol1)))
    @debug unsafe_string(ikey1)

    fpath = joinpath(dirname(@__FILE__), "../assets/test/demo.mol")
    sdf = read(open(fpath), String)
    mol2 = JSON.parse(unsafe_string(sdftomol(pointer(sdf))))
    mw2 = standard_weight(pointer(JSON.json(mol2)))
    @debug mw2

    res = has_exact_match(
        pointer(JSON.json(mol1)),
        pointer(JSON.json(mol2)),
        pointer(JSON.json(Dict())))
    @debug res
    res = has_substruct_match(
        pointer(JSON.json(mol1)),
        pointer(JSON.json(mol2)),
        pointer(JSON.json(Dict())))
    @debug res
    res = tcmcis(
        pointer(JSON.json(mol1)),
        pointer(JSON.json(mol2)),
        pointer(JSON.json(Dict("diameter" => 2, "tolerance" => 1))))
    @debug res
    res = tcmces(
        pointer(JSON.json(mol1)),
        pointer(JSON.json(mol2)),
        pointer(JSON.json(Dict("diameter" => 2, "tolerance" => 1))))
    @debug res
end

end

Script.run()