using MolecularGraph


function getmw(smiles::String)
    mol = smilestomol(smiles)
    precalculate!(mol)
    return standardweight(Float64, mol)
end


function getstruct(smiles::String)
    mol = smilestomol(smiles)
    precalculate!(mol)
    return drawsvg(mol, 200, 200)
end


teststr = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
getmw(teststr)
getstruct(teststr)