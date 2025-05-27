#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    moltosmiles, moltosmarts, to_rdkmol, to_rdkqmol

using RDKitMinimalLib: get_mol, get_qmol, get_smiles


function to_rdkmol(mol::MolGraph)
    return get_mol(to_json(Val{:rdkit}(), mol))
end


function to_rdkqmol(mol::MolGraph{Int,QueryTree,QueryTree})
    return get_qmol(to_json(Val{:rdkit}(), mol))
end


function moltosmiles(mol::MolGraph, details=nothing)
    return get_smiles(to_rdkmol(mol), details)
end


function moltosmarts(mol::MolGraph{Int,QueryTree,QueryTree}, details=nothing)
    return get_smiles(to_rdkmol(mol), details)
end
