#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    moltosmiles, moltosmarts

using RDKitMinimalLib


function moltosmiles(mol::MolGraph, details=nothing)
    return get_smiles(to_rdkmol(mol), details)


function moltosmarts(mol::QueryMol, details=nothing)
    return get_smiles(to_rdkmol(mol), details)
end
