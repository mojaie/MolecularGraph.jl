#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


export
    molecular_weight,
    H_acceptor_count,
    H_donor_count,
    rotatable_count


function molecular_weight(mol::VectorMol)
    reduce(+, mol.v[:MolWeight]; init=0)
end

function H_acceptor_count(mol::VectorMol)
    reduce(+, mol.v[:H_Acceptor]; init=0)
end

function H_donor_count(mol::VectorMol)
    reduce(+, mol.v[:H_Donor]; init=0)
end

function rotatable_count(mol::VectorMol)
    reduce(+, mol.v[:Rotatable]; init=0)
end
