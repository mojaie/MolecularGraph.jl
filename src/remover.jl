#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    remove_H,
    remove_all_H,
    remove_water,
    remove_salt



function remove_H(mol::MutableMolecule)
    mol = copy(mol)
    # Implicit hydrogens
    implicitHs = []
end

function remove_all_H(mol::MutableMolecule)
    # Implicit hydrogens
    implicitHs = []
    # Hydrogens connected to chiral center
    chiralHs = []
end

function remove_water(mol::MutableMolecule)

end
