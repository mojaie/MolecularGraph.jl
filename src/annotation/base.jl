#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    default_annotation!,
    clear_annotation!


function default_annotation!(mol::VectorMol)
    topology!(mol)
    elemental!(mol)
    rotatable!(mol)
    aromatic!(mol)
    return
end

function clear_annotation!(mol::VectorMol)
    empty!(mol.vector)
    empty!(mol.annotation)
    return
end
