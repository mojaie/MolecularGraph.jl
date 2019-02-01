#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    default_annotation!


function default_annotation!(mol::VectorMol; recalculate=false)
    if recalculate
        empty!(mol.v)
        empty!(mol.annotation)
    end
    topology!(mol)
    elemental!(mol)
    rotatable!(mol)
    aromatic!(mol)
    return
end
