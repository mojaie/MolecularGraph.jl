#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    default_annotation!


function default_annotation!(mol::VectorMol)
    molgraph_topology!(mol)
    default_annot!(mol)
    atom_annot!(mol)
    group_annot!(mol)
    rotatable!(mol)
    aromatic!(mol)
    return
end
