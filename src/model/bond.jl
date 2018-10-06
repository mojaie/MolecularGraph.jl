#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export Bond


mutable struct Bond
    order::Int
    is_lower_first::Bool
    bondtype::Int
    rotatable::Bool
    aromatic::Bool
    smiles_cis_trans::Bool
    visible::Bool

    function Bond()
        initialize!(new())
    end
end


function initialize!(bond::Bond)
    bond.visible = true
    bond
end
