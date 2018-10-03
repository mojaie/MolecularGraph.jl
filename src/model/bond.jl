

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
