#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export Bond


mutable struct Bond <: Edge
    """Bond

    * Notation
        * Single bond
            * 0: u - v
            * 1: u ◀ v (Up-arrow)
            * 2: u ▶ v
            * 3: u ◁ v (Down-arrow)
            * 4: u ▷ v
            * 5: u ~ v (Chiral)
        * Double bond
            * 0: v ニ u (clockwise, default)
            * 1: u ニ v (counter-clockwise)
            * 2: u ＝ v (equal length, for terminal bond by default)
            * 3: u × v (Cis-Trans Unknown)
    """
    u::UInt16
    v::UInt16
    order::UInt8

    rotatable::Bool
    aromatic::Bool
    notation::UInt8
    smiles_cis_trans::UInt8
    visible::Bool

    function Bond()
        bond = new()
        bond.order = 1
        bond.rotatable = false
        bond.aromatic = false
        bond.notation = 0
        bond.smiles_cis_trans = 0
        bond.visible = true
        bond
    end
end
