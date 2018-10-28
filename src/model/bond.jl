#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Bond,
    sdfbond,
    smilesbond


struct Bond <: AbstractEdge
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
    u::Int
    v::Int
    order::Int
    sdf_notation::Union{Int, Nothing}
    smiles_aromatic::Union{Bool, Nothing}
    smiles_cistrans::Union{Int, Nothing}
end

Bond(b::Bond, u::Integer, v::Integer) = Bond(
    u, v, b.order,
    b.sdf_notation, b.smiles_aromatic, b.smiles_cistrans
)


function sdfbond(u, v, order, notation)
    Bond(u, v, order, notation, nothing, nothing)
end


function smilesbond(u, v, order, aromatic, cistrans)
    Bond(u, v, order, nothing, aromatic, cistrans)
end
