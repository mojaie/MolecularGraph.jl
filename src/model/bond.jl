#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFileBond,
    SmilesBond,
    SmartsBond,
    connect


import ..Graph: similaredge


mutable struct SDFileBond <: Bond
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
    u::Union{Int, Nothing}
    v::Union{Int, Nothing}
    order::Int
    notation::Union{Int, Nothing}
end

SDFileBond(u, v) = SDFileBond(u, v, 1, 0)
SDFileBond(u, v, order) = SDFileBond(u, v, order, 0)
similaredge(b::SDFileBond, u, v) = SDFileBond(u, v, b.order, b.notation)


mutable struct SmilesBond <: Bond
    u::Union{Int, Nothing}
    v::Union{Int, Nothing}
    order::Int
    isaromatic::Union{Bool, Nothing}
    cistrans::Union{Int, Nothing}
end

SmilesBond(u, v) = SmilesBond(u, v, 1, false, nothing)
SmilesBond(u, v, order) = SmilesBond(u, v, order, false, nothing)
similaredge(b::SmilesBond, u, v) = SmilesBond(
    u, v, b.order, b.isaromatic, b.cistrans)


mutable struct SmartsBond <: QueryBond
    u::Union{Int, Nothing}
    v::Union{Int, Nothing}
    query::Pair
end

SmartsBond(u, v) = SmartsBond(u, v, :any=>true)
similaredge(b::SmartsBond, u, v) = SmartsBond(u, v, b.query)
