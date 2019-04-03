#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFileBond, SmilesBond, SmartsBond,
    setorder


struct SDFileBond <: Bond
    """Bond
    * Notation
        * Single bond
            * 0: u - v
            * 1: u ◀ v (Up-arrow)
            * 4: u ~ v (Up or down)
            * 6: u ◁ v (Down-arrow)
        * Double bond
            * 0: v ニ u (clockwise, default)
            * 1: u ニ v (counter-clockwise)
            * 2: u ＝ v (equal length, for terminal bond by default)
            * 3: u × v (Cis-Trans Unknown)
    """
    order::Int
    notation::Union{Int,Nothing}
end

SDFileBond() = SDFileBond(1, 0)
SDFileBond(order) = SDFileBond(order, 0)

setorder(edge::SDFileBond, order) = SDFileBond(order, edge.notation)


struct SmilesBond <: Bond
    order::Int
    isaromatic::Union{Bool, Nothing}
    cistrans::Union{Int, Nothing}
end

SmilesBond() = SmilesBond(1, false, nothing)
SmilesBond(order) = SmilesBond(order, false, nothing)

setorder(edge::SmilesBond, order
    ) = SmilesBond(order, edge.isaromatic, edge.cistrans)


struct SmartsBond <: QueryBond
    query::Pair
end

SmartsBond() = SmartsBond(:any => true)
