#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const SMILES_BONDS = Dict{}(
    '-' => (order=1,),
    '=' => (order=2,),
    '#' => (order=3,),
    ':' => (isaromatic=true,),
    '/' => (direction=:up,),
    '\\' => (direction=:down,)
)

const SMARTS_NON_AROMATIC_BONDS = Dict('-' => "1", '=' => "2", '#' => "3")

const SMARTS_OTHER_BONDS = Dict(
    '@' => qtrue(:is_in_ring),
    ':' => qtrue(:isaromatic),
    '/' => qeq(:direction, "up"),
    '\\' => qeq(:direction, "down")
)


"""
    bondsymbol!(state::SMILESParser) -> Union{NamedTuple,Nothing}

BondSymbol <- [-=#:/\\]
"""
function bondsymbol!(state::SMILESParser)
    sym1 = read(state)
    if sym1 in keys(SMILES_BONDS)
        forward!(state)
        return SMILES_BONDS[sym1]
    end
    return  # end token found (or implicit bond), no position move
end


"""
    bondsymbol!(state::SMARTSParser, qtree::QueryTree{T,V}) where {T,V} -> T

BondSymbol <- [-=#@:/\\] / '/?' / '\\?'
"""
function bondsymbol!(
        state::SMARTSParser, qtree::QueryTree{T,V}) where {T<:Integer,V<:QueryNode}
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == '/' && sym2 == '?'
        forward!(state, 2)
        node = add_qnode!(qtree, qnot())
        add_qnode!(qtree, node, qeq(:stereo, "down"))
        return node
    elseif sym1 == '\\' && sym2 == '?'
        forward!(state, 2)
        node = add_qnode!(qtree, qnot())
        add_qnode!(qtree, node, qeq(:stereo, "up"))
        return node
    elseif sym1 in keys(SMARTS_NON_AROMATIC_BONDS)
        forward!(state)
        node = add_qnode!(qtree, qand())
        add_qnode!(qtree, node, qeq(:order, SMARTS_NON_AROMATIC_BONDS[sym1]))
        c = add_qnode!(qtree, node, qnot())
        add_qnode!(qtree, c, qtrue(:isaromatic))
        return node
    elseif sym1 in keys(SMARTS_OTHER_BONDS)
        forward!(state)
        return add_qnode!(qtree, SMARTS_OTHER_BONDS[sym1])
    end
    return zero(T)  # end token found (or implicit bond), no position move
end


"""
    defaultbond(state::SMILESParser) -> SMILESBond

SMILES default bond (single or aromatic bond)
"""
function defaultbond(state::SMILESParser{T,V,E}) where {T,V,E}
    return E()
end


"""
    defaultbond(state::SMARTSParser) -> QueryBond

SMARTS default bond (single and non-aromatic, or aromatic)
"""
function defaultbond(state::SMARTSParser{T,V,E}) where {T,V,E}
    return E(
        [(1, 2), (2, 3), (2, 4), (4, 5), (1, 6)],
        [qor(), qand(), qeq(:order, "1"), qnot(), qtrue(:isaromatic), qtrue(:isaromatic)]
    )
end


"""
    bond!(state::SMILESParser) -> Union{QueryTree,Nothing}

Bond <- BondSymbol?
"""
function bond!(state::SMILESParser{T,V,E}) where {T,V,E}
    v = bondsymbol!(state)
    isnothing(v) && return
    return E(;v...)
end


"""
    bond!(state::SMARTSParser) -> Union{QueryTree,Nothing}

Bond <- '~' / (BondSymbol / LogicalCond)?
"""
function bond!(state::SMARTSParser{T,V,E}) where {T,V,E}
    if read(state) == '~'
        forward!(state)
        return E(Tuple{Int,Int}[], [qanytrue()])
    end
    qtree = E()
    v = lglowand!(state, qtree)
    v == 0 && return
    return qtree
end
