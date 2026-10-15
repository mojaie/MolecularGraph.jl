#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


const SMARTS_BOND_ORDER = Dict(
    '-' => "1",
    '/' => "1",
    '\\' => "1",
    '=' => "2",
    '#' => "3"
)


"""
    bondsymbol!(state::SMILESParser) -> Union{NamedTuple,Nothing}

BondSymbol <- [-=#:/\\]

Return all bond properties as a NamedTuple of a fixed type (unspecified fields are defaults)
so that the bond constructor can be statically dispatched.
"""
function bondsymbol!(state::SMILESParser)
    sym1 = read(state)
    order, isaromatic, direction = 1, false, :unspecified
    if sym1 == '-'
        order = 1
    elseif sym1 == '='
        order = 2
    elseif sym1 == '#'
        order = 3
    elseif sym1 == ':'
        isaromatic = true
    elseif sym1 == '/'
        direction = :up
    elseif sym1 == '\\'
        direction = :down
    else
        return  # end token found (or implicit bond), no position move
    end
    forward!(state)
    return (order=order, isaromatic=isaromatic, direction=direction)
end


"""
    bondsymbol!(state::SMARTSParser, qtree::QueryTree{T,V}) where {T,V} -> T

BondSymbol <- [-=#@:/\\] / '/?' / '\\?'
"""
function bondsymbol!(
        state::SMARTSParser, qtree::QueryTree{T,V}) where {T<:Integer,V<:QueryNode}
    sym = read(state)
    sym in "-=#@:/\\" || return zero(T)  # end token found (or implicit bond), no position move
    forward!(state)
    sym == ':' && return add_qnode!(qtree, qtrue(:isaromatic))
    sym == '@' && return add_qnode!(qtree, qtrue(:is_in_ring))
    node = add_qnode!(qtree, qand())
    add_qnode!(qtree, node, qeq(:order, SMARTS_BOND_ORDER[sym]))
    c = add_qnode!(qtree, node, qnot())
    add_qnode!(qtree, c, qtrue(:isaromatic))
    sym in "/\\" || return node
    if read(state) == '?'
        forward!(state)
        add_qnode!(qtree, node, qeq(:direction, sym == '/' ? "upuns" : "downuns"))
    else
        add_qnode!(qtree, node, qeq(:direction, sym == '/' ? "up" : "down"))
    end
    return node
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
    # aromatic=true OR (order=1 AND aromatic=false)
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

Bond <- '~' / (BondSymbol / LogicalCond)*
"""
function bond!(state::SMARTSParser{T,V,E}) where {T,V,E}
    if read(state) == '~'
        forward!(state)
        return E(Tuple{Int,Int}[], [qanytrue()])
    end
    qtree = E()
    vs = Int[]
    v = lglowand!(state, qtree)
    v == 0 && return
    while v != 0
        push!(vs, v)
        v = lglowand!(state, qtree)
    end
    length(vs) == 1 && return qtree
    node = add_qnode!(qtree, qand())
    for i in vs
        add_qedge!(qtree, node, i)
    end
    return qtree
end
