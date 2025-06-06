#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const SMILES_BOND_SYMBOL = Dict(
    '-' => QueryLiteral(:order, 1),
    '=' => QueryLiteral(:order, 2),
    '#' => QueryLiteral(:order, 3),
    '@' => QueryLiteral(:is_in_ring),
    ':' => QueryLiteral(:isaromatic),
    '/' => QueryLiteral(:direction, :up),
    '\\' => QueryLiteral(:direction, :down)
)


const SMARTS_BOND_SYMBOL = Dict(
    '-' => QueryOperator(:and, [
        QueryTree(QueryLiteral(:order, 1)),
        QueryTree(QueryOperator(:not, [QueryTree(QueryLiteral(:isaromatic))]))
    ]),
    '=' => QueryOperator(:and, [
        QueryTree(QueryLiteral(:order, 2)),
        QueryTree(QueryOperator(:not, [QueryTree(QueryLiteral(:isaromatic))]))
    ]),
    '#' => QueryOperator(:and, [
        QueryTree(QueryLiteral(:order, 3)),
        QueryTree(QueryOperator(:not, [QueryTree(QueryLiteral(:isaromatic))]))
    ]),
    '@' => QueryLiteral(:is_in_ring),
    ':' => QueryLiteral(:isaromatic),
    '/' => QueryLiteral(:direction, :up),
    '\\' => QueryLiteral(:direction, :down)
)


defaultbond(state::SMILESParser{T,V,E}) where {T,V,E} = E(
    Dict(
        :order => 1,
        :isaromatic => false,
        :direction => :unspecified
    ))
defaultbond(state::SMARTSParser{T,V,E}
    ) where {T,V,E} = E(QueryOperator(:or, [
        QueryOperator(:and, [
            QueryTree(QueryLiteral(:order, 1)),
            QueryTree(QueryOperator(:not, [QueryLiteral(:isaromatic)]))
        ]),
        QueryLiteral(:isaromatic)
    ]))


"""
    bondsymbol!(state::SmartsParserState) -> Union{Pair,Nothing}

BondSymbol <- [-=#@:/\\] / '/?' / '\\?'
"""
function bondsymbol!(state::T) where T <: AbstractSMARTSParser
    mapping = isa(state, SMILESParser) ? SMILES_BOND_SYMBOL : SMARTS_BOND_SYMBOL
    sym1 = readtoken(state)
    sym2 = lookahead(state, 1)
    if sym1 == '/' && sym2 == '?'
        forward!(state, 2)
        return QueryOperator(:not, [QueryTree(QueryLiteral(:stereo, :down))])
    elseif sym1 == '\\' && sym2 == '?'
        forward!(state, 2)
        return QueryOperator(:not, [QueryTree(QueryLiteral(:stereo, :up))])
    elseif sym1 in keys(mapping)
        forward!(state)
        return mapping[sym1]
    end
    return EndToken()  # Implicit bond
end


"""
    bond!(state::SmilesParser) -> Union{SmilesBond,Nothing}

Bond <- BondSymbol?
"""
function bond!(state::SMILESParser{T,V,E}) where {T,V,E}
    q = bondsymbol!(state)
    q isa EndToken && return
    qd = SMILESBondContainer()
    smiles_props!(qd, q)
    return E(qd)
end


"""
    bond!(state::SmartsParser) -> Union{SmartsBond,Nothing}

Bond <- '~' / (BondSymbol / LogicalCond)?
"""
function bond!(state::SMARTSParser{T,V,E}) where {T,V,E}
    if readtoken(state) == '~'
        forward!(state)
        return E(QueryAny(true))
    end
    q = lglowand!(state, bondsymbol!)
    return q isa EndToken ? nothing : E(q)
end
