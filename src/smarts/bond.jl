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
        QueryLiteral(:order, 1),
        QueryOperator(:not, [QueryLiteral(:isaromatic)])
    ]),
    '=' => QueryOperator(:and, [
        QueryLiteral(:order, 2),
        QueryOperator(:not, [QueryLiteral(:isaromatic)])
    ]),
    '#' => QueryOperator(:and, [
        QueryLiteral(:order, 3),
        QueryOperator(:not, [QueryLiteral(:isaromatic)])
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
            QueryLiteral(:order, 1),
            QueryOperator(:not, [QueryLiteral(:isaromatic)])
        ]),
        QueryLiteral(:isaromatic)
    ]))


"""
    bondsymbol!(state::SmartsParserState) -> Union{Pair,Nothing}

BondSymbol <- [-=#@:/\\] / '/?' / '\\?'
"""
function bondsymbol!(state::Union{SMILESParser,SMARTSParser})
    mapping = isa(state, SMILESParser) ? SMILES_BOND_SYMBOL : SMARTS_BOND_SYMBOL
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == '/' && sym2 == '?'
        forward!(state, 2)
        return QueryOperator(:not, [QueryLiteral(:stereo, :down)])
    elseif sym1 == '\\' && sym2 == '?'
        forward!(state, 2)
        return QueryOperator(:not, [QueryLiteral(:stereo, :up)])
    elseif sym1 in keys(mapping)
        forward!(state)
        return mapping[sym1]
    end
    # Implicit single bond returns nothing
end


"""
    bond!(state::SmilesParser) -> Union{SmilesBond,Nothing}

Bond <- BondSymbol?
"""
function bond!(state::SMILESParser{T,V,E}) where {T,V,E}
    q = bondsymbol!(state)
    q === nothing && return
    qd = smiles_dict(q)
    default_qd = Dict(
        :order => 1,
        :isaromatic => false,
        :direction => :unspecified
    )
    merge!(default_qd, qd)
    return E(default_qd)
end


"""
    bond!(state::SmartsParser) -> Union{SmartsBond,Nothing}

Bond <- '~' / (BondSymbol / LogicalCond)?
"""
function bond!(state::SMARTSParser{T,V,E}) where {T,V,E}
    if read(state) == '~'
        forward!(state)
        return E(QueryAny(true))
    end
    q = lglowand!(state, bondsymbol!)
    q === nothing && return
    return E(q)
end
