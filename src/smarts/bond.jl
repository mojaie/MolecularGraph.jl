#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const SMILES_BOND_SYMBOL = Dict(
    '-' => (v -> v[1], [(:order, 1)]),
    '=' => (v -> v[1], [(:order, 2)]),
    '#' => (v -> v[1], [(:order, 3)]),
    '@' => (v -> v[1], [(:isring,)]),
    ':' => (v -> v[1], [(:isaromatic,)]),
    '/' => (v -> v[1], [(:stereo, :up)]),
    '\\' => (v -> v[1], [(:stereo, :down)])
)


const SMARTS_BOND_SYMBOL = Dict(
    '-' => (v -> v[1] & ~v[2], [(:order, 1), (:isaromatic,)]),
    '=' => (v -> v[1] & ~v[2], [(:order, 2), (:isaromatic,)]),
    '#' => (v -> v[1] & ~v[2], [(:order, 3), (:isaromatic,)]),
    '@' => (v -> v[1], [(:isring,)]),
    ':' => (v -> v[1], [(:isaromatic,)]),
    '/' => (v -> v[1], [(:stereo, :up)]),
    '\\' => (v -> v[1], [(:stereo, :down)])
)

defaultbond(state::SMILESParser{T,V,E}) where {T,V,E} = E()
defaultbond(state::SMARTSParser{T,V,E}
    ) where {T,V,E} = E(v -> v[1] & ~v[2] | v[2], [(:order, 1), (:isaromatic,)])


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
        return (v -> ~v[1], [(:stereo, :down)])
    elseif sym1 == '\\' && sym2 == '?'
        forward!(state, 2)
        return (v -> ~v[1], [(:stereo, :up)])
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
    qd = Dict{Symbol,Any}()
    for fml in q[2]  # operator is fixed to :eq so far
        qd[fml[1]] = length(fml) == 1 ? q[1](trues(length(q[2]))) : fml[2]
    end
    return E(qd)
end


"""
    bond!(state::SmartsParser) -> Union{SmartsBond,Nothing}

Bond <- '~' / (BondSymbol / LogicalCond)?
"""
function bond!(state::SMARTSParser{T,V,E}) where {T,V,E}
    if read(state) == '~'
        forward!(state)
        return E(any_query(true)...)
    end
    q = lglowand!(state, bondsymbol!)
    q === nothing && return
    return E(q...)
end
