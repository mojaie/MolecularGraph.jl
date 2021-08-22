#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


const SMARTS_BOND_SYMBOL = Dict(
    '-' => QueryFormula(:bondorder, 1),
    '=' => QueryFormula(:bondorder, 2),
    '#' => QueryFormula(:bondorder, 3),
    '@' => QueryFormula(:isringbond, true),
    ':' => QueryFormula(:isaromaticbond, true),
    '/' => QueryFormula(:stereo, :up),
    '\\' => QueryFormula(:stereo, :down)
)


defaultbond(state::SmilesParser) = SmilesBond()
defaultbond(state::SmartsParser) = SmartsBond(QueryFormula(:defaultbond, true))


"""
    bondsymbol!(state::SmartsParserState) -> Union{Pair,Nothing}

BondSymbol <- [-=#@:/\\] / '/?' / '\\?'
"""
function bondsymbol!(state::SmartsParserState)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == '/' && sym2 == '?'
        forward!(state, 2)
        return QueryFormula(:not, QueryFormula(:stereo, :down))
    elseif sym1 == '\\' && sym2 == '?'
        forward!(state, 2)
        return QueryFormula(:not, QueryFormula(:stereo, :up))
    elseif sym1 in keys(SMARTS_BOND_SYMBOL)
        forward!(state)
        return SMARTS_BOND_SYMBOL[sym1]
    end
    # Implicit single bond returns nothing
end



"""
    bond!(state::SmilesParser) -> Union{SmilesBond,Nothing}

Bond <- BondSymbol?
"""
function bond!(state::SmilesParser)
    fml = bondsymbol!(state)
    if fml === nothing
        return
    elseif fml.key == :bondorder
        return SmilesBond(fml.value)
    elseif fml.key == :isaromaticbond
        return SmilesBond(1, true, :unspecified)
    elseif fml.key == :stereo
        return SmilesBond(1, false, fml.value)
    end
    # return nothing
end


"""
    bond!(state::SmartsParser) -> Union{SmartsBond,Nothing}

Bond <- '~' / (BondSymbol / LogicalCond)?
"""
function bond!(state::SmartsParser)
    if read(state) == '~'
        forward!(state)
        return SmartsBond()
    end
    fml = lglowand!(state, bondsymbol!)
    if fml !== nothing
        fml = tidyformula(fml)
        return SmartsBond(fml)
    end
    # return nothing: Invalid bond token or implicit single bond
end
