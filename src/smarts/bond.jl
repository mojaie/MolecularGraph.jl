#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


const SMARTS_BOND_SYMBOL = Dict(
    '-' => :bondorder => 1,
    '=' => :bondorder => 2,
    '#' => :bondorder => 3,
    '@' => :isringbond => true,
    ':' => :isaromaticbond => true,
    '/' => :stereo => :up,
    '\\' => :stereo => :down
)


defaultbond(state::SmilesParser) = SmilesBond()
defaultbond(state::SmartsParser) = SmartsBond(:defaultbond => true)


"""
    bondsymbol!(state::SmartsParserState) -> Union{Pair,Nothing}

BondSymbol <- [-=#@:/\\] / '/?' / '\\?'
"""
function bondsymbol!(state::SmartsParserState)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == '/' && sym2 == '?'
        forward!(state, 2)
        return :not => (:stereo => :down)
    elseif sym1 == '\\' && sym2 == '?'
        forward!(state, 2)
        return :not => (:stereo => :up)
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
    elseif fml[1] == :bondorder
        return SmilesBond(fml[2])
    elseif fml[1] == :isaromaticbond
        return SmilesBond(1, true, :unspecified)
    elseif fml[1] == :stereo
        return SmilesBond(1, false, fml[2])
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
        fml = associate_operations(fml)
        return SmartsBond(fml)
    end
    # return nothing: Invalid bond token or implicit single bond
end
