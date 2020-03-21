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
defaultbond(state::SmartsParser) = SmartsBond(
    :or => (:bondorder => 1, :isaromaticbond => true))


function bondsymbol!(state::SmartsParserState)
    """ BondSymbol <- [-=#@:/\\] / '/?' / '\\?'
    """
    c = read(state)
    if c in keys(SMARTS_BOND_SYMBOL)
        cond = SMARTS_BOND_SYMBOL[c]
        forward!(state)
        if c == '/' && read(state) == '?'
            forward!(state)
            return :not => (:stereo => :down)
        elseif c == '\\' && read(state) == '?'
            forward!(state)
            return :not => (:stereo => :up)
        else
            return cond
        end
    end
    # Implicit single bond returns nothing
    return
end


function bond!(state::SmilesParser)
    """ Bond <- BondSymbol?
    """
    b = bondsymbol!(state)
    if b === nothing
        return
    elseif b[1] == :bondorder
        return SmilesBond(b[2])
    elseif b[1] == :isaromaticbond
        return SmilesBond(1, true, :unspecified)
    elseif b[1] == :stereo
        return SmilesBond(1, false, b[2])
    end
    return
end


function bond!(state::SmartsParser)
    """ Bond <- '~' / (BondSymbol / LogicalCond)?
    """
    if read(state) == '~'
        forward!(state)
        return SmartsBond()
    else
        b = lglowand!(state, bondsymbol!)
        if b !== nothing
            return SmartsBond(b)
        end
    end
    # return nothing: Invalid bond token or implicit single bond
    return
end
