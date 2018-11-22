#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    defaultbond,
    bond!,
    bondsymbol!


defaultbond(state::SmilesParserState) = SmilesBond(1, false, nothing)
defaultbond(state::SmartsParserState) = SmartsBond(:BondOrder => 1)


function bond!(state::SmilesParserState)
    """ Bond <- BondSymbol?
    """
    b = bondsymbol!(state)
    if b === nothing
        return
    elseif b[1] == :BondOrder
        return SmilesBond(b[2], false, nothing)
    elseif b[1] == :Aromatic
        return SmilesBond(1, true, nothing)
    elseif b[1] == :stereo
        return SmilesBond(1, false, b[2])
    end
end


function bond!(state::SmartsParserState)
    """ Bond <- '~' / (BondSymbol / LogicalCond)?
    """
    if read(state) == '~'
        forward!(state)
        return SmartsBond(:any => true)
    else
        b = lglowand!(state, bondsymbol!)
        if b !== nothing
            return SmartsBond(b)
        end
    end
    # return nothing: Invalid bond token or implicit single bond
end


const SMARTS_BOND_SYMBOL = Dict(
    '-' => :BondOrder => 1,
    '=' => :BondOrder => 2,
    '#' => :BondOrder => 3,
    '@' => :RingBond => true,
    ':' => :Aromatic => true,
    '/' => :stereo => 1,
    '\\' => :stereo => 2
)


function bondsymbol!(state::AbstractSmartsParser)
    """ BondSymbol <- [-=#@:/\\] / '/?' / '\\?'
    """
    c = read(state)
    if c in keys(SMARTS_BOND_SYMBOL)
        cond = SMARTS_BOND_SYMBOL[c]
        forward!(state)
        if c in ('/', '\\') && read(state) == '?'
            forward!(state)
            # '/?' => 3, '\?' => 4
            return :stereo => cond.second + 2
        else
            return cond
        end
    end
    # Implicit single bond returns nothing
end
