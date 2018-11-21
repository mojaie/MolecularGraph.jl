#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    bond!,
    bondquery!,
    bondsymbol!


function bond!(state::SmilesParserState)
    """ Bond <- [-=#:/\\]?
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
    # raise Error
end


function bondquery!(state::SmartsParserState)
    """ BondQuery <- '~' / BondLogicalCond '?'?
    """
    if read(state) == '~'
        forward!(state)
        return SmartsBond(:any => true)
    else
        b = lglowand!(state, bondsymbol!)
        if b === nothing
            return
        else
            return SmartsBond(b)
        end
    end
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
