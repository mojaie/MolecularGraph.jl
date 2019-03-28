#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    defaultbond,
    bond!,
    bondsymbol!


defaultbond(state::SmartsParser{SMILES}) = SmilesBond(nothing, nothing)
defaultbond(state::SmartsParser{SMARTS}) = SmartsBond(
    nothing, nothing, :or => (:bondorder => 1, :isaromaticbond => true))


function bond!(state::SmartsParser{SMILES})
    """ Bond <- BondSymbol?
    """
    b = bondsymbol!(state)
    if b === nothing
        return
    elseif b[1] == :bondorder
        return SmilesBond(nothing, nothing, b[2])
    elseif b[1] == :isaromaticbond
        return SmilesBond(nothing, nothing, 1, true, nothing)
    elseif b[1] == :stereo
        return SmilesBond(nothing, nothing, 1, false, b[2])
    end
end


function bond!(state::SmartsParser{SMARTS})
    """ Bond <- '~' / (BondSymbol / LogicalCond)?
    """
    if read(state) == '~'
        forward!(state)
        return SmartsBond(nothing, nothing, :any => true)
    else
        b = lglowand!(state, bondsymbol!)
        if b !== nothing
            return SmartsBond(nothing, nothing, b)
        end
    end
    # return nothing: Invalid bond token or implicit single bond
end


const SMARTS_BOND_SYMBOL = Dict(
    '-' => :bondorder => 1,
    '=' => :bondorder => 2,
    '#' => :bondorder => 3,
    '@' => :bond_isringmem => true,
    ':' => :isaromaticbond => true,
    '/' => :stereo => 1,
    '\\' => :stereo => 2
)


function bondsymbol!(state::SmartsParser)
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
