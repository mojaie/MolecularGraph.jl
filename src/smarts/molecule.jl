#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    componentquery!,
    component!,
    connectedquery!,
    fragment!,
    branch!,
    chain!


function componentquery!(state::SmartsParserState)
    """ Start <- Component ('.' Component)*
    """
    if read(state) == '\0'
        return # Empty molecule
    end
    component!(state)
    while !state.done && read(state) == '.'
        forward!(state)
        component!(state)
    end
end


function component!(state::SmartsParserState)
    """ Component <- '(' ConnectedQuery ')' / Fragment
    """
    if read(state) == '('
        forward!(state)
        connectedquery!(state)
        c2 = read(state)
        @assert c2 == ')' "(component!) unexpected token: $(c2)"
        forward!(state)
    else
        fragment!(state, atomcount(state.mol), true)
    end
end


function connectedquery!(state::AbstractSmartsParser)
    """ Connected <- Fragment ('.' Fragment)*
    """
    if state isa SmilesParserState && read(state) == '\0'
        return # Empty molecule
    end
    pos = atomcount(state.mol)
    heads = [pos + 1]
    fragment!(state, pos, true)
    while !state.done && read(state) == '.'
        forward!(state)
        h = atomcount(state.mol)
        fragment!(state, h, true)
        push!(heads, h + 1)
    end
    if state isa SmartsParserState
        push!(state.mol.connectivity, heads)
    end
end


function fragment!(state::AbstractSmartsParser, root, ishead)
    """ Fragment <- Chain Branch?
    """
    chain!(state, root, ishead)
    if !state.done
        if read(state) == '('
            branch!(state)
        end
    end
end


function branch!(state::AbstractSmartsParser)
    """ Branch <- ('(' Fragment ')')+ Fragment
    """
    root = atomcount(state.mol)
    while read(state) == '('
        forward!(state)
        fragment!(state, root, false)
        c = read(state)
        @assert c == ')' "(branch!) unexpected token: $(c)"
        forward!(state)
    end
    fragment!(state, root, false)
end


function chain!(state::AbstractSmartsParser, root, ishead)
    """ Chain <- (Bond? Atom)+
    """
    atomf = state isa SmilesParserState ? atom! : atomquery!
    bondf = state isa SmilesParserState ? bond! : bondquery!
    pos = root
    if ishead
        # first atom
        pos += 1
        a = atomf(state)
        if a === nothing
            throw(MolParseError("unexpected token: $(read(state))"))
        end
        updateatom!(state.mol, a, pos)
    end
    while !state.done
        b = bondf(state)
        a = atomf(state)
        u = pos
        if a === nothing
            c = read(state)
            @assert c in "()." "(chain!) unexpected token: $(c)"
            break
        elseif a isa Pair && a.first == :ring
            # Ring
            if a.second in keys(state.ringlabel)
                (v, rb) = state.ringlabel[a.second]
                b = b === nothing ? rb : b
                delete!(state.ringlabel, a.second) # Ring number is reusable
            else
                state.ringlabel[a.second] = (u, b)
                continue
            end
        else
            v = atomcount(state.mol) + 1
            updateatom!(state.mol, a, v)
            pos = v
        end
        if b === nothing
            if state isa SmilesParserState
                b = SmilesBond(1, false, nothing)
            else
                b = SmartsBond(:BondOrder => 1)
            end
        end
        b = connect(b, u, v)
        updatebond!(state.mol, b, bondcount(state.mol) + 1)
    end
    if pos == root
        # Empty chain is prohibited
        throw(MolParseError("unexpected token: $(read(state))"))
    end
end
