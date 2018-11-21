#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    componentquery!,
    component!,
    fragment!,
    branch!,
    chain!


function componentquery!(state::SmartsParserState)
    """ Start <- Component ('.' Component)*
    """
    if read(state) == '\0'
        return # Empty query
    end
    component!(state)
    while !state.done && read(state) == '.'
        forward!(state)
        state.root = state.node + 1
        component!(state)
    end
end


function component!(state::SmartsParserState)
    """ Component <- '(' Fragment ')' / Fragment
    """
    if read(state) == '('
        # Connectivity restriction
        forward!(state)
        push!(state.mol.connectivity, [state.root])
        fragment!(state)
        c2 = read(state)
        @assert c2 == ')' "unexpected token: $(c2) at $(state.pos)"
        forward!(state)
    else
        fragment!(state)
    end
end


function fragment!(state::AbstractSmartsParser)
    """ Fragment <- Group
    """
    if state isa SmilesParserState && read(state) == '\0'
        return # Empty molecule
    end
    group!(state, nothing)
    # Validity check
    if state.node + 1 == state.root
        throw(MolParseError(
            "unexpected token: $(read(state)) at $(state.pos)"))
    elseif length(state.ringlabel) > 0
        label = collect(keys(state.ringlabel))[1]
        throw(MolParseError("unclosed ring: $(label)"))
    end
end


function group!(state::AbstractSmartsParser, bond)
    """ Group <- Atom ((Bond? Group) / Chain)* Chain
    """
    atomf = state isa SmilesParserState ? atom! : atomquery!
    bondf = state isa SmilesParserState ? bond! : bondquery!
    a = atomf(state)
    if a === nothing
        # Do not start with '(' ex. C((C)C)C is invalid
        throw(MolParseError(
            "unexpected token: $(read(state)) at $(state.pos)"))
    end
    state.node += 1
    updateatom!(state.mol, a, state.node)
    if bond !== nothing
        # Connect branch
        b = connect(bond, state.branch, state.node)
        updatebond!(state.mol, b, bondcount(state.mol) + 1)
    end
    state.branch = state.node
    while true
        if read(state) == '('
            forward!(state)
            buf = state.branch
            b = something(bondf(state), defaultbond(state))
            group!(state, b)
            state.branch = buf
            c = read(state)
            @assert c == ')' "unexpected token: $(c) at $(state.pos)"
            forward!(state)
            if state.done || read(state) in ")."
                # Do not finish with ')' ex, CC(C)(C) should be CC(C)C
                throw(MolParseError(
                    "unexpected token: $(read(state)) at $(state.pos)"))
            end
        else
            chain!(state)
            state.branch = state.node
            if state.done || read(state) in ")."
                break
            end
        end
    end
end


function chain!(state::AbstractSmartsParser)
    """ Chain <- (Bond? (Atom / RingLabel))+
    """
    atomf = state isa SmilesParserState ? atom! : atomquery!
    bondf = state isa SmilesParserState ? bond! : bondquery!
    u = state.branch
    while !state.done
        # Bond?
        if read(state) == '.' && lookahead(state, 1) != '('
            # Disconnected
            forward!(state)
            b = :disconn
        else
            b = bondf(state)
        end

        # RingLabel
        if isdigit(read(state))
            num = parse(Int, read(state))
            forward!(state)
            if num in keys(state.ringlabel)
                (v, rb) = state.ringlabel[num]
                b = something(b, rb, defaultbond(state))
                delete!(state.ringlabel, num) # Ring label is reusable
                b = connect(b, u, v)
                updatebond!(state.mol, b, bondcount(state.mol) + 1)
            else
                state.ringlabel[num] = (u, b)
            end
            continue
        end

        # Atom
        a = atomf(state)
        if a === nothing
            c = read(state)
            @assert c in "().\0" "unexpected token: $(c) at $(state.pos)"
            if b == :disconn
                throw(MolParseError(
                    "unexpected token: $(read(state)) at $(state.pos)"))
            end
            break
        else
            state.node += 1
            updateatom!(state.mol, a, state.node)
        end
        if b == :disconn
            if (state isa SmartsParserState)
                for conn in state.mol.connectivity
                    if conn[1] == state.root
                        push!(conn, state.node)
                    end
                end
            end
        else
            b = something(b, defaultbond(state))
            b = connect(b, u, state.node)
            updatebond!(state.mol, b, bondcount(state.mol) + 1)
        end
        u = state.node
    end
end
