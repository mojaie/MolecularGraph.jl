#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


function componentquery!(state::SMARTSParser)
    """ Start <- Component ('.' Component)*
    """
    read(state) == '\0' && return  # Empty query: smartstomol("")
    component!(state)
    while !state.done && read(state) == '.'
        forward!(state)
        state.root = state.node + 1
        component!(state)
    end
end


function component!(state::SMARTSParser)
    """ Component <- '(' Fragment ')' / Fragment
    """
    read(state) == '(' || (fragment!(state); return)
    forward!(state)
    push!(state.connectivity, [state.root])  # Set connectivity restriction
    fragment!(state)
    c2 = read(state)
    c2 == ')' || error("component not closed: $(c2) found at $(state.pos)")
    forward!(state)
end


function fragment!(state::Union{SMILESParser,SMARTSParser})
    """ Fragment <- Group
    """
    read(state) == '\0' && return  # Empty molecule
    group!(state, nothing)
    # Validity check
    state.node + 1 == state.root && error("empty group: $(read(state)) found at $(state.pos)")
    length(state.ringlabel) == 0 || error("unclosed ring: $(collect(keys(state.ringlabel))[1])")
end


function group!(state::Union{SMILESParser{T,V,E},SMARTSParser{T,V,E}}, bond) where {T,V,E}
    """ Group <- Atom ((Bond? Group) / Chain)* Chain
    """
    a = atom!(state)
    isempty(a) && error(
        "unexpected token: branch starts with '(' at $(state.pos)")  # ex. CC((CC)C)C  TODO: is still valid?
    state.node += 1
    push!(state.vprops, popfirst!(a))
    if bond !== nothing
        # Connect branch
        push!(state.edges, u_edge(T, state.branch, state.node))
        push!(state.eprops, bond)
    end
    center = state.node
    for h in a
        # hydrogens
        state.node += 1
        push!(state.vprops, h)
        push!(state.edges, u_edge(T, center, state.node))
        push!(state.eprops, defaultbond(state))
    end
    state.branch = center
    while true
        if read(state) == '('
            forward!(state)
            buf = state.branch
            b = something(bond!(state), defaultbond(state))
            group!(state, b)
            state.branch = buf
            c = read(state)
            @assert c == ')' "unexpected token: $(c) at $(state.pos)"
            forward!(state)

            # ex. CC(C)(C) should be CC(C)C but acceptable
            raw"""
            if state.done || read(state) in ")."
                error("unexpected token: branch ends with ')' at $(state.pos)")
            end
            """
        else
            ccen = chain!(state)
            state.branch = ccen
            if state.done || read(state) in ")."
                break
            end
        end
    end
end


function chain!(state::Union{SMILESParser{T,V,E},SMARTSParser{T,V,E}}) where {T,V,E}
    """ Chain <- (Bond? (Atom / RingLabel))+
    """
    u = state.branch
    while !state.done
        # Bond?
        if read(state) == '.' && lookahead(state, 1) != '('  # .( starts a new component
            # Disconnected
            forward!(state)
            b = nothing
        else
            b = something(bond!(state), defaultbond(state))
        end

        # RingLabel
        if isdigit(read(state)) || read(state) == '%'
            if read(state) == '%'
                forward!(state)
                start = state.pos
                while isdigit(lookahead(state, 1))
                    forward!(state)
                end
                num = parse(Int, SubString(state.input, start, state.pos))
                forward!(state)
            else
                num = parse(Int, read(state))
                forward!(state)
            end
            if num in keys(state.ringlabel)
                eidx = state.ringlabel[num]
                v = src(state.edges[eidx])
                rb = state.eprops[eidx]
                db = defaultbond(state)
                if b != db && rb != db && b != rb
                    error("Inconsistent bond props $(v): $(rb), $(u): $(b)")
                end
                delete!(state.ringlabel, num) # Ring label is reusable
                state.edges[eidx] = u_edge(T, u, v)
                state.eprops[eidx] = b == db ? rb : b
            else
                push!(state.edges, Edge{T}(u, u))  # placeholder
                push!(state.eprops, b)
                state.ringlabel[num] = length(state.edges)
            end
            continue
        end

        # Atom
        a = atom!(state)
        if isempty(a)
            c = read(state)
            c in "().\0" || error("unexpected token: $(c) at $(state.pos)")
            isnothing(b) && error("unexpected token: $(read(state)) at $(state.pos)")
            break
        else
            state.node += 1
            push!(state.vprops, popfirst!(a))
        end
        if isnothing(b)
            if isa(state, SMARTSParser)
                for conn in state.connectivity
                    conn[1] == state.root && push!(conn, state.node)
                end
            end
        else
            push!(state.edges, u_edge(T, u, state.node))
            push!(state.eprops, b)
        end
        center = state.node
        for h in a
            # hydrogens
            state.node += 1
            push!(state.edges, u_edge(T, center, state.node))
            push!(state.eprops, defaultbond(state))
            push!(state.vprops, h)
        end
        u = center
    end
    return u
end
