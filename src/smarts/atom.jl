#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# TODO: radical


const SMARTS_ATOM_COND_SYMBOL = Dict(
    'X' => :connectivity,
    'D' => :degree,
    'v' => :valence,
    'H' => :total_hydrogens,
    'r' => :smallest_ring,
    'R' => :ring_count
)

const SMARTS_CHARGE_SIGN = Dict(
    '+' => 1,
    '-' => -1
)


"""
    atomsymbol!(state::SMILESParser{T,V,E}) -> Vector{V}

Atomsymbol <- Br / Cl / [BCcNnOoPpSsFI]
"""
function atomsymbol!(state::SMILESParser{T,V,E}) where {T,V,E}
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 * sym2 in ("Br", "Cl")
        forward!(state, 2)
        return [V(symbol=Symbol(sym1 * sym2))]
    elseif sym1 in "BCNOPSFI"
        forward!(state)
        return [V(symbol=Symbol(sym1), isaromatic=false)]
    elseif sym1 in "bcnops"
        forward!(state)
        return [V(symbol=Symbol(uppercase(sym1)), isaromatic=true)]
    elseif sym1 in ('(', ')', '.', '\0')
        # end token found, no position move
        return V[]
    end
    error("unexpected token: '$(sym1)' at $(state.pos)")
end


"""
    atomsymbol!(state::SMARTSParser, qtree::QueryTree) -> Integer

Atomsymbol <- Br / Cl / [AaBCcNnOoPpSsFI*]
"""
function atomsymbol!(state::SMARTSParser, qtree::QueryTree{T,V}) where {T,V}
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 * sym2 in ("Br", "Cl")
        forward!(state, 2)
        return add_qnode!(qtree, qeq(:symbol, sym1 * sym2))
    elseif sym1 in "BCNOPSFI"
        forward!(state)
        node = add_qnode!(qtree, qand())
        c1 = add_qnode!(qtree, node, qeq(:symbol, string(sym1)))
        c2 = add_qnode!(qtree, node, qnot())
        c2c = add_qnode!(qtree, c2, qtrue(:isaromatic))
        return node
    elseif sym1 in "bcnops"
        forward!(state)
        node = add_qnode!(qtree, qand())
        c1 = add_qnode!(qtree, node, qeq(:symbol, uppercase(string(sym1))))
        c2 = add_qnode!(qtree, node, qtrue(:isaromatic))
        return node
    elseif sym1 == 'A'
        forward!(state)
        node = add_qnode!(qtree, qnot())
        c = add_qnode!(qtree, node, qtrue(:isaromatic))
        return node
    elseif sym1 == 'a'
        forward!(state)
        return add_qnode!(qtree, qtrue(:isaromatic))
    elseif sym1 == '*'
        forward!(state)
        return add_qnode!(qtree, qanytrue())
    elseif sym1 in ('(', ')', '.', '\0')
        # end token found, no position move
        return zero(T)
    end
    error("unexpected token: '$(sym1)' at $(state.pos)")
end


function atompropsymbol!(
        state::SMILESParser, sym::Symbol, slow::Bool, cond::NTuple)
    if sym in cond
        return (symbol=sym, isaromatic=slow)
    end
    slow && error("aromatic $(sym) is not supported")
    return (symbol=sym,)
end


function atompropsymbol!(
        state::SMARTSParser, qtree::QueryTree, sym::Symbol,
        slow::Bool, cond::NTuple)
    if sym in cond
        node = add_qnode!(qtree, qand())
        add_qnode!(qtree, node, qeq(:symbol, string(sym)))
        if slow
            add_qnode!(qtree, node, qtrue(:isaromatic))
        else
            c = add_qnode!(qtree, node, qnot())
            add_qnode!(qtree, c, qtrue(:isaromatic))
        end
        return node
    end
    slow && error("aromatic $(sym) is not supported")
    return add_qnode!(qtree, qeq(:symbol, string(sym)))
end


function atompropcharge!(state::AbstractSMARTSParser, sym1::Char, sym2::Char)
    forward!(state)
    if isdigit(sym2)
        chg = parse(Int, sym2)
        forward!(state)
    else
        chg = 1
        while read(state) == sym1
            forward!(state)
            chg += 1
        end
    end
    return chg
end


function atompropnumber!(state::AbstractSMARTSParser)
    start = state.pos
    while isdigit(lookahead(state, 1))
        forward!(state)
    end
    num = parse(Int, SubString(state.input, start, state.pos))
    forward!(state)
    return num
end


"""
    atompropcond!(state::SMILESParser) -> Union{NamedTuple,Nothing}

AtomProp <- Stereo / CHG / H
"""
function atompropcond!(state::SMILESParser)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == 'H' # Hydrogen count
        forward!(state)
        if isdigit(sym2)
            num = parse(Int, sym2)
            forward!(state)
        else
            num = 1
        end
        return (total_hydrogens=num,)
    elseif sym1 in keys(SMARTS_CHARGE_SIGN)  # Charge
        chg = atompropcharge!(state, sym1, sym2)
        return (charge=chg * SMARTS_CHARGE_SIGN[sym1],)
    elseif sym1 == '@'
        # Stereo: @ -> anticlockwise, @@ -> clockwise
        forward!(state)
        cw = :anticlockwise
        if read(state) == '@'
            forward!(state)
            cw = :clockwise
        end
        return (stereo=cw,)
    end
    return  # can be nothing
end


"""
    atompropcond!(state::SMARTSParser, qtree::QueryTree) -> Integer

AtomProp <- '\$(' RecursiveQuery ')' / Stereo / CHG / [DHRrvX]
"""
function atompropcond!(state::SMARTSParser, qtree::QueryTree)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if haskey(SMARTS_ATOM_COND_SYMBOL, sym1)
        # Neighbor and ring conditions
        forward!(state)
        cond = SMARTS_ATOM_COND_SYMBOL[sym1]
        if isdigit(sym2)
            num = parse(Int, sym2)
            forward!(state)
        elseif sym1 in ('r', 'R')
            node = add_qnode!(qtree, qnot())
            add_qnode!(qtree, node, qeq(:ring_count, "0"))
            return node
        else
            num = 1
        end
        return add_qnode!(qtree, qeq(cond, string(num)))
    elseif sym1 in keys(SMARTS_CHARGE_SIGN)  # Charge
        chg = atompropcharge!(state, sym1, sym2)
        return add_qnode!(qtree, qeq(:charge, string(chg * SMARTS_CHARGE_SIGN[sym1])))
    elseif sym1 == '@'
        # Stereo: @ -> anticlockwise, @@ -> clockwise, ? -> or not specified
        forward!(state)
        cw = false
        if read(state) == '@'
            forward!(state)
            cw = true
        end
        if read(state) == '?'
            forward!(state)
            node = add_qnode!(qtree, qnot())
            add_qnode!(qtree, node, qeq(:stereo, cw ? "anticlockwise" : "clockwise"))
            return node
        else
            return add_qnode!(qtree, qeq(:stereo, cw ? "clockwise" : "anticlockwise"))
        end
    elseif sym1 == '$' && sym2 == '('  # Recursive
        forward!(state, 2)
        start = state.pos
        toclose = 1
        while true
            if read(state) == ')'
                if toclose == 1
                    break
                else
                    toclose -= 1
                end
            elseif read(state) == '('
                toclose += 1
            end
            forward!(state)
        end
        q = SubString(state.input, start, state.pos - 1)
        forward!(state)
        return add_qnode!(qtree, qeq(:recursive, string(q)))
    end
    return 0  # End token found
end


"""
    atompropsymcond!(state::SMILESParser) -> Union{NamedTuple,Nothing}

AtomPropSymConds <- Symbol / AtomNum
"""
function atompropsymcond!(state::SMILESParser)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == '\0'  # End token
        return 
    elseif sym2 != '\0' && haskey(ATOMSYMBOLMAP, Symbol(uppercase(sym1), sym2))
        # Two letter atoms
        forward!(state, 2)
        return atompropsymbol!(
            state, Symbol(uppercase(sym1), sym2),
            islowercase(sym1), (:As, :Se)
        )
    elseif haskey(ATOMSYMBOLMAP, Symbol(uppercase(sym1))) && !isdigit(sym2)
        # Single letter atoms
        forward!(state)
        return atompropsymbol!(
            state, Symbol(uppercase(sym1)),
            islowercase(sym1), (:B, :C, :N, :O, :P, :S)
        )
    end
    return  # can be nothing
end


"""
    atompropsymcond!(state::SMARTSParser, qtree::QueryTree) -> Integer

AtomPropSymConds <- Symbol / AtomNum
"""
function atompropsymcond!(state::SMARTSParser, qtree::QueryTree)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    vs = Int[]
    if sym1 == '\0'  # End token
        return 0
    elseif sym1 == '#'  # Atomic number
        forward!(state)
        num = atompropnumber!(state)
        return add_qnode!(qtree, qeq(:symbol, string(atom_symbol(num))))
    elseif sym2 != '\0' && haskey(ATOMSYMBOLMAP, Symbol(uppercase(sym1), sym2))
        # Two letter atoms
        forward!(state, 2)
         return atompropsymbol!(
            state, qtree, Symbol(uppercase(sym1), sym2),
            islowercase(sym1), (:As, :Se)
        )
    elseif sym1 == 'A'
        forward!(state)
        node = add_qnode!(qtree, qnot())
        add_qnode!(qtree, node, qtrue(:isaromatic))
        return node
    elseif sym1 == 'a'
        forward!(state)
        return add_qnode!(qtree, qtrue(:isaromatic))
    elseif sym1 == '*'
        forward!(state)
        return add_qnode!(qtree, qanytrue())
    elseif sym1 in ('D', 'T')  # D is not a deuterium, but Degree
        return 0
    elseif haskey(ATOMSYMBOLMAP, Symbol(uppercase(sym1))) && !isdigit(sym2)  # `H2` is not a symbol, but HydrogenCount
        # Single letter atoms
        forward!(state)
        return atompropsymbol!(
            state, qtree, Symbol(uppercase(sym1)),
            islowercase(sym1), (:B, :C, :N, :O, :P, :S)
        )
    end
    return 0
end


"""
    atomprop!(state::SMILESParser) -> Union{NamedTuple,Nothing}

AtomProp <- ((Mass? / AtomPropSymConds) / AtomPropConds)+
"""
function atomprop!(state::SMILESParser{T,V,E}) where {T,V,E}
    sym1 = read(state)
    vs = NamedTuple[]
    hcnt = 0
    hasiso = isdigit(sym1)
    if hasiso  # Isotope
        num = atompropnumber!(state)
        push!(vs, (isotope=num,))
    end
    v = atompropsymcond!(state)
    if !isnothing(v)
        push!(vs, v)
    end
    hasiso && isnothing(v) && error("(atomprop!) atom symbol not specified")
    v = atompropcond!(state)
    while !isnothing(v)
        if haskey(v, :total_hydrogens)
            hcnt = v[:total_hydrogens]
        else
            push!(vs, v)
        end
        v = atompropcond!(state)
    end
    isempty(vs) && hcnt == 0 && error("(atomprop!) empty atomprop")
    # explicit hydrogen (e.g. [CH3]) -> hydrogen nodes
    merged = isempty(vs) ? NamedTuple() : merge(vs...)
    if !haskey(merged, :symbol)  # special case: [H2]
        merged = merge(merged, (symbol=:H,))
        hcnt -= 1
    end
    return [V(;merged...), (V(;symbol=:H) for _ in 1:hcnt)...]
end


"""
    atomprop!(state::SMARTSParser, qtree::QueryTree) -> Integer

AtomProp <- ((Mass? / AtomPropSymConds) / AtomPropConds)+
"""
function atomprop!(state::SMARTSParser, qtree::QueryTree)
    sym1 = read(state)
    vs = Int[]
    hasiso = isdigit(sym1)
    if hasiso  # Isotope
        num = atompropnumber!(state)
        push!(vs, add_qnode!(qtree, qeq(:isotope, string(num))))
    end
    v = atompropsymcond!(state, qtree)
    if v != 0
        push!(vs, v)
    end
    hasiso && v == 0 && error("(atomprop!) atom symbol not specified")
    v = atompropcond!(state, qtree)
    while v != 0
        push!(vs, v)
        v = atompropcond!(state, qtree)
    end
    isempty(vs) && error("(atomprop!) empty atomprop")
    length(vs) == 1 && return vs[1]
    node = add_qnode!(qtree, qand())
    for i in vs
        add_qedge!(qtree, node, i)
    end
    return node
end


"""
    atom!(state::SMILESParser) -> Vector{SMILESAtom}

Atom <- '[' AtomProp+ ']' / AtomSymbol
"""
function atom!(state::SMILESParser{T,V,E}) where {T,V,E}
    sym1 = read(state)
    if sym1 == '['
        forward!(state)
        props = atomprop!(state)
        read(state) == ']' || error("(atom)! unexpected token: $(read(state))")
        forward!(state)
    else
        props = atomsymbol!(state)
    end
    return props
end


"""
    atom!(state::SMARTSParser) -> Vector{QueryTree}

Atom <- '[' (AtomProp / LogicalOperator)+ ']' / AtomSymbol
"""
function atom!(state::SMARTSParser{T,V,E}) where {T,V,E}
    qtree = V()
    sym1 = read(state)
    if sym1 == '['
        forward!(state)
        v = lglowand!(state, qtree)
        v == 0 && error("(atom!) empty atomprop")
        read(state) == ']' || error("(atom!) unexpected token: $(read(state))")
        forward!(state)
    else
        v = atomsymbol!(state, qtree)
        v == 0 && return V[]  # termination token ().\0
    end
    return [qtree]
end
