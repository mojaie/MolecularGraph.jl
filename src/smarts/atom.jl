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
    error("smiles parser error - #$(state.pos) unexpected token '$(sym1)'")
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
    error("smarts parser error - #$(state.pos) unexpected token '$(sym1)'")
end


function atompropsymbol!(
        state::SMILESParser, sym::Symbol, slow::Bool, cond::NTuple)
    if sym in cond
        return (symbol=sym, isaromatic=slow)
    end
    slow && error("smiles parser error - #$(state.pos) aromatic '$(sym)' not supported")
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
    slow && error("smarts parser error - #$(state.pos) aromatic '$(sym)' not supported")
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

Return `(key=k, value=v, stereo=s)` where `k` is one of `:total_hydrogens`, `:charge` and `:stereo`
(a NamedTuple of a fixed type keeps the return type inferable).
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
        return (key=:total_hydrogens, value=num, stereo=:unspecified)
    elseif sym1 in keys(SMARTS_CHARGE_SIGN)  # Charge
        chg = atompropcharge!(state, sym1, sym2)
        return (key=:charge, value=chg * SMARTS_CHARGE_SIGN[sym1], stereo=:unspecified)
    elseif sym1 == '@'
        # Stereo: @ -> anticlockwise, @@ -> clockwise
        forward!(state)
        cw = :anticlockwise
        if read(state) == '@'
            forward!(state)
            cw = :clockwise
        end
        return (key=:stereo, value=0, stereo=cw)
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
            start = state.pos
            while isdigit(lookahead(state, 1))
                forward!(state)
            end
            num = parse(Int, SubString(state.input, start, state.pos))
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
            return add_qnode!(qtree, qeq(:stereo, cw ? "clockwiseuns" : "anticlockwiseuns"))
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
            if state.done
                error("smarts parser error - #$(state.pos) unclosed recursive")
            end
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
    atomprop!(state::SMILESParser{T,V,E}) -> Vector{V}

AtomProp <- ((Mass? / AtomPropSymConds) / AtomPropConds)+
"""
function atomprop!(state::SMILESParser{T,V,E}) where {T,V,E}
    # Properties are accumulated in typed local variables (later ones override earlier ones)
    # so that the atom constructor can be statically dispatched.
    sym1 = read(state)
    symbol = :H
    hassymbol = false
    charge = 0
    isotope = 0
    isaromatic = false
    stereo = :unspecified
    hasprop = false
    hcnt = 0
    hasiso = isdigit(sym1)
    if hasiso  # Isotope
        isotope = atompropnumber!(state)
        hasprop = true
    end
    v = atompropsymcond!(state)
    if !isnothing(v)
        symbol = v.symbol
        hassymbol = true
        haskey(v, :isaromatic) && (isaromatic = v.isaromatic)
        hasprop = true
    end
    hasiso && isnothing(v) && error("smiles parser error - #$(state.pos) missing atom symbol")
    while true
        c = atompropcond!(state)
        isnothing(c) && break
        if c.key === :total_hydrogens
            hcnt = c.value
        elseif c.key === :charge
            charge = c.value
            hasprop = true
        else  # :stereo
            stereo = c.stereo
            hasprop = true
        end
    end
    !hasprop && hcnt == 0 && error("smiles parser error - #$(state.pos) missing atom prop")
    # explicit hydrogen (e.g. [CH3]) -> hydrogen nodes
    if !hassymbol  # special case: [H2]
        hcnt -= 1
    end
    atom = V(; symbol=symbol, charge=charge, isotope=isotope, isaromatic=isaromatic, stereo=stereo)
    atoms = V[atom]
    for _ in 1:hcnt
        push!(atoms, V(;symbol=:H))
    end
    return atoms
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
    hasiso && v == 0 && error("smarts parser error - #$(state.pos) missing atom symbol")
    v = atompropcond!(state, qtree)
    while v != 0
        push!(vs, v)
        v = atompropcond!(state, qtree)
    end
    isempty(vs) && error("smarts parser error - #$(state.pos) missing atom prop")
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
        push!(state.bracket, state.node + 1)  # the first atom is the bracket atom
        props = atomprop!(state)
        read(state) == ']' || error("smiles parser error - #$(state.pos) unexpected token '$(read(state))'")
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
        v == 0 && error("smarts parser error - #$(state.pos) missing atom prop")
        read(state) == ']' || error("smarts parser error - #$(state.pos) unexpected token '$(read(state))'")
        forward!(state)
    else
        v = atomsymbol!(state, qtree)
        v == 0 && return V[]  # termination token ().\0
    end
    return [qtree]
end
