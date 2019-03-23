#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    atom!,
    atomsymbol!,
    atomprop!


# TODO: radical


function atom!(state::SmartsParser{SMILES})
    """ Atom <- '[' AtomProp+ ']' / AtomSymbol
    """
    c = read(state)
    if c == '['
        forward!(state)
        a = lghighand!(state, atomprop!)
        ca = read(state)
        @assert ca == ']' "(atom!) unexpected token: $(ca)"
        forward!(state)
        prop = Dict()
        # TODO: merge logical operator
        if a.first == :and
            for p in a.second
                if p.first == :and
                    for p2 in p.second
                        prop[p2.first] = p2.second
                    end
                else
                    prop[p.first] = p.second
                end
            end
        else
            prop[a.first] = a.second
        end
        return SmilesAtom(
            prop[:Symbol],
            get!(prop, :Charge, 0),
            1,
            get!(prop, :Mass, nothing),
            get!(prop, :Aromatic, false),
            get!(prop, :stereo, nothing)
        )
    else
        a = atomsymbol!(state)
        if a === nothing
            return
        else
            sym = a.second[1].second
            arom = a.second[2].second
            return SmilesAtom(sym, 0, 1, nothing, arom, nothing)
        end
    end
end


function atom!(state::SmartsParser{SMARTS})
    """ Atom <- '[' (AtomProp / LogicalOperator)+ ']' / AtomSymbol
    """
    c = read(state)
    if c == '['
        forward!(state)
        q = lglowand!(state, atomprop!)
        cq = read(state)
        @assert cq == ']' "(atomquery!) unexpected token: $(cq)"
        forward!(state)
        return SmartsAtom(q)
    else
        a = atomsymbol!(state)
        if a === nothing
            return
        else
            return SmartsAtom(a)
        end
    end
end


function atomsymbol!(state::SmartsParser)
    """ Atomsymbol <- Br / Cl / [AaBCcNnOoPpSsFI*]
    """
    if read(state) == 'C' && lookahead(state, 1) == 'l'
        forward!(state, 2)
        return :and => (:Symbol => :Cl, :Aromatic => false)
    elseif read(state) == 'B' && lookahead(state, 1) == 'r'
        forward!(state, 2)
        return :and => (:Symbol => :Br, :Aromatic => false)
    elseif read(state) in "BCNOPSFI"
        sym = Symbol(read(state))
        forward!(state)
        return :and => (:Symbol => sym, :Aromatic => false)
    elseif read(state) in "cnops"
        sym = Symbol(uppercase(read(state)))
        forward!(state)
        return :and => (:Symbol => sym, :Aromatic => true)
    elseif read(state) == 'A'
        forward!(state)
        return :Aromatic => false
    elseif read(state) == 'a'
        forward!(state)
        return :Aromatic => true
    elseif read(state) == '*'
        forward!(state)
        return :any => true
    end
end


const SMARTS_ATOM_COND_SYMBOL = Dict(
    'X' => :Connectivity,
    'D' => :Degree,
    'v' => :Valence,
    'H' => :H_Count,
    'r' => :RingSize,
    'R' => :RingMemCount
)

const SMARTS_CHARGE_SIGN = Dict(
    '+' => 1,
    '-' => -1
)


function atomprop!(state::SmartsParser)
    """ AtomProp <- '\$(' RecursiveQuery ')' / Mass / Symbol / AtomNum /
        Stereo / CHG / [DHRrvX]
    """
    c = read(state)
    c2 = lookahead(state, 1)
    atomsyms = keys(PERIODIC_TABLE)
    if isuppercase(c)
        # Non-organic atoms
        c3 = lookahead(state, 2)
        if string(c, c2, c3) in atomsyms
            # Note: three-letter atoms (U-series) are not supported yet
            forward!(state, 3)
            return :Symbol => Symbol(c, c2, c3)
        elseif string(c, c2) in atomsyms
            forward!(state, 2)
            return :Symbol => Symbol(c, c2)
        end
    end
    # Organic atoms
    a = atomsymbol!(state)
    if a !== nothing
        return a
    end
    # Hydrogen special cases
    if c == 'H'
        cb = lookahead(state, -1)
        if c2 == ']' && (isdigit(cb) || cb == '[')
            # Hydrogen atom
            forward!(state)
            return :Symbol => :H
        elseif c2 == '+'
            # Proton
            forward!(state)
            return :and => (:Symbol => :H, :Charge => 1)
        end
    end
    # Atom properties
    if c in keys(SMARTS_ATOM_COND_SYMBOL)
        forward!(state)
        if isdigit(c2)
            num = parse(Int, c2)
            forward!(state)
        else
            num = 1
        end
        return SMARTS_ATOM_COND_SYMBOL[c] => num
    elseif isuppercase(c) && string(c) in atomsyms
        # Single letter non-organic atoms
        forward!(state, 1)
        return :Symbol => Symbol(c)
    elseif c == '#'
        # Atomic number
        forward!(state)
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = parse(Int, SubString(state.input, start, state.pos))
        forward!(state)
        return :Symbol => atomsymbol(num)
    elseif c in keys(SMARTS_CHARGE_SIGN)
        # Charge
        forward!(state)
        c2 = read(state)
        if isdigit(c2)
            chg = parse(Int, c2)
            forward!(state)
        else
            chg = 1
            while read(state) == c
                forward!(state)
                chg += 1
            end
        end
        return :Charge => chg * SMARTS_CHARGE_SIGN[c]
    elseif c == '@'
        # Stereo
        # @ => 1, @@ => 2, @? => 3, @@? => 4
        s = 1
        if lookahead(state, 1) == '@'
            forward!(state)
            s = 2
        end
        if lookahead(state, 1) == '?'
            forward!(state)
            s += 2
        end
        forward!(state)
        return :stereo => s
    elseif isdigit(c)
        # Isotope
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = SubString(state.input, start, state.pos)
        forward!(state)
        return :Mass => parse(Int, num)
    elseif c == '$' && lookahead(state, 1) == '('
        # Recursive
        forward!(state, 2)
        start = state.pos
        while read(state) != ')'
            forward!(state)
        end
        q = SubString(state.input, start, state.pos - 1)
        forward!(state)
        return :recursive => q
    end
end
