#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# TODO: radical


const SMARTS_ATOM_COND_SYMBOL = Dict(
    'X' => :connectivity,
    'D' => :nodedegree,
    'v' => :valence,
    'H' => :hydrogenconnected,
    'r' => :sssrsizes,
    'R' => :sssrcount
)

const SMARTS_CHARGE_SIGN = Dict(
    '+' => 1,
    '-' => -1
)


"""
    atomsymbol!(state::SmartsParserState) -> Union{Pair,Nothing}

Atomsymbol <- Br / Cl / [AaBCcNnOoPpSsFI*]
"""
function atomsymbol!(state::SmartsParserState)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == 'C' && sym2 == 'l'
        forward!(state, 2)
        return :and => (:atomsymbol => :Cl, :isaromatic => false)
    elseif sym1 == 'B' && sym2 == 'r'
        forward!(state, 2)
        return :and => (:atomsymbol => :Br, :isaromatic => false)
    elseif sym1 in "BCNOPSFI"
        forward!(state)
        return :and => (:atomsymbol => Symbol(sym1), :isaromatic => false)
    elseif sym1 in "bcnops"
        forward!(state)
        return :and => (:atomsymbol => Symbol(uppercase(sym1)), :isaromatic => true)
    elseif sym1 == 'A'
        forward!(state)
        return :isaromatic => false
    elseif sym1 == 'a'
        forward!(state)
        return :isaromatic => true
    elseif sym1 == '*'
        forward!(state)
        return :any => true
    end
    # return nothing
end



"""
    atomprop!(state::SmartsParserState) -> Union{Pair,Nothing}

AtomProp <- '\$(' RecursiveQuery ')' / Mass / Symbol / AtomNum / Stereo / CHG / [DHRrvX]
"""
function atomprop!(state::SmartsParserState)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if haskey(ATOMSYMBOLMAP, string(uppercase(sym1), sym2))
        # Two letter atoms
        forward!(state, 2)
        if string(uppercase(sym1), sym2) in ("As", "Se")
            return :and => (
                :atomsymbol => Symbol(uppercase(sym1), sym2),
                :isaromatic => islowercase(sym1)
            )
        end
        @assert isuppercase(sym1)
        return :atomsymbol => Symbol(sym1, sym2)
    elseif haskey(SMARTS_ATOM_COND_SYMBOL, sym1)
        # Neighbor and ring conditions
        forward!(state)
        if isdigit(sym2)
            num = parse(Int, sym2)
            forward!(state)
        else
            sym1 in ('r', 'R') && return :not => (:sssrcount => 0)
            num = 1
        end
        return SMARTS_ATOM_COND_SYMBOL[sym1] => num
    elseif sym1 in ('A', 'a', '*')
        forward!(state)
        sym1 == 'A' && return :isaromatic => false
        sym1 == 'a' && return :isaromatic => true
        sym1 == '*' && return :any => true
    elseif haskey(ATOMSYMBOLMAP, string(uppercase(sym1)))
        # Single letter atoms
        forward!(state)
        if uppercase(sym1) in "BCNOPS"
            return :and => (
                :atomsymbol => Symbol(uppercase(sym1)),
                :isaromatic => islowercase(sym1)
            )
        end
        @assert isuppercase(sym1)
        return :atomsymbol => Symbol(sym1)
    elseif sym1 == '#'
        # Atomic number
        forward!(state)
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = parse(Int, SubString(state.input, start, state.pos))
        forward!(state)
        return :atomsymbol => atomsymbol(num)
    elseif sym1 in keys(SMARTS_CHARGE_SIGN)
        # Charge
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
        return :charge => chg * SMARTS_CHARGE_SIGN[sym1]
    elseif sym1 == '@'
        # Stereo
        # @ -> anticlockwise, @@ -> clockwise, ? -> or not specified
        forward!(state)
        cw = false
        if read(state) == '@'
            forward!(state)
            cw = true
        end
        if read(state) == '?'
            forward!(state)
            return :not => (:stereo => cw ? :anticlockwise : :clockwise)
        else
            return :stereo => cw ? :clockwise : :anticlockwise
        end
    elseif isdigit(sym1)
        # Isotope
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = SubString(state.input, start, state.pos)
        forward!(state)
        return :mass => parse(Int, num)
    elseif sym1 == '$' && sym2 == '('
        # Recursive
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
        return :recursive => q
    end
    # return nothing
end



"""
    atom!(state::SmilesParser) -> Vector{SmilesAtom}

Atom <- '[' AtomProp+ ']' / AtomSymbol
"""
function atom!(state::SmilesParser)
    sym1 = read(state)
    if sym1 == '['
        forward!(state)
        fml = lghighand!(state, atomprop!)
        fml = associate_operations(fml)
        symcls = read(state)
        @assert symcls == ']' "(atom!) unexpected token: $(symcls)"
        forward!(state)
        prop = Dict(fml.first === :and ? fml.second : fml)
        @assert haskey(prop, :atomsymbol) || prop[:hydrogenconnected] >= 1
        atoms = [SmilesAtom(
            get(prop, :atomsymbol, :H),
            get(prop, :charge, 0),
            1,
            get(prop, :mass, nothing),
            get(prop, :isaromatic, false),
            get(prop, :stereo, :unspecified)
        )]
        hcnt = get(prop, :hydrogenconnected, 0)
        if !haskey(prop, :atomsymbol)
            hcnt -= 1
        end
        for h in 1:hcnt
            push!(atoms, SmilesAtom(:H))
        end
        return atoms
    else
        fml = atomsymbol!(state)
        if fml === nothing
            return SmilesAtom[]
        else
            asym = fml.second[1].second
            arom = fml.second[2].second
            return [SmilesAtom(asym, 0, 1, nothing, arom, :unspecified)]
        end
    end
end


"""
    atom!(state::SmartsParser) -> Vector{SmartsAtom}

Atom <- '[' (AtomProp / LogicalOperator)+ ']' / AtomSymbol
"""
function atom!(state::SmartsParser)
    sym1 = read(state)
    if sym1 == '['
        forward!(state)
        fml = lglowand!(state, atomprop!)
        @assert fml !== nothing "(atom!) empty atomprop"
        fml = associate_operations(fml)
        symcls = read(state)
        @assert symcls == ']' "(atom!) unexpected token: $(symcls)"
        forward!(state)
        return [SmartsAtom(fml)]
    else
        fml = atomsymbol!(state)
        if fml === nothing
            return SmartsAtom[]
        else
            return [SmartsAtom(fml)]
        end
    end
end
