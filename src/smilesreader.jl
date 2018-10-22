#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    smilestomol


function smilestomol(smiles::AbstractString; precalc=true)
    mol = MolecularGraph()
    tokens = tokenize(smiles)
    ringclose = Dict()
    parsegroup!(mol, ringclose, tokens)
    # println(Int64[a.index for a in atomvector(mol)])
    # println(Tuple{Int64, Int64}[(b.u, b.v) for b in bondvector(mol)])
    if precalc
        assign_descriptors!(mol)
    end
    mol
end


function tokenize(smiles::AbstractString)
    tokens = []
    i = 1
    while i != lastindex(smiles) + 1
        m = match(
            r"^[\.\(\)]|^([=#:/\\]?(\[.+?\]|Br|Cl|[BCcNnOoPpSsFI])([=#:/\\]?[1-9])?)",
            smiles[i:end]
        )
        if m === nothing
            throw(OperationError("Invalid SMILES"))
        end
        i += length(m.match)
        push!(tokens, m.match)
    end
    tokens
end


function parsegroup!(mol, ringclose, tokens, base)
    pos = base
    nobond = base == 0
    while length(tokens) > 0
        token = popfirst!(tokens)
        if token == "("
            parsegroup!(mol, ringclose, tokens, pos)
            continue
        elseif token == ")"
            break
        elseif token == "."
            nobond = true
            continue
        end

        newpos = length(atomvector(mol)) + 1
        # Atom
        # println(objectid(mol), " ", pos, " ", newpos, " ", token)
        (atok, btok, ringc) = parsetoken(token)
        atom = Atom()
        atom.index = newpos
        parseatom!(atom, atok)
        newatom!(mol, atom)

        # Ring bond
        if ringc !== nothing
            num = ringc[end]
            bondtype = length(ringc) == 2 ? SubString(ringc, 1, 1) : nothing
            if num in keys(ringclose)
                (upos, ubd) = ringclose[num]
                rbond = Bond()
                rbond.u = upos
                rbond.v = newpos
                rtype = ubd === nothing ? bondtype : ubd
                parsebond!(rbond, rtype)
                updatebond!(mol, rbond)
                delete!(ringclose, num)
            else
                ringclose[num] = (newpos, bondtype)
            end
        end

        # Bond
        if nobond
            nobond = false
        else
            bond = Bond()
            bond.u = pos
            bond.v = newpos
            parsebond!(bond, btok)
            updatebond!(mol, bond)
        end
        pos = newpos
    end
    return
end

parsegroup!(mol, ringclose, tokens) = parsegroup!(mol, ringclose, tokens, 0)


function parsetoken(token::AbstractString)
    m = match(
        r"^([=#:/\\])?(\[.+?\]|Br|Cl|[BCcNnOoPpSsFI])([=#:/\\]?[1-9])?",
        token
    )
    bond = m.captures[1]
    atom = replace(m.captures[2], r"\[(.+?)\]" => s"\1")
    ringc = m.captures[3]
    (atom, bond, ringc)
end


function parseatom!(atom::Atom, token::AbstractString)
    m = match(
        r"^([0-9]*)([cnposA-Z][a-z]*)(@*)(H[2-4]?)?([\+\-]+[1-4]?)?",
        token
    )
    atom.symbol = uppercasefirst(m.captures[2])
    if islowercase(m.captures[2][1])
        atom.smiles_aromatic = true
    end
    chgconv = Dict(
        nothing => 0, "+" => 1, "++" => 2, "+++" => 3, "++++" => 4,
        "-" => -1, "--" => -2, "---" => -3, "----" => -4,
        "+" => 1, "+2" => 2, "+3" => 3, "+4" => 4,
        "-" => -1, "-2" => -2, "-3" => -3, "-4" => -4
    )
    atom.charge = chgconv[m.captures[5]]
    hcconv = Dict{Any, UInt8}(nothing => 0, "H" => 1, "H2" => 2, "H3" => 3, "H4" => 4)
    sethydrogen!(atom, hcconv[m.captures[4]])
    mass = m.captures[1]
    atom.mass = mass == "" ? nothing : parse(UInt8, mass)
    atom.smiles_stereo = m.captures[3]
    # TODO: atom initialization
    atom.visible = atom.symbol != "C" || atom.mass !== nothing
    atom.Hacceptor = atom.symbol in ("N", "O", "F")
    return
end


function parsebond!(bond::Bond, token::Union{AbstractString, Nothing})
    orderconv = Dict(
        nothing => 1, "=" => 2, "#" => 3, ":" => 1, "/" => 1, "\\" => 1
    )
    bond.order = orderconv[token]
    if token == ":"
        bond.aromatic = true
    elseif token == "/"
        bond.smiles_cis_trans = 1
    elseif token == "\\"
        bond.smiles_cis_trans = 2
    end
    return
end
