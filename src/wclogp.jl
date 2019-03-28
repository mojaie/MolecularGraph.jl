#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    wclogptype,
    wclogphydrogentype,
    wclogpcontrib


const WCLOGP_TABLE = YAML.load(open(
    joinpath(dirname(@__FILE__), "..", "assets", "const", "wclogp.yaml")
))["logP"]

const ALIPH_HETERO = Set([:H, :C, :N, :O, :P, :S, :F, :Cl, :Br, :I])
const AROM_HETERO = Dict(
    :H  => :C18, :C => :C21, :N => :C22, :O => :C23, :S => :C24,
    :F => :C14, :Cl => :C15, :Br => :C16, :I => :C17
)
# TODO:
const P_BLOCK = (:Al, :B, :Si, :Ga, :Ge, :As, :Se, :Sn, :Te, :Pb,
                 :Ne, :Ar, :Kr, :Xe, :Rn)
const D_BLOCK = (:Fe, :Co, :Cu, :Zn, :Tc, :Cd, :Pt, :Au, :Hg, :Gd)


"""
    wclogptype(mol::VectorMol)

Return Wildman-Crippen LogP atom types.
"""
@cache function wclogptype(mol::VectorMol)
    atomtypes = fill(:undef, atomcount(mol))
    symbols = atomsymbol(mol)
    charges = charge(mol)
    order = bondorder(mol)
    degrees = nodedegree(mol)
    hydcnt = hcount(mol)
    pies = pielectron(mol)
    isarom = isaromatic(mol)
    isarombond = isaromaticbond(mol)
    # Carbons
    for i in findall(symbols .== :C)
        if pies[i] == 0
            # Aliphatic (C1-4,8-12,27)
            nbrs = collect(adjacencies(mol, i))
            if !isempty(setdiff(symbols[nbrs], ALIPH_HETERO))
                atomtypes[i] = :C27 # Adjacent to inorganic atoms
            elseif all(isarom[nbrs] .== false)
                if all(symbols[nbrs] .== :C)
                    # Aliphatic carbon (C1,2)
                    if degrees[i] <= 2
                        atomtypes[i] = :C1
                    else
                        atomtypes[i] = :C2
                    end
                else
                    # Adjacent to heteroatoms (C3,4)
                    if degrees[i] <= 2
                        atomtypes[i] = :C3
                    else
                        atomtypes[i] = :C4
                    end
                end
            else
                # Adjacent to aromatic atoms (C8-12)
                if degrees[i] == 1
                    nbr = nbrs[1]
                    if symbols[nbr] == :C
                        atomtypes[i] = :C8
                    else
                        atomtypes[i] = :C9
                    end
                elseif degrees[i] == 2
                    atomtypes[i] = :C10
                elseif degrees[i] == 3
                    atomtypes[i] = :C11
                elseif degrees[i] == 4
                    atomtypes[i] = :C12
                end
            end
        elseif isarom[i]
            # Aromatic (C13-25)
            if degree(mol, i) == 2
                atomtypes[i] = :C18 # Aromatic non-substituted
                continue
            end
            subst = nothing
            substbond = nothing
            aromcnt = 0
            for (nbr, b) in neighbors(mol, i)
                if isarombond[b]
                    aromcnt += 1
                else
                    subst = nbr
                    substbond = b
                end
            end
            if aromcnt == 3
                atomtypes[i] = :C19 # Bridgehead
            elseif !haskey(AROM_HETERO, symbols[subst])
                atomtypes[i] = :C13 # Inorganic substituent
            elseif isarom[subst]
                atomtypes[i] = :C20 # Aromatic substituent
            elseif order[substbond] == 2
                atomtypes[i] = :C25 # Double bond substituent
            else
                # Typical substituent (C14-17,21-24)
                atomtypes[i] = AROM_HETERO[symbols[subst]]
            end
        else
            # Aliphatic multiple bond (C5-7,26)
            nbrs = collect(adjacencies(mol, i))
            bonds = collect(incidences(mol, i))
            if any(order[bonds] .== 3)
                atomtypes[i] = :C7 # Alkyne, Nitrile
            elseif any((pies[nbrs] .> 0) .* (symbols[nbrs] .!= :C))
                atomtypes[i] = :C5 # Double bond to non-C
            elseif any(isarom[nbrs] .== true)
                atomtypes[i] = :C26 # Adjacent to aromatic
            else
                atomtypes[i] = :C6 # Double bond to C (including allene)
            end
        end
        if atomtypes[i] == :undef
            atomtypes[i] = :CS # Not found
        end
    end

    # Nitrogens
    for i in findall(symbols .== :N)
        if charges[i] > 0
            nbrs = collect(adjacencies(mol, i))
            if isarom[i]
                atomtypes[i] = :N12 # Charged aromatic nitrogen
            elseif hydcnt[i] == 0
                if degrees[i] == 2 && all(symbols[nbrs] .== :N)
                    atomtypes[i] = :N14 # Azide is exceptionally N14
                else
                    atomtypes[i] = :N13 # Quart-ammonium
                end
            elseif pies[i] == 0
                atomtypes[i] = :N10 # Protonated amine
            else
                atomtypes[i] = :N14 # Other charged
            end
        elseif charges[i] < 0
            atomtypes[i] = :N14 # Neg charged
        elseif isarom[i]
            atomtypes[i] = :N11 # Neutral aromatic
        elseif pies[i] == 2
            atomtypes[i] = :N9 # Nitrile
        elseif pies[i] == 1
            # Imine (N5,6)
            if degrees[i] == 1
                atomtypes[i] = :N5
            elseif degrees[i] == 2
                atomtypes[i] = :N6
            end
        else
            nbrs = collect(adjacencies(mol, i))
            if all(isarom[nbrs] .== false)
                # Aliphatic amine (N1,2,7)
                if degrees[i] == 1
                    atomtypes[i] = :N1
                elseif degrees[i] == 2
                    atomtypes[i] = :N2
                elseif degrees[i] == 3
                    atomtypes[i] = :N7
                end
            else
                # Aromatic amine (N3,4,8)
                if degrees[i] == 1
                    atomtypes[i] = :N3
                elseif degrees[i] == 2
                    atomtypes[i] = :N4
                elseif degrees[i] == 3
                    atomtypes[i] = :N8
                end
            end
        end
        if atomtypes[i] == :undef
            # Not found
            atomtypes[i] = :NS
        end
    end

    # Oxygens
    for i in findall(symbols .== :O)
        if isarom[i]
            # Aromatic oxygen (O1)
            atomtypes[i] = :O1
            continue
        end
        nbrs = collect(adjacencies(mol, i))
        bonds = collect(incidences(mol, i))
        if degrees[i] == 1
            # Hydroxyl (O2,5-12)
            nbr = nbrs[1]
            bond = bonds[1]
            if atomtypes[nbr] == :C5 && charges[i] < 0
                # Acid (O12)
                atomtypes[i] = :O12
            elseif order[bond] == 2 || charges[i] < 0
                # Carbonyl or oxide (O5-11)
                if atomtypes[nbr] == :C25
                    atomtypes[i] = :O8 # Aromatic carbonyl
                elseif atomtypes[nbr] == :C5
                    # Aliphatic carbonyl
                    cnbrs = adjacencies(mol, nbr)
                    pop!(cnbrs, i)
                    cnbrs = collect(cnbrs)
                    if all(symbols[cnbrs] .!= :C)
                        atomtypes[i] = :O11 # Carbonyl heteroatom
                    elseif any(isarom[cnbrs] .== true)
                        atomtypes[i] = :O10 # Carbonyl aromatic
                    else
                        atomtypes[i] = :O9 # Carbonyl aliphatic
                    end
                elseif symbols[nbr] in (:O, :N)
                    atomtypes[i] = :O5 # O2? or N-oxide
                elseif symbols[nbr] == :S
                    atomtypes[i] = :O6 # S-oxide
                else
                    atomtypes[i] = :O7 # Other oxide
                end
            else
                # Alcohol (O2)
                atomtypes[i] = :O2
            end
        elseif degrees[i] == 2
            # Ether (O3,4)
            if all(isarom[nbrs] .== false)
                atomtypes[i] = :O3 # Aliphatic ether
            else
                atomtypes[i] = :O4 # Aromatic ether
            end
        end
    end

    # Others
    for i in findall(atomtypes .== :undef)
        if symbols[i] in (:F, :Cl, :Br, :I)
            if charges[i] == 0
                atomtypes[i] = symbols[i]
            else
                atomtypes[i] = :Hal
            end
        elseif symbols[i] == :P
            atomtypes[i] = :P
        elseif symbols[i] == :S
            if isarom[i]
                atomtypes[i] = :S3
            elseif charges[i] == 0
                atomtypes[i] = :S1
            else
                atomtypes[i] = :S2
            end
        elseif symbols[i] == :H
            nbr = pop!(adjacencies(mol, i))
            atomtypes[i] = wclogphydrogentype(mol, nbr)
        elseif symbols[i] in P_BLOCK
            atomtypes[i] = :Me1
        elseif symbols[i] in D_BLOCK
            atomtypes[i] = :Me2
        end
    end
    return atomtypes
end


function wclogphydrogentype(mol::VectorMol, i)
    symbols = atomsymbol(mol)
    pies = pielectron(mol)
    isarom = isaromatic(mol)
    if symbols[i] == :C
        return :H1 # Hydrocarbon
    elseif symbols[i] == :N
        return :H3 # Amine
    elseif symbols[i] == :O
        nbr = pop!(adjacencies(mol, i))
        if symbols[nbr] == :N
            return :H3 # Hydroxyamine
        elseif symbols[nbr] in (:O, :S)
            return :H4 # Peroxide, sulfoxide
        elseif symbols[nbr] == :C
            if pies[nbr] == 1 && !isarom[nbr]
                return :H4 # Acid
            else
                return :H2 # Alcohol
            end
        else
            return :H2 # Other hydroxyl
        end
    else
        return :H2 # Other hydrate
    end
    return :HS # HS may be impossible
end


function wclogpcontrib(mol::VectorMol)
    contrib = zeros(atomcount(mol))
    symbols = atomsymbol(mol)
    hydcnt = hcount(mol)
    atomtypes = wclogptype(mol)
    for i in nodekeys(mol)
        hcnt = hydcnt[i]
        cont = get(WCLOGP_TABLE, string(atomtypes[i]), 0)
        if symbols[i] == :H
            contrib[i] = 0  # avoid double count
        elseif hcnt == 0
            contrib[i] = cont
        else
            htype = wclogphydrogentype(mol, i)
            hcont = WCLOGP_TABLE[string(htype)] * hcnt
            contrib[i] = cont + hcont
        end
    end
    return contrib
end
