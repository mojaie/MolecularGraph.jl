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
    wclogptype(mol::GraphMol)

Return Wildman-Crippen LogP atom types.
"""
@cache function wclogptype(mol::GraphMol)
    atomtypes = fill(:undef, nodecount(mol))
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    bondorder_ = bondorder(mol)
    nodedegree_ = nodedegree(mol)
    hcount_ = hcount(mol)
    pielectron_ = pielectron(mol)
    isaromatic_ = isaromatic(mol)
    isaromaticbond_ = isaromaticbond(mol)
    # Carbons
    for i in findall(atomsymbol_ .== :C)
        if pielectron_[i] == 0
            # Aliphatic (C1-4,8-12,27)
            nbrs = collect(adjacencies(mol, i))
            if !isempty(setdiff(atomsymbol_[nbrs], ALIPH_HETERO))
                atomtypes[i] = :C27 # Adjacent to inorganic atoms
            elseif all(isaromatic_[nbrs] .== false)
                if all(atomsymbol_[nbrs] .== :C)
                    # Aliphatic carbon (C1,2)
                    if nodedegree_[i] <= 2
                        atomtypes[i] = :C1
                    else
                        atomtypes[i] = :C2
                    end
                else
                    # Adjacent to heteroatoms (C3,4)
                    if nodedegree_[i] <= 2
                        atomtypes[i] = :C3
                    else
                        atomtypes[i] = :C4
                    end
                end
            else
                # Adjacent to aromatic atoms (C8-12)
                if nodedegree_[i] == 1
                    nbr = nbrs[1]
                    if atomsymbol_[nbr] == :C
                        atomtypes[i] = :C8
                    else
                        atomtypes[i] = :C9
                    end
                elseif nodedegree_[i] == 2
                    atomtypes[i] = :C10
                elseif nodedegree_[i] == 3
                    atomtypes[i] = :C11
                elseif nodedegree_[i] == 4
                    atomtypes[i] = :C12
                end
            end
        elseif isaromatic_[i]
            # Aromatic (C13-25)
            if degree(mol, i) == 2
                atomtypes[i] = :C18 # Aromatic non-substituted
                continue
            end
            subst = nothing
            substbond = nothing
            aromcnt = 0
            for (inc, adj) in neighbors(mol, i)
                if isaromaticbond_[inc]
                    aromcnt += 1
                else
                    subst = adj
                    substbond = inc
                end
            end
            if aromcnt == 3
                atomtypes[i] = :C19 # Bridgehead
            elseif !haskey(AROM_HETERO, atomsymbol_[subst])
                atomtypes[i] = :C13 # Inorganic substituent
            elseif isaromatic_[subst]
                atomtypes[i] = :C20 # Aromatic substituent
            elseif bondorder_[substbond] == 2
                atomtypes[i] = :C25 # Double bond substituent
            else
                # Typical substituent (C14-17,21-24)
                atomtypes[i] = AROM_HETERO[atomsymbol_[subst]]
            end
        else
            # Aliphatic multiple bond (C5-7,26)
            nbrs = collect(adjacencies(mol, i))
            bonds = collect(incidences(mol, i))
            if any(bondorder_[bonds] .== 3)
                atomtypes[i] = :C7 # Alkyne, Nitrile
            elseif any((pielectron_[nbrs] .> 0) .* (atomsymbol_[nbrs] .!= :C))
                atomtypes[i] = :C5 # Double bond to non-C
            elseif any(isaromatic_[nbrs] .== true)
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
    for i in findall(atomsymbol_ .== :N)
        if charge_[i] > 0
            nbrs = collect(adjacencies(mol, i))
            if isaromatic_[i]
                atomtypes[i] = :N12 # Charged aromatic nitrogen
            elseif hcount_[i] == 0
                if nodedegree_[i] == 2 && all(atomsymbol_[nbrs] .== :N)
                    atomtypes[i] = :N14 # Azide is exceptionally N14
                else
                    atomtypes[i] = :N13 # Quart-ammonium
                end
            elseif pielectron_[i] == 0
                atomtypes[i] = :N10 # Protonated amine
            else
                atomtypes[i] = :N14 # Other charged
            end
        elseif charge_[i] < 0
            atomtypes[i] = :N14 # Neg charged
        elseif isaromatic_[i]
            atomtypes[i] = :N11 # Neutral aromatic
        elseif pielectron_[i] == 2
            atomtypes[i] = :N9 # Nitrile
        elseif pielectron_[i] == 1
            # Imine (N5,6)
            if nodedegree_[i] == 1
                atomtypes[i] = :N5
            elseif nodedegree_[i] == 2
                atomtypes[i] = :N6
            end
        else
            nbrs = collect(adjacencies(mol, i))
            if all(isaromatic_[nbrs] .== false)
                # Aliphatic amine (N1,2,7)
                if nodedegree_[i] == 1
                    atomtypes[i] = :N1
                elseif nodedegree_[i] == 2
                    atomtypes[i] = :N2
                elseif nodedegree_[i] == 3
                    atomtypes[i] = :N7
                end
            else
                # Aromatic amine (N3,4,8)
                if nodedegree_[i] == 1
                    atomtypes[i] = :N3
                elseif nodedegree_[i] == 2
                    atomtypes[i] = :N4
                elseif nodedegree_[i] == 3
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
    for i in findall(atomsymbol_ .== :O)
        if isaromatic_[i]
            # Aromatic oxygen (O1)
            atomtypes[i] = :O1
            continue
        end
        nbrs = collect(adjacencies(mol, i))
        bonds = collect(incidences(mol, i))
        if nodedegree_[i] == 1
            # Hydroxyl (O2,5-12)
            nbr = nbrs[1]
            bond = bonds[1]
            if atomtypes[nbr] == :C5 && charge_[i] < 0
                # Acid (O12)
                atomtypes[i] = :O12
            elseif bondorder_[bond] == 2 || charge_[i] < 0
                # Carbonyl or oxide (O5-11)
                if atomtypes[nbr] == :C25
                    atomtypes[i] = :O8 # Aromatic carbonyl
                elseif atomtypes[nbr] == :C5
                    # Aliphatic carbonyl
                    cnbrs = adjacencies(mol, nbr)
                    pop!(cnbrs, i)
                    cnbrs = collect(cnbrs)
                    if all(atomsymbol_[cnbrs] .!= :C)
                        atomtypes[i] = :O11 # Carbonyl heteroatom
                    elseif any(isaromatic_[cnbrs] .== true)
                        atomtypes[i] = :O10 # Carbonyl aromatic
                    else
                        atomtypes[i] = :O9 # Carbonyl aliphatic
                    end
                elseif atomsymbol_[nbr] in (:O, :N)
                    atomtypes[i] = :O5 # O2? or N-oxide
                elseif atomsymbol_[nbr] == :S
                    atomtypes[i] = :O6 # S-oxide
                else
                    atomtypes[i] = :O7 # Other oxide
                end
            else
                # Alcohol (O2)
                atomtypes[i] = :O2
            end
        elseif nodedegree_[i] == 2
            # Ether (O3,4)
            if all(isaromatic_[nbrs] .== false)
                atomtypes[i] = :O3 # Aliphatic ether
            else
                atomtypes[i] = :O4 # Aromatic ether
            end
        end
    end

    # Others
    for i in findall(atomtypes .== :undef)
        if atomsymbol_[i] in (:F, :Cl, :Br, :I)
            if charge_[i] == 0
                atomtypes[i] = atomsymbol_[i]
            else
                atomtypes[i] = :Hal
            end
        elseif atomsymbol_[i] == :P
            atomtypes[i] = :P
        elseif atomsymbol_[i] == :S
            if isaromatic_[i]
                atomtypes[i] = :S3
            elseif charge_[i] == 0
                atomtypes[i] = :S1
            else
                atomtypes[i] = :S2
            end
        elseif atomsymbol_[i] == :H
            nbr = iterate(adjacencies(mol, i))[1]
            atomtypes[i] = wclogphydrogentype(mol, nbr)
        elseif atomsymbol_[i] in P_BLOCK
            atomtypes[i] = :Me1
        elseif atomsymbol_[i] in D_BLOCK
            atomtypes[i] = :Me2
        end
    end
    return atomtypes
end


function wclogphydrogentype(mol::GraphMol, i)
    atomsymbol_ = atomsymbol(mol)
    pielectron_ = pielectron(mol)
    isaromatic_ = isaromatic(mol)
    if atomsymbol_[i] == :C
        return :H1 # Hydrocarbon
    elseif atomsymbol_[i] == :N
        return :H3 # Amine
    elseif atomsymbol_[i] == :O
        nbr = iterate(adjacencies(mol, i))[1]
        if atomsymbol_[nbr] == :N
            return :H3 # Hydroxyamine
        elseif atomsymbol_[nbr] in (:O, :S)
            return :H4 # Peroxide, sulfoxide
        elseif atomsymbol_[nbr] == :C
            if pielectron_[nbr] == 1 && !isaromatic_[nbr]
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


function wclogpcontrib(mol::GraphMol)
    contrib = zeros(nodecount(mol))
    atomsymbol_ = atomsymbol(mol)
    hcount_ = hcount(mol)
    atomtypes = wclogptype(mol)
    for i in 1:nodecount(mol)
        hcnt = hcount_[i]
        cont = get(WCLOGP_TABLE, string(atomtypes[i]), 0)
        if atomsymbol_[i] == :H
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
