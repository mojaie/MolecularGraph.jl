#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# TODO: re-design annotation

export
    default_annotation!,
    intrinsic_annot!,
    valence_annot!,
    atom_annot!,
    group_annot!,
    rotatable!,
    aromatic!


struct IntrinsicAnnot <: Annotation end
struct ValenceAnnot <: Annotation end
struct AtomAnnot <: Annotation end
struct GroupAnnot <: Annotation end


function default_annotation!(mol::VectorMol)
    molgraph_topology!(mol)
    intrinsic_annot!(mol)
    valence_annot!(mol)
    return
end


# TODO: :RingBond
# TODO: :Degree (without hydrogen)
# TODO: :Connectivity (with hydrogen)
# TODO: :RingSize
# TODO: :RingCount (SSSR?)


function intrinsic_annot!(mol::VectorMol)
    # Symbol
    mol.v[:Symbol] = [atom.symbol for atom in mol.graph.nodes]
    # Charge
    mol.v[:Charge] = [atom.charge for atom in mol.graph.nodes]
    # Radical
    mol.v[:Multiplicity] = [atom.multiplicity for atom in mol.graph.nodes]
    # Bond order
    mol.v[:BondOrder] = [bond.order for bond in mol.graph.edges]
    mol.annotation[:Intrinsic] = IntrinsicAnnot()
    return
end


function valence_annot!(mol::VectorMol)
    required_annotation(mol, :Intrinsic)
    nbrords = Vector{Int}[]
    for (i, nbr) in enumerate(mol.graph.adjacency)
        push!(nbrords, [mol.v[:BondOrder][b] for b in values(nbr)])
    end
    # Number of bonds TODO: rename to "Degree" (exclude explicit H)
    mol.v[:NumBonds] = length.(nbrords)
    # Valence
    sm = arr -> reduce(+, arr; init=0)  # TODO: sum of empty arr is invalid
    mol.v[:Valence] = sm.(nbrords)
    # Number of pi electrons
    mol.v[:Pi] = mol.v[:Valence] - mol.v[:NumBonds]
    # Number of lone pairs
    mol.v[:LonePair] = lonepair.(mol.v[:Symbol], mol.v[:Charge])
    # Hydrogen count
    mol.v[:H_Count] = h_count.(mol.v[:Valence], mol.v[:LonePair])
    # Hydrogen bond donor count
    mol.v[:H_Donor] = h_donor.(mol.v[:Symbol], mol.v[:H_Count])
    # Hydrogen bond acceptor count
    mol.v[:H_Acceptor] = h_acceptor.(mol.v[:Symbol], mol.v[:LonePair])
    # Molecular weight including neighbor hydrogens
    mol.v[:MolWeight] = molweight.(mol.graph.nodes, mol.v[:H_Count])
    mol.annotation[:Valence] = ValenceAnnot()
    return
end


function lonepair(symbol, charge)
    defs = Dict(
        :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
        :Si => 0, :P => 1, :S => 2, :Cl => 3,
        :As => 1, :Se => 2, :Br => 3, :I => 3
    )
    num = get(defs, symbol, nothing)
    num === nothing ? nothing : num - charge
end


function h_count(valence, lonepair)
    lonepair === nothing ? 0 : max(0, 4 - abs(lonepair) - valence)
end


function h_donor(symbol, h_count)
    symbol in (:N, :O) && h_count > 0
end


function h_acceptor(symbol, lonepair)
    symbol in (:N, :O, :F) && lonepair > 0
end


function molweight(atom, h_count)
    atomweight(atom) + H_WEIGHT * h_count
end


function atom_annot!(mol)
    required_annotation(mol, :Valence)

    # Oxygen type
    # 0:atom, 1:hydroxyl, 2:ether, 3:oxonium,
    # 4:oxo, 5:oxocarbenium
    mol.v[:Oxygen] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for i in 1:atomcount(mol)
        if mol.v[:Symbol][i] != :O
            continue
        end
        if mol.v[:Pi][i] == 1
            mol.v[:Oxygen][i] = mol.v[:NumBonds][i] + 3
        else
            mol.v[:Oxygen][i] = mol.v[:NumBonds][i]
        end
    end

    # Nitrogen type
    # 0:atom, 1:primary amine, 2:sec, 3:tert, 4:quart(ammonium),
    # 5:primary imine, 6:sec, 7:iminium, 8:nitrile
    mol.v[:Nitrogen] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for i in 1:atomcount(mol)
        if mol.v[:Symbol][i] != :N
            continue
        end
        if mol.v[:Pi][i] == 2
            mol.v[:Nitrogen][i] = 8
        elseif mol.v[:Pi][i] == 1
            mol.v[:Nitrogen][i] = mol.v[:NumBonds][i] + 4
        else
            mol.v[:Nitrogen][i] = mol.v[:NumBonds][i]
        end
    end

    # Sulfer type
    # 0:atom, 1:thiol, 2:sulfide, 3:sulfonium,
    # 4:thio, 5:thiocarbenium, 6:higher valent
    mol.v[:Sulfur] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for i in 1:atomcount(mol)
        if mol.v[:Symbol][i] != :N
            continue
        end
        if mol.v[:Valence][i] > 3
            mol.v[:Sulfur][i] = 6
        elseif mol.v[:Pi][i] == 1
            mol.v[:Sulfur][i] = mol.v[:NumBonds][i] + 4
        else
            mol.v[:Sulfur][i] = mol.v[:NumBonds][i]
        end
    end

    mol.annotation[:AtomType] = AtomAnnot()
end


function group_annot!(mol)
    required_annotation(mol, :AtomType)

    # Alchohol
    # 0: methanol, 1: primary, 2: sec, 3:tert, 4: vinyl
    mol.v[:Alcohol] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Oxygen] .== 1]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] != :C
            continue
        elseif mol.v[:Pi][n] == 1
            mol.v[:Alcohol][i] = 4
        else
            mol.v[:Alcohol][n] = mol.v[:NumBonds][n] - 1
        end
    end

    # Thiol
    # 0: methanethiol, 1: primary, 2: sec, 3:tert, 4: vinyl
    mol.v[:Thiol] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Sulfur] .== 1]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] != :C
            continue
        elseif mol.v[:Pi][n] == 1
            mol.v[:Thiol][i] = 4
        else
            mol.v[:Thiol][n] = mol.v[:NumBonds][n] - 1
        end
    end

    # Carbonyl
    # 0:formaldehyde, 1: aldehyde, 2: ketone
    mol.v[:Carbonyl] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Oxygen] .== 4]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] != :C
            continue
        else
            mol.v[:Carbonyl][n] = mol.v[:NumBonds][n] - 1
        end
    end

    # Imine
    # 0:methaneimine, 1:aldimine, 2:ketimine
    mol.v[:Imine] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Nitrogen] .== 5]
        n = [i for (i, b) in nbrs if mol.v[:BondOrder][b] == 2][1]
        if mol.v[:Symbol][n] != :C
            continue
        else
            mol.v[:Imine][n] = mol.v[:NumBonds][n] - 1
        end
    end

    # Sulfoxide
    # 1:Sulfoxide 2:Sulfone
    mol.v[:Sulfoxide] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for nbrs in mol.graph.adjacency[mol.v[:Oxygen] .== 4]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] != :S
            continue
        else
            mol.v[:Sulfoxide][n] = mol.v[:Sulfoxide][n] === nothing ? 1 : 2
        end
    end

    mol.annotation[:AtomGroup] = GroupAnnot()
end


function rotatable!(mol)
    required_annotation(mol, :Topology)
    required_annotation(mol, :Valence)
    pred = (i, b) -> (
        mol.v[:BondOrder][i] == 1
        && mol.v[:NumBonds][b.u] != 1
        && mol.v[:NumBonds][b.v] != 1
        && isempty(intersect(mol.v[:Cycle][b.u], mol.v[:Cycle][b.v]))
    )
    mol.v[:Rotatable] = pred.(collect(1:bondcount(mol)), mol.graph.edges)
    return
end


function aromatic!(mol)
    required_annotation(mol, :Topology)
    required_annotation(mol, :AtomGroup)
    mol.v[:Aromatic] = falses(atomcount(mol))
    for ring in mol.annotation[:Topology].cycles
        if satisfyHuckel(mol, ring)
            mol.v[:Aromatic][ring] .= true
        end
    end
    return
end


function satisfyHuckel(mol::VectorMol, ring)
    cnt = 0
    for r in ring
        if mol.v[:Carbonyl][r] == 2
            continue
        elseif mol.v[:Pi][r] == 1
            cnt += 1
        elseif mol.v[:LonePair][r] > 0
            cnt += 2
        elseif mol.v[:LonePair][r] < 0
            continue
        else
            return false
        end
    end
    cnt % 4 == 2
end


# Graph based stereoisomer determination
#:BondEorZ
#    -> SMILES_CisTrans
#:CIP_Rule  # TODO Dont use CIP for chirality info. use Index like SMILES
#    -> SDFile_BondType
#    -> SMILES_Stereo
