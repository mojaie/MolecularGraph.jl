#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    default_annotation!


struct IntrinsicAnnot <: Annotation end
struct ValenceAnnot <: Annotation end


function default_annotation!(mol::Molecule)
    molgraph_topology!(mol)
    intrinsic_annot!(mol)
    valence_annot!(mol)
    group_annot!(mol)
    rotatable!(mol)
    aromatic!(mol)
    # TODO: standardize chirality flag
    return
end


function remove_hydrogen(mol::MutableMolecule)
    # TODO: ignore hydrogen
    # TODO: stash stereo hydrogen (annotation :SDF_Stereo_Hydrogen)
end


function intrinsic_annot!(mol::Molecule)
    # Symbol
    mol.v[:Symbol] = [atom.symbol for atom in mol.graph.nodes]
    # Charge
    mol.v[:Charge] = [atom.charge for atom in mol.graph.nodes]
    # Radical
    mol.v[:Multiplicity] = [atom.multiplicity for atom in mol.graph.nodes]
    # Mass
    mol.v[:Mass] = [atom.mass for atom in mol.graph.nodes]
    # Bond order
    mol.v[:BondOrder] = [bond.order for bond in mol.graph.edges]
    mol.annotation[:Intrinsic] = IntrinsicAnnot()
    return
end


function valence_annot!(mol::Molecule)
    required_annotation(mol, :Intrinsic)
    nbrords = Vector{Int}[]
    for (i, nbr) in enumerate(mol.graph.adjacency)
        push!(nbrords, [mol.v[:BondOrder][b] for b in values(nbr)])
    end
    # Number of bonds
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
    mol.v[:MolWeight] = molweight.(mol.v[:Symbol], mol.v[:H_Count])
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


function molweight(symbol, h_count)
    atomweight(symbol) + H_WEIGHT * h_count
end


function group_annot!(mol)
    required_annotation(mol, :Valence)
    mol.v[:Carbonyl] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    carbonylO = Vector{Bool}(
        (mol.v[:Symbol] .== :O)
        .* (mol.v[:NumBonds] .== 1)
        .* (mol.v[:Valence] .== 2)
    )
    for nbrs in mol.graph.adjacency[carbonylO]
        n = collect(keys(nbrs))[1]
        if mol.v[:Symbol][n] == :C
            mol.v[:Carbonyl][n] = mol.v[:NumBonds][n] - 1
        end
    end
    # TODO: fragments
    # mol.annotmap[:F_Alcohol] C-O
    # mol.annotmap[:F_Thiol] C-S
    # search NumBond=1 O, pi=0
    # neighbor[1]=C C.nbr=2->primary  3->sec  4->tert

    # mol.annotmap[:F_Amine] N
    # search pi=0 N
    # NumBond=1->primary  2->sec  3->tert 4->quart
end


function rotatable!(mol)
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
    mol.v[:Aromatic] = falses(atomcount(mol))
    for ring in mol.annotation[:Topology].cycles
        if satisfyHuckel(mol, ring)
            mol.v[:Aromatic][ring] .= true
        end
    end
    return
end


function satisfyHuckel(mol::Molecule, ring)
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
