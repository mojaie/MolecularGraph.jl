#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    sssr, atom_ringmem, bond_ringmem,
    atom_ringsizes, bond_ringsizes,
    atom_isringmem, bond_isringmem,
    atom_ringcount, bond_ringcount,
    scaffoldmem, componentmem,
    atomsymbol, charge, multiplicity, bondorder, nodedegree,
    apparentvalence, valence, lonepair, heavyatomcount,
    explicithcount, implicithcount, hcount,
    connectivity, pielectron, ishdonor, ishacceptor, stdweight,
    isrotatable, isaromatic, isaromaticbond,
    molweight,
    hydrogen_acceptor_count,
    hydrogen_donor_count,
    wildman_crippen_logp,
    rotatable_count



sssr(mol) = mol.graph[:mincycles]
atom_ringmem(mol) = mol.graph[:node_cyclemem]
bond_ringmem(mol) = mol.graph[:edge_cyclemem]

function atom_ringsizes(mol)
    return [Set(length.(sssr(mol)[cs])) for cs in atom_ringmem(mol)]
end

function bond_ringsizes(mol)
    return [Set(length.(sssr(mol)[cs])) for cs in bond_ringmem(mol)]
end

# TODO: waiting for fix #28992
# @cache atom_isringmem(mol) = .!isempty.(node_cyclemem(mol))
atom_isringmem(mol) = [!i for i in isempty.(atom_ringmem(mol))]

# @cache bond_isringmem(mol) = .!isempty.(node_cyclemem(mol))
bond_isringmem(mol) = [!i for i in isempty.(bond_ringmem(mol))]

atom_ringcount(mol) = length.(atom_ringmem(mol))
bond_ringcount(mol) = length.(bond_ringmem(mol))

scaffoldmem(mol) = twoedge_membership(mol.graph)
componentmem(mol) = connected_membership(mol.graph)


"""
    atomsymbol(mol::VectorMol) -> Vector{Int}

Return the size `n` vector of atom symbols within the molecule that have `n`
atoms.
"""
@cache atomsymbol(mol::VectorMol) = getproperty.(nodevalues(mol), :symbol)


"""
    multiplicity(mol::VectorMol) -> Vector{Int}

Return the size `n` vector of atom charges within the molecule that have `n`
atoms.
"""
@cache charge(mol::VectorMol) = getproperty.(nodevalues(mol), :charge)


"""
    multiplicity(mol::VectorMol) -> Vector{Int}

Return the size `n` vector of the atom multiplicity within the molecule that
have `n` atoms (1: non-radical, 2: radical, 3: biradical).
"""
@cache multiplicity(
    mol::VectorMol) = getproperty.(nodevalues(mol), :multiplicity)


"""
    bondorder(mol::VectorMol) -> Vector{Int}

Return the size `n` vector of the bond order within the molecule that have `n`
bonds.
"""
@cache bondorder(mol::VectorMol) = getproperty.(edgevalues(mol), :order)


"""
    nodedegree(mol::VectorMol) -> Vector{Int}

Return the size `n` vector of the node degree within the molecular graph that
have `n` atom nodes.
"""
@cache nodedegree(mol::VectorMol) = [degree(mol, n) for n in nodekeys(mol)]


"""
    apparentvalence(mol::VectorMol) -> Vector{Int}

Return the size `n` vector of the apparent valence (the sum of incident bond
order) within the molecule that have `n` atoms.
"""
@cache function apparentvalence(mol::VectorMol)
    vec = zeros(Int, nodecount(mol))
    bondorder = mol[:bondorder]
    for n in nodekeys(mol)
        for b in incidences(mol, n)
            vec[n] += bondorder[b]
        end
    end
    return vec
end


"""
    valence(mol::VectorMol) -> Vector{Union{Int,Nothing}}

Return the size `n` vector of the intrinsic valence (with considering implicit
hydrogens) within the molecule that have `n` atoms.

Note that implicit hydrogens are available for only organic atoms. The intrinsic
valence value of inorganic atoms would be `nothing`.
"""
@cache function valence(mol::VectorMol)
    defs = Dict(
        :H => 1, :B => 3, :C => 4, :N => 3, :O => 2, :F => 1,
        :Si => 4, :P => 3, :S => 2, :Cl => 1,
        :As => 3, :Se => 2, :Br => 1, :I => 1
    )
    vec = Union{Int,Nothing}[]
    for (sym, chg) in zip(mol[:atomsymbol], mol[:charge])
        num = get(defs, sym, nothing)
        v = num === nothing ? nothing : num + chg
        push!(vec, v)
    end
    return vec
end


"""
    lonepair(mol::VectorMol) -> Vector{Union{Int,Nothing}}

Return the size `n` vector of the number of lone pairs within the molecule that
have `n` atoms.

Note that implicit hydrogens are available for only organic atoms. The lonepair
value of inorganic atoms would be `nothing`.
"""
@cache function lonepair(mol::VectorMol)
    defs = Dict(
        :H => 0, :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
        :Si => 0, :P => 1, :S => 2, :Cl => 3,
        :As => 1, :Se => 2, :Br => 3, :I => 3
    )
    vec = Union{Int,Nothing}[]
    for (sym, chg) in zip(mol[:atomsymbol], mol[:charge])
        num = get(defs, sym, nothing)
        v = num === nothing ? nothing : num - chg
        push!(vec, v)
    end
    return vec
end


"""
    heavyatomcount(mol::VectorMol) -> Vector{Int}

Return the size `n` vector of the number of adjacent non-hydrogen atoms within
the molecule that have `n` atoms.
"""
@cache function heavyatomcount(mol::VectorMol)
    vec = zeros(Int, nodecount(mol))
    symbols = mol[:atomsymbol]
    for n in nodekeys(mol)
        for nbr in adjacencies(mol, n)
            if symbols[nbr] != :H
                vec[n] += 1
            end
        end
    end
    return vec
end


@cache function implicithcount(mol::VectorMol)
    hcnt = (v, av) -> v === nothing ? 0 : max(0, v - av)
    return hcnt.(mol[:valence], mol[:apparentvalence])
end

@cache explicithcount(mol::VectorMol) = mol[:nodedegree] - mol[:heavyatomcount]
@cache hcount(mol::VectorMol) = mol[:explicithcount] + mol[:implicithcount]
@cache connectivity(mol::VectorMol) = mol[:nodedegree] + mol[:implicithcount]
@cache pielectron(mol::VectorMol) = mol[:apparentvalence] - mol[:nodedegree]


@cache function ishdonor(mol::VectorMol)
    dc = (sym, h) -> sym in (:N, :O) && h > 0
    return dc.(mol[:atomsymbol], mol[:hcount])
end


@cache function ishacceptor(mol::VectorMol)
    ac = (sym, lp) -> lp === nothing ? false : sym in (:N, :O, :F) && lp > 0
    return ac.(mol[:atomsymbol], mol[:lonepair])
end


function stdweight(mol::VectorMol)
    weight = (atom, imh) -> atomweight(atom) + H_WEIGHT * h
    return weight.(nodevalues(mol), mol[:implicithcount])
end


"""
    isrotatable(mol::VectorMol)

Return whether the bonds are rotatable or not.
"""
@cache function isrotatable(mol::VectorMol)
    ringmem = mol[:atom_ringmem]
    order = mol[:bondorder]
    deg = mol[:nodedegree]
    vec = Bool[]
    for (i, bond) in edgesiter(mol)
        (u, v) = (bond.u, bond.v)
        rot = (order[i] == 1 && deg[u] != 1 && deg[v] != 1
            && isempty(intersect(ringmem[u], ringmem[v])))
        push!(vec, rot)
    end
    return vec
end


@cache function isaromaticring(mol::VectorMol)
    vec = falses(circuitrank(mol.graph))
    for (i, ring) in enumerate(sssr(mol))
        if satisfyhuckel(mol, ring)
            vec[i] = true
        elseif nodetype(mol) === SmilesAtom # SMILES aromatic atom
            sub = nodesubgraph(mol.graph, Set(ring))
            if all(n.isaromatic for n in nodevalues(sub))
                vec[i] = true
            end
        end
    end
    return vec
end


"""
    isaromatic(mol::VectorMol)

Return the vector whether the atom belongs to an aromatic ring.

Note that aromaticity described here means simplified binary descriptor
(aromatic or not) based on classical Huckel's rule. This is intended for use in
some kind of pharmaceutical researches. Non-classical aromaticity such as
Moebius aromaticity is not considered.
"""
@cache function isaromatic(mol::VectorMol)
    aromatic = falses(nodecount(mol))
    for ring in sssr(mol)[isaromaticring(mol)]
        sub = nodesubgraph(mol.graph, Set(ring))
        for i in nodekeys(sub)
            aromatic[i] = true
        end
    end
    return aromatic
end


@cache function isaromaticbond(mol::VectorMol)
    aromaticbond = falses(edgecount(mol))
    for ring in sssr(mol)[isaromaticring(mol)]
        sub = nodesubgraph(mol.graph, Set(ring))
        for i in edgekeys(sub)
            aromaticbond[i] = true
        end
    end
    return aromaticbond
end


function satisfyhuckel(mol::VectorMol, ring)
    cnt = 0
    symbols = atomsymbol(mol)
    degrees = nodedegree(mol)
    pies = pielectron(mol)
    lps = lonepair(mol)
    carbonylO = findall((symbols .== :O) .* (degrees .== 1) .* (pies .== 1))
    carbonylC = Int[]
    for o in carbonylO
        c = collect(adjacencies(mol.graph, o))[1]
        if symbols[c] == :C
            push!(carbonylC, c)
        end
    end
    for r in ring
        if r in carbonylC
            continue
        elseif pies[r] == 1
            cnt += 1
        elseif lps[r] === nothing
            return false
        elseif lps[r] > 0
            cnt += 2
        elseif lps[r] < 0
            continue
        else
            return false
        end
    end
    return cnt % 4 == 2
end


"""
    molweight(mol::VectorMol; digits=2) -> Float64

Return standard molecular weight.
"""
function molweight(mol::VectorMol; digits=2)
    return round(reduce(+, stdweight(mol); init=0), digits=digits)
end


"""
    hydrogen_acceptor_count(mol::VectorMol) -> Int

Return the number of hydrogen bond acceptors (N, O and F).
"""
hydrogen_acceptor_count(mol::VectorMol) = reduce(+, ishacceptor(mol); init=0)


"""
    hydrogen_donor_count(mol::VectorMol) -> Int

Return the number of hydrogen bond donors (O and N attached to hydrogens).
"""
hydrogen_donor_count(mol::VectorMol) = reduce(+, ishdonor(mol); init=0)


"""
    wildman_crippen_logp(mol::VectorMol) -> Float64

Return predicted logP value calculated by using Wildman and Crippen method.

# Reference

1. Wildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical
   Parameters by Atomic Contributions. Journal of Chemical Information and
   Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l
"""
function wildman_crippen_logp(mol::VectorMol; digits=2)
    return round(reduce(+, wclogpcontrib(mol); init=0), digits=digits)
end


"""
    rotatable_count(mol::VectorMol) -> Int

Return the number of rotatable bonds.
"""
rotatable_count(mol::VectorMol) = reduce(+, isrotatable(mol); init=0)
