#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    sssr, atom_sssrmem, bond_sssrmem,
    atom_sssrsizes, bond_sssrsizes,
    isringatom, isringbond,
    atom_sssrcount, bond_sssrcount,
    scaffoldmem, componentmem,
    atomsymbol, charge, multiplicity, bondorder, nodedegree,
    apparentvalence, valence, lonepair, heavyatomcount,
    explicithcount, implicithcount, hcount,
    connectivity, pielectron, ishdonor, ishacceptor, stdweight,
    isrotatable, isaromatic, isaromaticbond,
    molweight,
    hacceptorcount,
    hdonorcount,
    wclogp,
    rotatablecount



# Molecular graph topology

sssr = mincyclenodes
atom_sssrmem = node_cyclemem
bond_sssrmem = edge_cyclemem
scaffoldmem = two_edge_membership
componentmem = connected_membership


function atom_sssrsizes(mol::UndirectedGraph)
    vec = [Set{Int}() for i in 1:nodecount(mol)]
    rings = sssr(mol)
    for (i, cs) in enumerate(atom_sssrmem(mol))
        for c in cs
            push!(vec[i], length(rings[c]))
        end
    end
    return vec
end

atom_sssrsizes(view::SubgraphView) = atom_sssrsizes(view.graph)


function bond_sssrsizes(mol::UndirectedGraph)
    vec = [Set{Int}() for i in 1:edgecount(mol)]
    rings = sssr(mol)
    for (i, cs) in enumerate(bond_sssrmem(mol))
        for c in cs
            push!(vec[i], length(rings[c]))
        end
    end
    return vec
end

bond_sssrsizes(view::SubgraphView) = bond_sssrsizes(view.graph)


# TODO: waiting for fix #28992
# @cache atom_isringmem(mol) = .!isempty.(node_cyclemem(mol))
isringatom(mol::UndirectedGraph) = [!i for i in isempty.(atom_sssrmem(mol))]
isringatom(view::SubgraphView) = isringatom(view.graph)

# @cache bond_isringmem(mol) = .!isempty.(node_cyclemem(mol))
isringbond(mol::UndirectedGraph) = [!i for i in isempty.(bond_sssrmem(mol))]
isringbond(view::SubgraphView) = isringbond(view.graph)

atom_sssrcount(mol::UndirectedGraph) = length.(atom_sssrmem(mol))
atom_sssrcount(view::SubgraphView) = atom_sssrcount(view.graph)

bond_sssrcount(mol::UndirectedGraph) = length.(bond_sssrmem(mol))
bond_sssrcount(view::SubgraphView) = bond_sssrcount(view.graph)



# Elemental properties

"""
    atomsymbol(mol::GraphMol) -> Vector{Int}

Return the size `n` vector of atom symbols within the molecule that have `n`
atoms.
"""
@cache atomsymbol(mol::GraphMol) = getproperty.(nodeattrs(mol), :symbol)
atomsymbol(view::SubgraphView) = atomsymbol(view.graph)


"""
    charge(mol::GraphMol) -> Vector{Int}

Return the size `n` vector of atom charges within the molecule that have `n`
atoms.
"""
@cache charge(mol::GraphMol) = getproperty.(nodeattrs(mol), :charge)
charge(view::SubgraphView) = charge(view.graph)


"""
    multiplicity(mol::GraphMol) -> Vector{Int}

Return the size `n` vector of the atom multiplicity within the molecule that
have `n` atoms (1: non-radical, 2: radical, 3: biradical).
"""
@cache multiplicity(mol::GraphMol) = getproperty.(nodeattrs(mol), :multiplicity)
multiplicity(view::SubgraphView) = multiplicity(view.graph)


"""
    bondorder(mol::GraphMol) -> Vector{Int}

Return the size `n` vector of the bond order within the molecule that have `n`
bonds.
"""
@cache bondorder(mol::GraphMol) = getproperty.(edgeattrs(mol), :order)
bondorder(view::SubgraphView) = bondorder(view.graph)


"""
    nodedegree(mol::GraphMol) -> Vector{Int}

Return the size `n` vector of the node degree within the molecular graph that
have `n` atom nodes.
"""
@cache nodedegree(mol::GraphMol) = [degree(mol, n) for n in 1:nodecount(mol)]
nodedegree(view::SubgraphView) = nodedegree(view.graph)


"""
    apparentvalence(mol::GraphMol) -> Vector{Int}

Return the size `n` vector of the apparent valence (the sum of incident bond
order) within the molecule that have `n` atoms.
"""
@cache function apparentvalence(mol::GraphMol)
    vec = zeros(Int, nodecount(mol))
    bondorder_ = bondorder(mol)
    for i in 1:nodecount(mol)
        for inc in incidences(mol, i)
            vec[i] += bondorder_[inc]
        end
    end
    return vec
end

apparentvalence(view::SubgraphView) = apparentvalence(view.graph)


"""
    valence(mol::GraphMol) -> Vector{Union{Int,Nothing}}

Return the size `n` vector of the intrinsic valence (with considering implicit
hydrogens) within the molecule that have `n` atoms.

Note that implicit hydrogens are available for only organic atoms. The intrinsic
valence value of inorganic atoms would be `nothing`.
"""
@cache function valence(mol::GraphMol)
    defs = Dict(
        :H => 1, :B => 3, :C => 4, :N => 3, :O => 2, :F => 1,
        :Si => 4, :P => 3, :S => 2, :Cl => 1,
        :As => 3, :Se => 2, :Br => 1, :I => 1
    )
    vec = Union{Int,Nothing}[]
    for (sym, chg) in zip(atomsymbol(mol), charge(mol))
        num = get(defs, sym, nothing)
        v = num === nothing ? nothing : num + chg
        push!(vec, v)
    end
    return vec
end

valence(view::SubgraphView) = valence(view.graph)


"""
    lonepair(mol::GraphMol) -> Vector{Union{Int,Nothing}}

Return the size `n` vector of the number of lone pairs within the molecule that
have `n` atoms.

Note that implicit hydrogens are available for only organic atoms. The lonepair
value of inorganic atoms would be `nothing`.
"""
@cache function lonepair(mol::GraphMol)
    defs = Dict(
        :H => 0, :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
        :Si => 0, :P => 1, :S => 2, :Cl => 3,
        :As => 1, :Se => 2, :Br => 3, :I => 3
    )
    vec = Union{Int,Nothing}[]
    for (sym, chg) in zip(atomsymbol(mol), charge(mol))
        num = get(defs, sym, nothing)
        v = num === nothing ? nothing : num - chg
        push!(vec, v)
    end
    return vec
end

lonepair(view::SubgraphView) = lonepair(view.graph)


"""
    heavyatomcount(mol::GraphMol) -> Vector{Int}

Return the size `n` vector of the number of adjacent non-hydrogen atoms within
the molecule that have `n` atoms.
"""
@cache function heavyatomcount(mol::GraphMol)
    vec = zeros(Int, nodecount(mol))
    atomsymbol_ = atomsymbol(mol)
    for i in 1:nodecount(mol)
        for adj in adjacencies(mol, i)
            if atomsymbol_[adj] != :H
                vec[i] += 1
            end
        end
    end
    return vec
end

heavyatomcount(view::SubgraphView) = heavyatomcount(view.graph)


@cache function implicithcount(mol::GraphMol)
    hcnt = (v, av) -> v === nothing ? 0 : max(0, v - av)
    return hcnt.(valence(mol), apparentvalence(mol))
end

implicithcount(view::SubgraphView) = implicithcount(view.graph)


@cache explicithcount(mol::GraphMol) = nodedegree(mol) - heavyatomcount(mol)
explicithcount(view::SubgraphView) = explicithcount(view.graph)

@cache hcount(mol::GraphMol) = explicithcount(mol) + implicithcount(mol)
hcount(view::SubgraphView) = hcount(view.graph)

@cache connectivity(mol::GraphMol) = nodedegree(mol) + implicithcount(mol)
connectivity(view::SubgraphView) = connectivity(view.graph)

@cache pielectron(mol::GraphMol) = apparentvalence(mol) - nodedegree(mol)
pielectron(view::SubgraphView) = pielectron(view.graph)


@cache function ishdonor(mol::GraphMol)
    dc = (sym, h) -> sym in (:N, :O) && h > 0
    return dc.(atomsymbol(mol), hcount(mol))
end

ishdonor(view::SubgraphView) = ishdonor(view.graph)


@cache function ishacceptor(mol::GraphMol)
    ac = (sym, lp) -> lp === nothing ? false : sym in (:N, :O, :F) && lp > 0
    return ac.(atomsymbol(mol), lonepair(mol))
end

ishacceptor(view::SubgraphView) = ishacceptor(view.graph)


function stdweight(mol::GraphMol)
    weight = (atom, imh) -> atomweight(atom) + H_WEIGHT * imh
    return weight.(nodeattrs(mol), implicithcount(mol))
end

stdweight(view::SubgraphView) = stdweight(view.graph)



# Rotatable bonds

"""
    isrotatable(mol::GraphMol)

Return whether the bonds are rotatable or not.
"""
@cache function isrotatable(mol::GraphMol)
    nodedegree_ = nodedegree(mol)
    isringbond_ = isringbond(mol)
    bondorder_ = bondorder(mol)
    vec = Bool[]
    for (i, (u, v)) in enumerate(edgesiter(mol))
        rot = (!isringbond_[i] && bondorder_[i] == 1
            && nodedegree_[u] != 1 && nodedegree_[v] != 1)
        push!(vec, rot)
    end
    return vec
end

isrotatable(view::SubgraphView) = isrotatable(view.graph)



# Aromatic rings

@cache function isaromaticring(mol::GraphMol)
    vec = falses(circuitrank(mol))
    for (i, ring) in enumerate(sssr(mol))
        if satisfyhuckel(mol, ring)
            vec[i] = true
        elseif nodeattrtype(mol) === SmilesAtom # SMILES aromatic atom
            sub = nodesubgraph(mol, Set(ring))
            if all(nodeattr(mol, n).isaromatic for n in nodeset(sub))
                vec[i] = true
            end
        end
    end
    return vec
end


"""
    isaromatic(mol::GraphMol)

Return the vector whether the atom belongs to an aromatic ring.

Note that aromaticity described here means simplified binary descriptor
(aromatic or not) based on classical Huckel's rule. This is intended for use in
some kind of pharmaceutical researches. Non-classical aromaticity such as
Moebius aromaticity is not considered.
"""
@cache function isaromatic(mol::GraphMol)
    aromatic = falses(nodecount(mol))
    for ring in sssr(mol)[isaromaticring(mol)]
        sub = nodesubgraph(mol, Set(ring))
        for n in nodeset(sub)
            aromatic[n] = true
        end
    end
    return aromatic
end

isaromatic(view::SubgraphView) = isaromatic(view.graph)


@cache function isaromaticbond(mol::GraphMol)
    aromaticbond = falses(edgecount(mol))
    for ring in sssr(mol)[isaromaticring(mol)]
        sub = nodesubgraph(mol, Set(ring))
        for e in edgeset(sub)
            aromaticbond[e] = true
        end
    end
    return aromaticbond
end

isaromaticbond(view::SubgraphView) = isaromaticbond(view.graph)


function satisfyhuckel(mol::GraphMol, ring)
    cnt = 0
    atomsymbol_ = atomsymbol(mol)
    nodedegree_ = nodedegree(mol)
    pielectron_ = pielectron(mol)
    lonepair_ = lonepair(mol)
    carbonylO = findall(
        (atomsymbol_ .== :O) .* (nodedegree_ .== 1) .* (pielectron_ .== 1))
    carbonylC = Int[]
    for o in carbonylO
        c = iterate(adjacencies(mol, o))[1]
        if atomsymbol_[c] == :C
            push!(carbonylC, c)
        end
    end
    for r in ring
        if r in carbonylC
            continue
        elseif pielectron_[r] == 1
            cnt += 1
        elseif lonepair_[r] === nothing
            return false
        elseif lonepair_[r] > 0
            cnt += 2
        elseif lonepair_[r] < 0
            continue
        else
            return false
        end
    end
    return cnt % 4 == 2
end



# Molecular properties

"""
    molweight(mol::GraphMol; digits=2) -> Float64

Return standard molecular weight.
"""
molweight(mol::GraphMol; digits=2
    ) = round(reduce(+, stdweight(mol); init=0), digits=digits)


"""
    hacceptorcount(mol::GraphMol) -> Int

Return the number of hydrogen bond acceptors (N, O and F).
"""
hacceptorcount(mol::GraphMol) = reduce(+, ishacceptor(mol); init=0)


"""
    hdonorcount(mol::GraphMol) -> Int

Return the number of hydrogen bond donors (O and N attached to hydrogens).
"""
hdonorcount(mol::GraphMol) = reduce(+, ishdonor(mol); init=0)


"""
    wclogp(mol::GraphMol) -> Float64

Return predicted logP value calculated by using Wildman and Crippen method.

# Reference

1. Wildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical
   Parameters by Atomic Contributions. Journal of Chemical Information and
   Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l
"""
wclogp(mol::GraphMol; digits=2
    ) = round(reduce(+, wclogpcontrib(mol); init=0), digits=digits)


"""
    rotatablecount(mol::GraphMol) -> Int

Return the number of rotatable bonds.
"""
rotatablecount(mol::GraphMol) = reduce(+, isrotatable(mol); init=0)
