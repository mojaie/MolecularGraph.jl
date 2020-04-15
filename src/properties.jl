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
    atomsymbol, charge, multiplicity, bondorder,
    nodedegree, valence, lonepair, heavyatomcount,
    explicithcount, implicithcount, hcount, connectivity,
    countatoms, molecularformula, empiricalformula,
    ishdonor, ishacceptor, stdweight,
    isrotatable, pielectron, hybridization,
    isaromaticring, isaromatic, isaromaticbond,
    molweight, hacceptorcount, hdonorcount,
    wclogp, rotatablecount



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

Return the size ``n`` vector of atom symbols within the molecule that have ``n``
atoms.
"""
@cache atomsymbol(mol::GraphMol) = getproperty.(nodeattrs(mol), :symbol)
atomsymbol(view::SubgraphView) = atomsymbol(view.graph)


"""
    charge(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of atom charges within the molecule that have ``n``
atoms.
"""
@cache charge(mol::GraphMol) = getproperty.(nodeattrs(mol), :charge)
charge(view::SubgraphView) = charge(view.graph)


"""
    multiplicity(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of atom multiplicities within the molecule that
have ``n`` atoms (1: non-radical, 2: radical, 3: biradical).
"""
@cache multiplicity(mol::GraphMol) = getproperty.(nodeattrs(mol), :multiplicity)
multiplicity(view::SubgraphView) = multiplicity(view.graph)


"""
    bondorder(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of bond orders within the molecule that have ``n``
bonds.
"""
@cache bondorder(mol::GraphMol) = getproperty.(edgeattrs(mol), :order)
bondorder(view::SubgraphView) = bondorder(view.graph)


"""
    nodedegree(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of node degrees within the molecular graph that
have ``n`` atom nodes.

`nodedegree` values correspond to SMARTS `D` property.
"""
@cache nodedegree(mol::GraphMol) = [degree(mol, n) for n in 1:nodecount(mol)]
nodedegree(view::SubgraphView) = nodedegree(view.graph)



# Valence

"""
    lonepair(mol::GraphMol) -> Vector{Union{Int,Nothing}}

Return the size ``n`` vector of the number of lone pairs within the molecule that
have ``n`` atoms.

Note that implicit hydrogens are available for only organic atoms. The lonepair
value of inorganic atoms would be `nothing`. The result can take negative value if the atom has empty valence shells.
"""
function lonepair(mol::GraphMol)
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


function apparentvalence(mol::GraphMol)
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

Return the size ``n`` vector of intrinsic valences (with considering implicit
hydrogens) within the molecule that have ``n`` atoms.

`valence` values correspond to SMARTS `v` property. `valence` values of inorganic atoms would be `nothing`.
"""
@cache function valence(mol::GraphMol)
    vec = Union{Int,Nothing}[]
    atomsymbol_ = atomsymbol(mol)
    lonepair_ = lonepair(mol)
    for i in 1:nodecount(mol)
        if lonepair_[i] === nothing
            push!(vec, nothing)
        elseif atomsymbol_[i] === :H
            push!(vec, 1)
        else
            push!(vec, 4 - abs(lonepair_[i]))
        end
    end
    return vec
end

valence(view::SubgraphView) = valence(view.graph)


"""
    explicithcount(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of the number of adjacent explicit hydrogen nodes within the molecule that have ``n`` atoms.
"""
@cache function explicithcount(mol::GraphMol)
    vec = zeros(Int, nodecount(mol))
    atomsymbol_ = atomsymbol(mol)
    for i in 1:nodecount(mol)
        for adj in adjacencies(mol, i)
            if atomsymbol_[adj] == :H
                vec[i] += 1
            end
        end
    end
    return vec
end

explicithcount(view::SubgraphView) = explicithcount(view.graph)


"""
    implicithcount(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of the number of implicit hydrogens within
the molecule that have ``n`` atoms.
"""
@cache function implicithcount(mol::GraphMol)
    hcnt = (v, av) -> v === nothing ? 0 : max(0, v - av)
    return hcnt.(valence(mol), apparentvalence(mol))
end

implicithcount(view::SubgraphView) = implicithcount(view.graph)


"""
    heavyatomcount(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of the number of adjacent non-hydrogen atoms within
the molecule that have ``n`` atoms.
"""
heavyatomcount(mol::GraphMol) = nodedegree(mol) - explicithcount(mol)
heavyatomcount(view::SubgraphView) = heavyatomcount(view.graph)


"""
    hcount(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of the total number of hydrogens attached to the atom within the molecule that have ``n`` atoms.

`hcount` values correspond to SMARTS `H` property.
"""
@cache hcount(mol::GraphMol) = explicithcount(mol) + implicithcount(mol)
hcount(view::SubgraphView) = hcount(view.graph)


"""
    connectivity(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of the total number of adjacent atoms (including implicit hydrogens) within the molecule that have ``n`` atoms.

`connectivity` values correspond to SMARTS `X` property.
"""
@cache connectivity(mol::GraphMol) = nodedegree(mol) + implicithcount(mol)
connectivity(view::SubgraphView) = connectivity(view.graph)



# Composition

"""
    countatoms(mol::GraphMol) -> Dict{Symbol,Int}

Count the number of atoms and return symbol => count dict.
"""
function countatoms(mol::GraphMol)
    counter = Dict{Symbol,Int}()
    for sym in atomsymbol(mol)
        if !haskey(counter, sym)
            counter[sym] = 1
        else
            counter[sym] += 1
        end
    end
    hcnt = reduce(+, implicithcount(mol); init=0)
    if hcnt > 0
        if !haskey(counter, :H)
            counter[:H] = hcnt
        else
            counter[:H] += hcnt
        end
    end
    return counter
end


function writeformula(counter::Dict{Symbol,Int})
    strs = []
    if haskey(counter, :C)
        push!(strs, "C")
        c = pop!(counter, :C)
        c > 1 && push!(strs, string(c))
        if haskey(counter, :H)
            push!(strs, "H")
            h = pop!(counter, :H)
            h > 1 && push!(strs, string(h))
        end
    end
    for sym in sort(collect(keys(counter)))
        push!(strs, string(sym))
        cnt = counter[sym]
        cnt > 1 && push!(strs, string(cnt))
    end
    return join(strs)
end


"""
    molecularformula(mol::GraphMol) -> String

Return the molecular formula in Hill system.
"""
function molecularformula(mol::GraphMol)
    counter = countatoms(mol)
    return writeformula(counter)
end


"""
    empiricalformula(mol::GraphMol) -> String

Return the empirical formula in Hill system.
"""
function empiricalformula(mol::GraphMol)
    counter = countatoms(mol)
    if length(counter) > 1
        dv = reduce(gcd, values(counter))
        for k in keys(counter)
            counter[k] = div(counter[k], dv)
        end
    end
    return writeformula(counter)
end



# Hybridization

"""
    pielectron(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of the number of ``\\pi`` electrons (including implicit hydrogens) within the molecule that have ``n`` atoms.

The ``\\pi`` electron enumelation is based on the following simple rules.

- Any atom incident to a double bond -> +1
- Any atom incident to two double bond -> +2
- Any atom incident to a triple bond -> +1
- Uncharged N and O not incident to any multiple bonds and adjacent to atoms described above -> +2

These rules are applied for only typical organic atoms. The values for inorganic atoms will be 0.
"""
@cache function pielectron(mol::GraphMol)
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    pie_ = apparentvalence(mol) - nodedegree(mol)
    vec = zeros(Int, nodecount(mol))
    for i in 1:nodecount(mol)
        vec[i] = pie_[i]
        atomsymbol_[i] in (:N, :O) || continue
        for adj in adjacencies(mol, i)
            if pie_[i] == 0 && pie_[adj] > 0 && charge_[i] == 0
                vec[i] = 2
                break
            end
        end
    end
    return vec
end

pielectron(view::SubgraphView) = pielectron(view.graph)


"""
    hybridization(mol::GraphMol) -> Vector{Int}

Return the size ``n`` vector of orbital hybridization symbols (`:sp3`, `:sp2`, `:sp`, `:none`) within the molecule that have ``n`` atoms.

These hybridizations are for only typical organic atoms. Inorganic atoms and other orbitals like s, sp3d and sp3d2 will be `:none`.
"""
@cache function hybridization(mol::GraphMol)
    atomsymbol_ = atomsymbol(mol)
    connectivity_ = connectivity(mol)
    lonepair_ = lonepair(mol)
    pielectron_ = pielectron(mol)
    organic_symbols = [:B, :C, :N, :O, :F, :Si, :P, :S, :Cl, :As, :Se, :Br, :I]
    vec = Symbol[]
    for i in 1:nodecount(mol)
        if atomsymbol_[i] in organic_symbols
            val = connectivity_[i] + lonepair_[i]
            if val == 4
                if (atomsymbol_[i] in [:N, :O] && pielectron_[i] == 2)
                    push!(vec, :sp2)  # adjacent to conjugated bonds
                else
                    push!(vec, :sp3)
                end
            elseif val == 3
                push!(vec, :sp2)
            elseif val == 2
                push!(vec, :sp)
            else
                push!(vec, :none)
            end
        else
            push!(vec, :none)
        end
    end
    return vec
end

hybridization(view::SubgraphView) = hybridization(view.graph)



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


"""
    rotatablecount(mol::GraphMol) -> Int

Return the number of rotatable bonds.
"""
rotatablecount(mol::GraphMol) = reduce(+, isrotatable(mol); init=0)



# Aromatic rings

function satisfyhuckel(mol::GraphMol, ring)
    # TODO: need improvement
    cnt = 0
    atomsymbol_ = atomsymbol(mol)
    nodedegree_ = nodedegree(mol)
    pie_ = apparentvalence(mol) - nodedegree(mol)
    lonepair_ = lonepair(mol)
    carbonylO = findall(
        (atomsymbol_ .== :O) .* (nodedegree_ .== 1) .* (pie_ .== 1))
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
        elseif pie_[r] == 1
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



# Molecular properties

function stdweight(mol::GraphMol)
    weight = (atom, imh) -> atomweight(atom) + H_WEIGHT * imh
    return weight.(nodeattrs(mol), implicithcount(mol))
end

stdweight(view::SubgraphView) = stdweight(view.graph)


"""
    molweight(mol::GraphMol; digits=2) -> Float64

Return standard molecular weight.
"""
molweight(mol::GraphMol; digits=2
    ) = round(reduce(+, stdweight(mol); init=0), digits=digits)


function ishacceptor(mol::GraphMol)
    ac = (sym, lp) -> lp === nothing ? false : sym in (:N, :O, :F) && lp > 0
    return ac.(atomsymbol(mol), lonepair(mol))
end

ishacceptor(view::SubgraphView) = ishacceptor(view.graph)


"""
    hacceptorcount(mol::GraphMol) -> Int

Return the number of hydrogen bond acceptors (N, O and F).
"""
hacceptorcount(mol::GraphMol) = reduce(+, ishacceptor(mol); init=0)


function ishdonor(mol::GraphMol)
    dc = (sym, h) -> sym in (:N, :O) && h > 0
    return dc.(atomsymbol(mol), hcount(mol))
end

ishdonor(view::SubgraphView) = ishdonor(view.graph)


"""
    hdonorcount(mol::GraphMol) -> Int

Return the number of hydrogen bond donors (O and N attached to hydrogens).
"""
hdonorcount(mol::GraphMol) = reduce(+, ishdonor(mol); init=0)
