#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    nodedegree,
    sssr, sssrmembership, sssrbondmembership,
    fusedrings, fusedringmembership, 
    smallestsssr, sssrcount,
    isringatom, isringbond,
    atomsymbol, charge, multiplicity, bondorder,
    valence, lonepair,
    heavyatomconnected, explicithconnected, implicithconnected,
    hydrogenconnected, connectivity,
    ishdonor, hacceptorcount, ishacceptor, hdonorcount,
    isrotatable, rotatablecount,
    atomcounter, heavyatomcount, molecularformula, empiricalformula,
    pielectron, hybridization,
    isaromaticring, isaromatic, isaromaticbond,
    precalculate!


const LONEPAIR_COUNT = Dict(
    :H => 0, :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
    :Si => 0, :P => 1, :S => 2, :Cl => 3,
    :As => 1, :Se => 2, :Br => 3, :I => 3
)



# Molecular graph topology

"""
    nodedegree(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing the node degree of the molecular graph
of 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `D` query.
"""
@cachefirst nodedegree(mol::GraphMol) = [degree(mol, n) for n in 1:nodecount(mol)]
nodedegree(view::SubgraphView) = nodedegree(view.graph)


"""
    sssr(mol::GraphMol) -> Vector{Vector{Int}}

Return vectors of ring nodes representing small set of smallest rings (SSSR).

See [`Graph.minimumcyclebasis`](@ref).
"""
sssr = minimumcyclebasisnodes


"""
    sssrmembership(mol::GraphMol) -> Vector{Set{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th atoms of the given molecule.

SSSR membership is represented as a set of SSSR indices assigned to each rings.
This means atoms that have the same SSSR index belong to the same SSSR.
"""
@cachefirst function sssrmembership(mol::GraphMol)
    nodes = [Set{Int}() for n in 1:nodecount(mol)]
    for (i, cyc) in enumerate(sssr(mol))
        for n in cyc
            push!(nodes[n], i)
        end
    end
    return nodes
end


"""
    sssrbondmembership(mol::GraphMol) -> Vector{Set{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th bonds of the given molecule.

SSSR membership is represented as a set of SSSR indices assigned to each rings.
This means bonds that have the same SSSR index belong to the same SSSR.
"""
@cachefirst function sssrbondmembership(mol::GraphMol)
    edges = [Set{Int}() for n in 1:edgecount(mol)]
    for (i, cyc) in enumerate(minimumcyclebasis(mol))
        for e in cyc
            push!(edges[e], i)
        end
    end
    return edges
end


"""
    fusedrings(mol::UndirectedGraph) -> Vector{Set{Int}}

Return vectors of fused ring node sets.

A fused ring is defined as a 2-edge connected components in terms of graph theory.
Spirocyclic structures are considered to be part of a fused ring.
"""
@cachefirst function fusedrings(mol::GraphMol)
    cobr = setdiff(edgeset(mol), bridges(mol))
    subg = plaingraph(edgesubgraph(mol, cobr))
    return connectedcomponents(subg)
end


"""
    fusedringmembership(mol::UndirectedGraph) -> Vector{Int}

Return a vector of size ``n`` representing [`fusedrings`](@ref) membership of
1 to ``n``th atoms of the given molecule.

Fused ring membership is represented as a set of fused ring indices assigned to each fused rings.
This means atoms that have the same fused ring index belong to the same fused ring.
"""
@cachefirst function fusedringmembership(mol::GraphMol)
    arr = [Set{Int}() for i in 1:nodecount(mol)]
    for (i, conn) in enumerate(fusedrings(mol))
        for n in conn
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    smallestsssr(mol::UndirectedGraph) -> Vector{Int}

Return a vector of size ``n`` representing the size of the smallest [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to.

This property corresponds to SMARTS `r` query.
"""
@cachefirst function smallestsssr(mol::GraphMol)
    sssr_ = sssr(mol)
    arr = []
    for ks in sssrmembership(mol)
        val = isempty(ks) ? 0 : minimum([length(sssr_[k]) for k in ks])
        push!(arr, val)
    end
    return arr
end
smallestsssr(view::SubgraphView) = smallestsssr(view.graph)


"""
    sssrcount(mol::UndirectedGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to.

This property corresponds to SMARTS `R` query.
"""
@cachefirst sssrcount(mol::UndirectedGraph) = length.(sssrmembership(mol))
sssrcount(view::SubgraphView) = sssrcount(view.graph)


"""
    isringatom(mol::UndirectedGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th atoms of
the given molecule belong to a ring or not.
"""
@cachefirst isringatom(mol::UndirectedGraph) = .!isempty.(sssrmembership(mol))
isringatom(view::SubgraphView) = isringatom(view.graph)


"""
    isringbond(mol::UndirectedGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th bonds of
the given molecule belong to a ring or not.
"""
@cachefirst isringbond(mol::UndirectedGraph) = .!isempty.(sssrbondmembership(mol))
isringbond(view::SubgraphView) = isringbond(view.graph)




# Elemental properties

"""
    atomsymbol(mol::GraphMol) -> Vector{Symbol}

Return a vector of size ``n`` representing atom symbols of 1 to ``n``th atoms of
the given molecule.
"""
@cachefirst atomsymbol(mol::UndirectedGraph) = getproperty.(nodeattrs(mol), :symbol)


"""
    charge(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing atom charges of 1 to ``n``th atoms of
the given molecule.
"""
@cachefirst charge(mol::UndirectedGraph) = getproperty.(nodeattrs(mol), :charge)


"""
    multiplicity(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing atom multiplicities of 1 to ``n``th atoms of
the given molecule (1: non-radical, 2: radical, 3: biradical).
"""
@cachefirst multiplicity(mol::UndirectedGraph) = getproperty.(nodeattrs(mol), :multiplicity)


"""
    bondorder(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing bond order of 1 to ``n``th bonds of
the given molecule.
"""
@cachefirst bondorder(mol::UndirectedGraph) = getproperty.(edgeattrs(mol), :order)




# Valence

"""
    lonepair(mol::GraphMol) -> Vector{Union{Int,Nothing}}

Return a vector of size ``n`` representing the number of lone pairs of
1 to ``n``th atoms of the given molecule.

The number of lone pair in inorganic atoms would be `nothing`.
The result can take negative value if the atom has empty shells (e.g. B).
"""
@cachefirst function lonepair(mol::GraphMol)
    vec = Union{Int,Nothing}[]
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    for i in 1:nodecount(mol)
        num = get(LONEPAIR_COUNT, atomsymbol_[i], nothing)
        v = num === nothing ? nothing : num - charge_[i]
        push!(vec, v)
    end
    return vec
end
lonepair(view::SubgraphView) = lonepair(view.graph)


"""
    apparentvalence(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing the number of total bond order
incident to 1 to ``n``th atoms of the given molecule.
"""
@cachefirst function apparentvalence(mol::GraphMol)
    vec = zeros(Int, nodecount(mol))
    bondorder_ = bondorder(mol)
    for i in 1:nodecount(mol)
        for (inc, adj) in neighbors(mol, i)
            vec[i] += bondorder_[inc]
        end
    end
    return vec
end

apparentvalence(view::SubgraphView) = apparentvalence(view.graph)


"""
    valence(mol::GraphMol) -> Vector{Union{Int,Nothing}}

Return a vector of size ``n`` representing the intrinsic valence of
1 to ``n``th atoms of the given molecule.

The number of implicit hydrogens would be calculated based on the valence.
The valence value in inorganic atoms would be `nothing`.
This property corresponds to SMARTS `v` query.
"""
@cachefirst function valence(mol::GraphMol)
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
    explicithconnected(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing the number of explicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Explicit" means hydrogens are connected to the heavy atom as atom nodes.
"""
@cachefirst function explicithconnected(mol::GraphMol)
    vec = zeros(Int, nodecount(mol))
    atomsymbol_ = atomsymbol(mol)
    for i in 1:nodecount(mol)
        for (inc, adj) in neighbors(mol, i)
            if atomsymbol_[adj] === :H
                vec[i] += 1
            end
        end
    end
    return vec
end
explicithconnected(view::SubgraphView) = explicithconnected(view.graph)


"""
    implicithconnected(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing the number of implicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Implicit" means hydrogens are not represented as graph nodes,
but it can be infered from the intrinsic valence of typical organic atoms.
"""
@cachefirst function implicithconnected(mol::UndirectedGraph)
    hcnt = (v, av) -> v === nothing ? 0 : max(0, v - av)
    return hcnt.(valence(mol), apparentvalence(mol))
end


"""
    heavyatomconnected(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing the number of non-hydrogen atoms
connected to 1 to ``n``th atoms of the given molecule.
"""
@cachefirst heavyatomconnected(mol::UndirectedGraph) = nodedegree(mol) - explicithconnected(mol)


"""
    hydrogenconnected(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing the number of total hydrogens
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `H` query.
"""
@cachefirst hydrogenconnected(mol::UndirectedGraph) = explicithconnected(mol) + implicithconnected(mol)


"""
    connectivity(mol::GraphMol) -> Vector{Int}

Return a vector of size ``n`` representing the number of total atoms
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `X` query.
"""
@cachefirst connectivity(mol::UndirectedGraph) = nodedegree(mol) + implicithconnected(mol)




# Hydrogen bond donor/acceptor

@cachefirst function ishacceptor(mol::UndirectedGraph)
    ac = (sym, lp) -> lp === nothing ? false : sym in (:N, :O, :F) && lp > 0
    return ac.(atomsymbol(mol), lonepair(mol))
end


"""
    hacceptorcount(mol::GraphMol) -> Int

Return the total number of hydrogen bond acceptors (N, O and F).
"""
hacceptorcount(mol::GraphMol) = reduce(+, ishacceptor(mol); init=0)


@cachefirst function ishdonor(mol::UndirectedGraph)
    dc = (sym, h) -> sym in (:N, :O) && h > 0
    return dc.(atomsymbol(mol), hydrogenconnected(mol))
end


"""
    hdonorcount(mol::GraphMol) -> Int

Return the total number of hydrogen bond donors (O and N attached to hydrogens).
"""
hdonorcount(mol::GraphMol) = reduce(+, ishdonor(mol); init=0)




# Rotatable bonds

"""
    isrotatable(mol::GraphMol)

Return a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule are rotatable or not.
"""
@cachefirst function isrotatable(mol::GraphMol)
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

Return the total number of rotatable bonds.
"""
rotatablecount(mol::GraphMol) = reduce(+, isrotatable(mol); init=0)




# Composition

"""
    atomcounter(mol::GraphMol) -> Dict{Symbol,Int}

Count the number of atoms and return symbol => count dict.
"""
@cachefirst function atomcounter(mol::GraphMol)
    counter = Dict{Symbol,Int}()
    for sym in atomsymbol(mol)
        if !haskey(counter, sym)
            counter[sym] = 1
        else
            counter[sym] += 1
        end
    end
    hcnt = reduce(+, implicithconnected(mol); init=0)
    if hcnt > 0
        if !haskey(counter, :H)
            counter[:H] = hcnt
        else
            counter[:H] += hcnt
        end
    end
    return counter
end


"""
    heavyatomcount(mol::GraphMol) -> Int

Return the total number of non-hydrogen atoms.
"""
function heavyatomcount(mol::GraphMol)
    counter = atomcounter(mol)
    delete!(counter, :H)
    return reduce(+, values(counter); init=0)
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
    counter = atomcounter(mol)
    return writeformula(counter)
end


"""
    empiricalformula(mol::GraphMol) -> String

Return the empirical formula in Hill system.
"""
function empiricalformula(mol::GraphMol)
    counter = atomcounter(mol)
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

Returns a vector of size ``n`` representing the number of ``\\pi`` electrons
of 1 to ``n``th atoms of the given molecule.

The counting of ``\\pi`` electrons is based on the following rules.

- Any atom incident to a double bond -> +1
- Any atom incident to two double bond -> +2
- Any atom incident to a triple bond -> +2
- Any other uncharged N, O or S that are neighbor of multiple bonds -> +2

These rules are applied for only typical organic atoms. The values for inorganic atoms will be 0.
"""
@cachefirst function pielectron(mol::GraphMol)
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    pie_ = apparentvalence(mol) - nodedegree(mol)
    vec = zeros(Int, nodecount(mol))
    for i in 1:nodecount(mol)
        vec[i] = pie_[i]
        atomsymbol_[i] in (:N, :O, :S) || continue
        for (inc, adj) in neighbors(mol, i)
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

Returns a vector of size ``n`` representing the orbital hybridization symbols
(`:sp3`, `:sp2`, `:sp` or `:none`) of 1 to ``n``th atoms of the given molecule.

The hybridization value in inorganic atoms and non-typical organic atoms will be `:none`
(e.g. s, sp3d and sp3d2 orbitals).
"""
@cachefirst function hybridization(mol::GraphMol)
    atomsymbol_ = atomsymbol(mol)
    pielectron_ = pielectron(mol)
    connectivity_ = connectivity(mol)
    lonepair_ = lonepair(mol)
    vec = Symbol[]
    for i in 1:nodecount(mol)
        if lonepair_[i] === nothing || atomsymbol_[i] === :H
            push!(vec, :none)
            continue
        end
        orbitals = connectivity_[i] + lonepair_[i]
        if orbitals == 4
            if (atomsymbol_[i] in [:N, :O] && pielectron_[i] == 2)
                push!(vec, :sp2)  # adjacent to conjugated bonds
            else
                push!(vec, :sp3)
            end
        elseif orbitals == 3
            push!(vec, :sp2)
        elseif orbitals == 2
            push!(vec, :sp)
        else
            push!(vec, :none)
        end
    end
    return vec
end
hybridization(view::SubgraphView) = hybridization(view.graph)




# Aromaticity

@cachefirst function isaromaticring(mol::GraphMol)
    atomsymbol_ = atomsymbol(mol)
    nodedegree_ = nodedegree(mol)
    pie_ = apparentvalence(mol) - nodedegree(mol)
    lonepair_ = lonepair(mol)
    carbonylO = findall(
        (atomsymbol_ .=== :O) .* (nodedegree_ .== 1) .* (pie_ .== 1))
    carbonylC = falses(nodecount(mol))
    for o in carbonylO
        c = pop!(adjacencies(mol, o))
        if atomsymbol_[c] === :C
            carbonylC[c] = true
        end
    end
    arr = falses(circuitrank(mol))
    for (i, ring) in enumerate(sssr(mol))
        cnt = 0
        for r in ring
            if carbonylC[r]
                continue
            elseif pie_[r] == 1
                cnt += 1
            elseif lonepair_[r] === nothing
                cnt = 0
                break
            elseif lonepair_[r] > 0
                cnt += 2
            elseif lonepair_[r] < 0
                continue
            else
                cnt = 0
                break
            end
        end
        if cnt % 4 == 2
            # Huckel rule check
            arr[i] = true
            continue
        end
        if nodeattrtype(mol) === SmilesAtom
            # SMILES aromatic atom
            sub = nodesubgraph(mol, Set(ring))
            if all(nodeattr(mol, n).isaromatic for n in nodeset(sub))
                arr[i] = true
            end
        end
    end
    return arr
end


"""
    isaromatic(mol::GraphMol) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th atoms
of the given molecule belong to an aromatic ring or not.

Some kind of aromaticity resulting from long conjugated chains and charge
delocalization may be unrecognizable. Also, non-classical aromaticity
such as Moebius aromaticity is not considered.
"""
@cachefirst function isaromatic(mol::GraphMol)
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


"""
    isaromaticbond(mol::GraphMol) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule belong to an aromatic ring or not.

See [`isaromatic`](@ref).
"""
@cachefirst function isaromaticbond(mol::GraphMol)
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



"""
    precalculate!(mol::GraphMol)

Convenient method to pre-calculate and cache performance bottleneck descriptors.
"""
function precalculate!(mol)
    setcache!(mol, :minimumcyclebasis)
    setcache!(mol, :sssr)
    setcache!(mol, :lonepair)
    setcache!(mol, :apparentvalence)
    setcache!(mol, :valence)
    setcache!(mol, :isaromaticring)
    nodeattrtype(mol) === SmilesAtom && setcache!(mol, :coordgen)
end
