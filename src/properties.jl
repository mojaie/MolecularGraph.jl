#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    sssr, which_ring, edge_which_ring, fused_rings, which_fused_ring, 
    smallest_ring, ring_count, is_in_ring, is_edge_in_ring,
    atom_symbol, charge, multiplicity, bond_order,
    valence, lone_pair,
    heavy_atoms, explicit_hydrogens, implicit_hydrogens,
    total_hydrogens, connectivity,
    is_hydrogen_donor, hydrogen_donor_count,
    is_hydrogen_acceptor, hydrogen_acceptor_count,
    is_rotatable, rotatable_count,
    atom_counter, heavy_atom_count, molecular_formula, empirical_formula,
    pi_electron, hybridization,
    is_ring_aromatic, is_aromatic, is_edge_aromatic,
    precalculate!


const LONEPAIR_COUNT = Dict(
    :H => 0, :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
    :Si => 0, :P => 1, :S => 2, :Cl => 3,
    :As => 1, :Se => 2, :Br => 3, :I => 3
)



# Molecular graph topology

"""
    degree(mol::AbstractMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the node degree of the molecular graph
of 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `D` query.
"""
Graphs.degree


"""
    sssr(mol::AbstractMolGraph) -> Vector{Vector{Int}}

Return vectors of ring nodes representing small set of smallest rings (SSSR).

See [`Graph.minimumcyclebasis`](@ref).
"""
sssr(mol::AbstractMolGraph) = get(mol.gprops, :sssr, mincyclebasis(mol.graph))


"""
    which_ring(mol::MolGraph) -> Vector{Set{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th nodes of the given graph.

SSSR membership is represented as a vector of SSSR indices assigned to each rings.
This means nodes that have the same SSSR index belong to the same SSSR.
"""
function which_ring(mol::MolGraph)
    nodes = [Int[] for _ in vertices(mol)]
    for (i, cyc) in enumerate(sssr(mol))
        for n in cyc
            push!(nodes[n], i)
        end
    end
    return nodes
end


"""
    edge_which_sssr(mol::MolGraph) -> Vector{Set{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th bonds of the given molecule.

SSSR membership is represented as a set of SSSR indices assigned to each rings.
This means bonds that have the same SSSR index belong to the same SSSR.
"""
function edge_which_ring(mol::MolGraph)
    edges = [Int[] for n in 1:ne(mol)]
    for (i, cyc) in enumerate(sssr(mol))
        for j in 1:(length(cyc) - 1)
            push!(edges[edge_rank(mol, cyc[j], cyc[j + 1])], i)
        end
        push!(edges[edge_rank(mol, cyc[1], cyc[end])], i)
    end
    return edges
end


"""
    fused_rings(mol::MolGraph{T,V,E}) -> Vector{Vector{T}}

Return vectors of fused ring node sets.

A fused ring is defined as a 2-edge connected components in terms of graph theory.
Spirocyclic structures are considered to be part of a fused ring.
"""
function fused_rings(mol::MolGraph)
    has_prop(mol, :fused_rings) && return get_prop(mol, :fused_rings)
    cobr = setdiff(Set(edges(mol)), bridges(mol.graph))
    subg, vmap = induced_subgraph(mol.graph, cobr)
    return  [vmap[c] for c in connected_components(subg)]
end


"""
    which_fused_ring(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing [`fusedrings`](@ref) membership of
1 to ``n``th atoms of the given molecule.

Fused ring membership is represented as a set of fused ring indices assigned to each fused rings.
This means atoms that have the same fused ring index belong to the same fused ring.
"""
function which_fused_ring(mol::MolGraph)
    has_prop(mol, :node_which_fused_ring) && return get_prop(mol, :node_which_fused_ring)
    arr = [T[] for _ in vertices(mol)]
    for (i, conn) in enumerate(fused_rings(mol))
        for n in conn
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    smallest_ring(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the size of the smallest [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to. 

If the node is not in a ring, the value would be 0.
This property corresponds to SMARTS `r` query.
"""
function smallest_ring(mol::MolGraph)
    sssr_ = sssr(mol)
    return [(isempty(ks) ? 0 : minimum(length.(sssr_[ks]))) for ks in which_ring(mol)]
end


"""
    ring_count(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to.

This property corresponds to SMARTS `R` query.
"""
ring_count(mol::MolGraph) = length.(which_ring(mol))


"""
    is_in_ring(mol::MolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th atoms of
the given molecule belong to a ring or not.
"""
is_in_ring(mol::MolGraph) = .!isempty.(which_ring(mol))


"""
    is_edge_in_ring(mol::MolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th bonds of
the given molecule belong to a ring or not.
"""
is_edge_in_ring(mol::MolGraph) = .!isempty.(edge_which_ring(mol))




# Elemental properties

"""
    atom_symbol(mol::MolGraph) -> Vector{Symbol}

Return a vector of size ``n`` representing atom symbols of 1 to ``n``th atoms of
the given molecule.
"""
atom_symbol(mol::MolGraph) = getproperty.(mol.vprops, :symbol)


"""
    charge(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom charges of 1 to ``n``th atoms of
the given molecule.
"""
charge(mol::MolGraph) = getproperty.(mol.vprops, :charge)


"""
    multiplicity(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom multiplicities of 1 to ``n``th atoms of
the given molecule (1: non-radical, 2: radical, 3: biradical).
"""
multiplicity(mol::MolGraph) = getproperty.(mol.vprops, :multiplicity)


"""
    bond_order(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing bond order of 1 to ``n``th bonds of
the given molecule.
"""
bond_order(mol::MolGraph) = getproperty.(mol.eprops, :order)
bond_order_map(mol::MolGraph) = Dict(k => v[:order] for (k, v) in zip(edges(mol), mol.eprops))



# Valence

"""
    lone_pair(mol::MolGraph) -> Vector{Union{Int,Nothing}}

Return a vector of size ``n`` representing the number of lone pairs of
1 to ``n``th atoms of the given molecule.

The number of lone pair in inorganic atoms would be `nothing`.
The result can take negative value if the atom has empty shells (e.g. B).
"""
function lone_pair(mol::MolGraph)
    has_prop(mol, :node_lone_pair) && return get_prop(mol, :node_lone_pair)
    vec = Union{Int,Nothing}[]
    atomsymbol_ = atom_symbol(mol)
    charge_ = charge(mol)
    for i in vertices(mol)
        num = get(LONEPAIR_COUNT, atomsymbol_[i], nothing)
        v = num === nothing ? nothing : num - charge_[i]
        push!(vec, v)
    end
    return vec
end


"""
    apparent_valence(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total bond order
incident to 1 to ``n``th atoms of the given molecule.
"""
function apparent_valence(mol::MolGraph{T,V,E}) where {T,V,E}
    has_prop(mol, :node_apparent_valence) && return get_prop(mol, :node_apparent_valence)
    vec = zeros(Int, nv(mol))
    bomap = bond_order_map(mol)
    for e in edges(mol)
        vec[src(e)] += bomap[e]
        vec[dst(e)] += bomap[e]
    end
    return vec
end


"""
    valence(mol::MolGraph) -> Vector{Union{Int,Nothing}}

Return a vector of size ``n`` representing the intrinsic valence of
1 to ``n``th atoms of the given molecule.

The number of implicit hydrogens would be calculated based on the valence.
The valence value in inorganic atoms would be `nothing`.
This property corresponds to SMARTS `v` query.
"""
function valence(mol::MolGraph)
    has_prop(mol, :node_valence) && return get_prop(mol, :node_valence)
    vec = Union{Int,Nothing}[]
    atomsymbol_ = atom_symbol(mol)
    lonepair_ = lone_pair(mol)
    for i in vertices(mol)
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


"""
    explicit_hydrogens(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of explicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Explicit" means hydrogens are explicitly represented as graph nodes.
"""
function explicit_hydrogens(mol::MolGraph)
    vec = zeros(Int, nv(mol))
    for (i, sym) in enumerate(atom_symbol(mol))
        sym === :H || continue
        for nbr in neighbors(mol, i)
            vec[nbr] += 1
        end
    end
    return vec
end


"""
    implicit_hydrogens(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of implicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Implicit" means hydrogens are not represented as graph nodes,
but it can be infered from the intrinsic valence of typical organic atoms.
"""
function implicit_hydrogens(mol::MolGraph)
    hcnt = (v, av) -> v === nothing ? 0 : max(0, v - av)
    return hcnt.(valence(mol), apparent_valence(mol))
end


"""
    heavyatoms(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of non-hydrogen atoms
connected to 1 to ``n``th atoms of the given molecule.
"""
heavy_atoms(mol::MolGraph) = degree(mol) - explicit_hydrogens(mol)


"""
    hydrogenconnected(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total hydrogens
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `H` query.
"""
total_hydrogens(mol::MolGraph) = explicit_hydrogens(mol) + implicit_hydrogens(mol)


"""
    connectivity(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total atoms
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `X` query.
"""
connectivity(mol::MolGraph) = degree(mol) + implicit_hydrogens(mol)




# Hydrogen bond donor/acceptor

function is_hydrogen_acceptor(mol::MolGraph)
    ac = (sym, lp) -> lp === nothing ? false : sym in (:N, :O, :F) && lp > 0
    return ac.(atom_symbol(mol), lone_pair(mol))
end


"""
    hacceptorcount(mol::MolGraph) -> Int

Return the total number of hydrogen bond acceptors (N, O and F).
"""
hydrogen_acceptor_count(mol::MolGraph) = reduce(+, is_hydrogen_acceptor(mol); init=0)


function is_hydrogen_donor(mol::MolGraph)
    dc = (sym, h) -> sym in (:N, :O) && h > 0
    return dc.(atom_symbol(mol), total_hydrogens(mol))
end


"""
    hdonorcount(mol::MolGraph) -> Int

Return the total number of hydrogen bond donors (O and N attached to hydrogens).
"""
hydrogen_donor_count(mol::MolGraph) = reduce(+, is_hydrogen_donor(mol); init=0)




# Rotatable bonds

"""
    isrotatable(mol::MolGraph)

Return a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule are rotatable or not.
"""
function is_rotatable(mol::MolGraph)
    degree_ = degree(mol)
    edgeinring_ = is_edge_in_ring(mol)
    bondorder_ = bond_order(mol)
    vec = Bool[]
    for (i, e) in enumerate(edges(mol))
        rot = (!edgeinring_[i] && bondorder_[i] == 1
            && degree_[src(e)] != 1 && degree_[dst(e)] != 1)
        push!(vec, rot)
    end
    return vec
end


"""
    rotatablecount(mol::MolGraph) -> Int

Return the total number of rotatable bonds.
"""
rotatable_count(mol::MolGraph) = reduce(+, is_rotatable(mol); init=0)




# Composition

"""
    atomcounter(mol::MolGraph) -> Dict{Symbol,Int}

Count the number of atoms and return symbol => count dict.
"""
function atom_counter(mol::MolGraph)
    counter = Dict{Symbol,Int}()
    for sym in atom_symbol(mol)
        if !haskey(counter, sym)
            counter[sym] = 1
        else
            counter[sym] += 1
        end
    end
    hcnt = reduce(+, implicit_hydrogens(mol); init=0)
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
    heavyatomcount(mol::MolGraph) -> Int

Return the total number of non-hydrogen atoms.
"""
function heavy_atom_count(mol::MolGraph)
    counter = atom_counter(mol)
    delete!(counter, :H)
    return reduce(+, values(counter); init=0)
end


function write_formula(counter::Dict{Symbol,Int})
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
    molecular_formula(mol::MolGraph) -> String

Return the molecular formula in Hill system.
"""
function molecular_formula(mol::MolGraph)
    counter = atom_counter(mol)
    return write_formula(counter)
end


"""
    empirical_formula(mol::MolGraph) -> String

Return the empirical formula in Hill system.
"""
function empirical_formula(mol::MolGraph)
    counter = atom_counter(mol)
    if length(counter) > 1
        dv = reduce(gcd, values(counter))
        for k in keys(counter)
            counter[k] = div(counter[k], dv)
        end
    end
    return write_formula(counter)
end




# Hybridization

"""
    pi_electron(mol::MolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the number of ``\\pi`` electrons
of 1 to ``n``th atoms of the given molecule.

The counting of ``\\pi`` electrons is based on the following rules.

- Any atom incident to a double bond -> +1
- Any atom incident to two double bond -> +2
- Any atom incident to a triple bond -> +2
- Any other uncharged N, O or S that are neighbor of multiple bonds -> +2

These rules are applied for only typical organic atoms. The values for inorganic atoms will be 0.
"""
function pi_electron(mol::MolGraph)
    atomsymbol_ = atom_symbol(mol)
    charge_ = charge(mol)
    pie_ = apparent_valence(mol) - degree(mol)
    vec = zeros(Int, nv(mol))
    for i in vertices(mol)
        vec[i] = pie_[i]
        atomsymbol_[i] in (:N, :O, :S) || continue
        for nbr in neighbors(mol, i)
            if pie_[i] == 0 && pie_[nbr] > 0 && charge_[i] == 0
                vec[i] = 2
                break
            end
        end
    end
    return vec
end


"""
    hybridization(mol::MolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the orbital hybridization symbols
(`:sp3`, `:sp2`, `:sp` or `:none`) of 1 to ``n``th atoms of the given molecule.

The hybridization value in inorganic atoms and non-typical organic atoms will be `:none`
(e.g. s, sp3d and sp3d2 orbitals).
"""
function hybridization(mol::MolGraph)
    atomsymbol_ = atom_symbol(mol)
    pielectron_ = pi_electron(mol)
    connectivity_ = connectivity(mol)
    lonepair_ = lone_pair(mol)
    vec = Symbol[]
    for i in vertices(mol)
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




# Aromaticity

function is_ring_aromatic(mol::MolGraph)
    has_prop(mol, :is_ring_aromatic) && return get_prop(mol, :is_ring_aromatic)
    atomsymbol_ = atom_symbol(mol)
    nodedegree_ = degree(mol)
    pie_ = apparent_valence(mol) - degree(mol)
    lonepair_ = lone_pair(mol)
    carbonyl_o = findall(
        (atomsymbol_ .=== :O) .* (nodedegree_ .== 1) .* (pie_ .== 1))
    carbonyl_c = falses(nv(mol))
    for o in carbonyl_o
        c = neighbors(mol, o)[1]
        if atomsymbol_[c] === :C
            carbonyl_c[c] = true
        end
    end
    sssr_ = sssr(mol)
    arr = falses(length(sssr_))
    for (i, ring) in enumerate(sssr_)
        cnt = 0
        for r in ring
            carbonyl_c[r] && continue
            if pie_[r] == 1
                cnt += 1
                continue
            elseif lonepair_[r] !== nothing
                if lonepair_[r] > 0
                    cnt += 2
                    continue
                elseif lonepair_[r] < 0
                    continue  # :B
                end  # lonepair == 0
            end
            cnt = 0
            break
        end
        if cnt % 4 == 2  # Huckel rule check
            arr[i] = true
        elseif vproptype(mol) <: SMILESAtom  # SMILES aromatic atom
            subg, vmap = induced_subgraph(mol.graph, ring)
            if all(get_prop(mol, vmap[n], :isaromatic) for n in vertices(subg))
                arr[i] = true
            end
        end
    end
    return arr
end


"""
    is_aromatic(mol::MolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th atoms
of the given molecule belong to an aromatic ring or not.

Some kind of aromaticity resulting from long conjugated chains and charge
delocalization may be unrecognizable. Also, non-classical aromaticity
such as Moebius aromaticity is not considered.
"""
function is_aromatic(mol::MolGraph)
    aromatic = falses(nv(mol))
    for ring in sssr(mol)[findall(is_ring_aromatic(mol))]
        aromatic[ring] .= true
    end
    return aromatic
end


"""
    is_edge_aromatic(mol::MolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule belong to an aromatic ring or not.

See [`isaromatic`](@ref).
"""
function is_edge_aromatic(mol::MolGraph)
    aromatic = falses(ne(mol))
    for ring in sssr(mol)[findall(is_ring_aromatic(mol))]
        for i in 1:(length(ring) - 1)
            aromatic[edge_rank(mol, p[i], p[i + 1])] = true
        end
        aromatic[edge_rank(mol, p[1], p[end])] = true
    end
    return aromatic
end



"""
    precalculate!(mol::MolGraph)

Convenient method to pre-calculate and cache performance bottleneck descriptors.
"""
function precalculate!(mol)
    mol.gprops[:sssr] = sssr(mol)
    mol.gprops[:node_lone_pair] = lone_pair(mol)
    mol.gprops[:node_apparent_valence] = apparent_valence(mol)
    mol.gprops[:node_valence] = valence(mol)
    mol.gprops[:is_ring_aromatic] = is_ring_aromatic(mol)
    # if vproptype(mol) <: SMILESAtom
    #     mol.gprops[:coordgen] = coordgen(mol)
end


# deprecated function names

nodedegree = degree
sssrmembership = which_ring
sssrbondmembership = edge_which_ring
fusedrings = fused_rings
fusedringmembership = which_fused_ring
smallestsssr = smallest_ring
sssrcount = ring_count
isringatom = is_in_ring
isringbond = is_edge_in_ring
atomsymbol(mol::MolGraph) = atom_symbol(mol)
bondorder = bond_order
lonepair = lone_pair
heavyatomconnected = heavy_atoms
explicithconnected = explicit_hydrogens
implicithconnected = implicit_hydrogens
hydrogenconnected = total_hydrogens
ishdonor = is_hydrogen_donor
hdonorcount = hydrogen_donor_count
ishacceptor = is_hydrogen_acceptor
hacceptorcount = hydrogen_acceptor_count
isrotatable = is_rotatable
rotatablecount = rotatable_count
atomcounter = atom_counter
heavyatomcount = heavy_atom_count
molecularformula = molecular_formula
empiricalformula = empirical_formula
pielectron = pi_electron
isaromaticring = is_ring_aromatic
isaromatic = is_aromatic
isaromaticbond = is_edge_aromatic