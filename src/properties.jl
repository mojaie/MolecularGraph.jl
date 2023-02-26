#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    sssr, sssr!,
    which_ring, edge_which_ring, fused_rings, which_fused_ring,
    smallest_ring, ring_count, is_in_ring, is_edge_in_ring,
    atom_symbol, charge, multiplicity, bond_order,
    lone_pair, lone_pair!, apparent_valence, apparent_valence!, valence, valence!,
    explicit_hydrogens, implicit_hydrogens, heavy_atoms,
    total_hydrogens, connectivity,
    is_hydrogen_donor, hydrogen_donor_count,
    is_hydrogen_acceptor, hydrogen_acceptor_count,
    is_rotatable, rotatable_count,
    atom_counter, heavy_atom_count, molecular_formula, empirical_formula,
    pi_electron, hybridization,
    is_ring_aromatic, is_ring_aromatic!, is_aromatic, is_edge_aromatic,
    precalculate!


const LONEPAIR_COUNT = Dict(
    :H => 0, :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
    :Si => 0, :P => 1, :S => 2, :Cl => 3,
    :As => 1, :Se => 2, :Br => 3, :I => 3
)



# Molecular graph topology descriptors

"""
    degree(mol::SimpleMolGraph{T}) -> Vector{T}

Return a vector of size ``n`` representing the node degree of the molecular graph
of 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `D` query.
"""
Graphs.degree(mol::SimpleMolGraph) = get(descriptors(mol), :v_degree, degree(mol.graph))


"""
    sssr(mol::SimpleMolGraph{T}) -> Vector{Vector{T}}

Return vectors of ring nodes representing small set of smallest rings (SSSR).

See [`Graph.minimumcyclebasis`](@ref).
"""
sssr(mol::SimpleMolGraph) = get(descriptors(mol), :sssr, mincyclebasis(mol.graph))
sssr!(mol::MolGraph) = set_descriptor!(mol, :sssr, mincyclebasis(mol.graph))


"""
    which_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th nodes of the given graph.

SSSR membership is represented as a vector of SSSR indices assigned to each rings.
This means nodes that have the same SSSR index belong to the same SSSR.
"""
function which_ring(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    has_descriptor(mol, :v_which_ring) && return get_descriptor(mol, :v_which_ring)
    arr = init_node_descriptor(Vector{Int}, mol)
    for (i, cyc) in enumerate(sssr(mol))
        for n in cyc
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    edge_which_sssr(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th bonds of the given molecule.

SSSR membership is represented as a set of SSSR indices assigned to each rings.
This means bonds that have the same SSSR index belong to the same SSSR.
"""
function edge_which_ring(mol::SimpleMolGraph)
    has_descriptor(mol, :e_which_ring) && return get_descriptor(mol, :e_which_ring)
    arr = init_edge_descriptor(Vector{Int}, mol)
    for (i, cyc) in enumerate(sssr(mol))
        for j in 1:(length(cyc) - 1)
            push!(arr[edge_rank(mol, cyc[j], cyc[j + 1])], i)
        end
        push!(arr[edge_rank(mol, cyc[1], cyc[end])], i)
    end
    return arr
end


"""
    fused_rings(mol::SimpleMolGraph{T}) -> Vector{Vector{T}}

Return vectors of fused ring node sets.

A fused ring is defined as a 2-edge connected components in terms of graph theory.
Spirocyclic structures are considered to be part of a fused ring.
"""
function fused_rings(mol::SimpleMolGraph)
    has_descriptor(mol, :fused_rings) && return get_descriptor(mol, :fused_rings)
    cobr = setdiff(Set(edges(mol)), bridges(mol.graph))
    subg, vmap = induced_subgraph(mol.graph, cobr)
    return  [vmap[c] for c in connected_components(subg)]
end


"""
    which_fused_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`fusedrings`](@ref) membership of
1 to ``n``th atoms of the given molecule.

Fused ring membership is represented as a set of fused ring indices assigned to each fused rings.
This means atoms that have the same fused ring index belong to the same fused ring.
"""
function which_fused_ring(mol::SimpleMolGraph)
    has_descriptor(mol, :v_which_fused_ring) && return get_descriptor(mol, :v_which_fused_ring)
    arr = init_node_descriptor(Vector{Int}, mol)
    for (i, conn) in enumerate(fused_rings(mol))
        for n in conn
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    smallest_ring(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the size of the smallest [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to. 

If the node is not in a ring, the value would be 0.
This property corresponds to SMARTS `r` query.
"""
function smallest_ring(mol::SimpleMolGraph)
    has_descriptor(mol, :v_smallest_ring) && return get_descriptor(mol, :v_smallest_ring)
    sssr_ = sssr(mol)
    whichring_ = which_ring(mol)
    arr = init_node_descriptor(Int, mol)
    for i in vertices(mol)
        rs = whichring_[i]
        isempty(rs) && continue
        arr[i] = minimum(length, sssr_[rs])
    end
    return arr
end


"""
    ring_count(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to.

This property corresponds to SMARTS `R` query.
"""
ring_count(mol::MolGraph) = get(descriptors(mol), :v_ring_count, length.(which_ring(mol)))
ring_count(mol::EditableMolGraph) = Dict(i => length(v) for (i, v) in which_ring(mol))


"""
    is_in_ring(mol::MolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th atoms of
the given molecule belong to a ring or not.
"""
is_in_ring(mol::MolGraph) = get(descriptors(mol), :v_is_in_ring, .!isempty.(which_ring(mol)))
is_in_ring(mol::EditableMolGraph) = Dict(i => !isempty(v) for (i, v) in which_ring(mol))


"""
    is_edge_in_ring(mol::MolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th bonds of
the given molecule belong to a ring or not.
"""
is_edge_in_ring(mol::MolGraph) = get(descriptors(mol), :e_is_in_ring, .!isempty.(edge_which_ring(mol)))
is_edge_in_ring(mol::EditableMolGraph) = Dict(i => !isempty(e) for (i, e) in edge_which_ring(mol))




# Primary properties

"""
    atom_symbol(mol::MolGraph) -> Vector{Symbol}
    atom_symbol(mol::EditableMolGraph) -> Vector{Symbol}

Return a vector of size ``n`` representing atom symbols of 1 to ``n``th atoms of
the given molecule.
"""
atom_symbol(mol::MolGraph) = get(descriptors(mol), :v_symbol, getproperty.(vprops(mol), :symbol))
atom_symbol(mol::EditableMolGraph) = Dict(i => get_prop(mol, i, :symbol) for i in vertices(mol))


"""
    charge(mol::MolGraph) -> Vector{Int}
    charge(mol::EditableMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom charges of 1 to ``n``th atoms of
the given molecule.
"""
charge(mol::MolGraph) = get(descriptors(mol), :v_charge, getproperty.(vprops(mol), :charge))
charge(mol::EditableMolGraph) = Dict(i => get_prop(mol, i, :charge) for i in vertices(mol))


"""
    multiplicity(mol::MolGraph) -> Vector{Int}
    multiplicity(mol::EditableMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom multiplicities of 1 to ``n``th atoms of
the given molecule (1: non-radical, 2: radical, 3: biradical).
"""
multiplicity(mol::MolGraph) = get(descriptors(mol), :v_multiplicity, getproperty.(vprops(mol), :multiplicity))
multiplicity(mol::EditableMolGraph) = Dict(i => get_prop(mol, i, :multiplicity) for i in vertices(mol))


"""
    bond_order(mol::MolGraph) -> Vector{Int}
    bond_order(mol::EditableMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing bond order of 1 to ``n``th bonds of
the given molecule.
"""
bond_order(mol::MolGraph) = get(descriptors(mol), :e_order, getproperty.(eprops(mol), :order))
bond_order(mol::EditableMolGraph) = Dict(i => get_prop(mol, i, :order) for i in edges(mol))


# mass -> src/mass.jl
# coords -> src/coords.jl


# Secondary properties (descriptors)

# Valence

"""
    lone_pair(mol::MolGraph) -> Vector{Union{Int,Nothing}}

Return a vector of size ``n`` representing the number of lone pairs of
1 to ``n``th atoms of the given molecule.

The number of lone pair in inorganic atoms would be `nothing`.
The result can take negative value if the atom has empty shells (e.g. B).
"""
function lone_pair(mol::SimpleMolGraph)
    has_descriptor(mol, :v_lone_pair) && return get_descriptor(mol, :v_lone_pair)
    arr = init_node_descriptor(Union{Int,Nothing}, mol)
    atomsymbol_ = atom_symbol(mol)
    charge_ = charge(mol)
    for i in vertices(mol)
        haskey(LONEPAIR_COUNT, atomsymbol_[i]) || continue
        arr[i] = LONEPAIR_COUNT[atomsymbol_[i]] - charge_[i]
    end
    return arr
end
lone_pair!(mol::MolGraph) = set_descriptor!(mol, :v_lone_pair, lone_pair(mol))


"""
    apparent_valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total bond order
incident to 1 to ``n``th atoms of the given molecule.
"""
function apparent_valence(mol::SimpleMolGraph)
    has_descriptor(mol, :v_apparent_valence) && return get_descriptor(mol, :v_apparent_valence)
    arr = init_node_descriptor(Int, mol)
    bondorder_ = bond_order(mol)
    for e in edges(mol)
        arr[src(e)] += bondorder_[edge_rank(mol, e)]
        arr[dst(e)] += bondorder_[edge_rank(mol, e)]
    end
    return arr
end
apparent_valence!(mol::MolGraph) = set_descriptor!(mol, :v_apparent_valence, apparent_valence(mol))


"""
    valence(mol::SimpleMolGraph) -> Vector{Union{Int,Nothing}}

Return a vector of size ``n`` representing the intrinsic valence of
1 to ``n``th atoms of the given molecule.

The number of implicit hydrogens would be calculated based on the valence.
The valence value in inorganic atoms would be `nothing`.
This property corresponds to SMARTS `v` query.
"""
function valence(mol::SimpleMolGraph)
    has_descriptor(mol, :v_valence) && return get_descriptor(mol, :v_valence)
    arr = init_node_descriptor(Union{Int,Nothing}, mol)
    atomsymbol_ = atom_symbol(mol)
    lonepair_ = lone_pair(mol)
    for i in vertices(mol)
        lonepair_[i] === nothing && continue
        if atomsymbol_[i] === :H
            arr[i] = 1
        else
            arr[i] = 4 - abs(lonepair_[i])
        end
    end
    return arr
end
valence!(mol::MolGraph) = set_descriptor!(mol, :v_valence, valence(mol))


"""
    explicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of explicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Explicit" means hydrogens are explicitly represented as graph nodes.
"""
function explicit_hydrogens(mol::SimpleMolGraph)
    arr = init_node_descriptor(Int, mol)
    atomsymbol_ = atom_symbol(mol)
    for i in vertices(mol)
        atomsymbol_[i] === :H || continue
        for nbr in neighbors(mol, i)
            arr[nbr] += 1
        end
    end
    return arr
end


"""
    implicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of implicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Implicit" means hydrogens are not represented as graph nodes,
but it can be infered from the intrinsic valence of typical organic atoms.
"""
function implicit_hydrogens(mol::SimpleMolGraph)
    arr = init_node_descriptor(Int, mol)
    valence_ = valence(mol)
    apparent_ = apparent_valence(mol)
    for i in vertices(mol)
        valence_[i] === nothing && continue
        arr[i] = max(0, valence_[i] - apparent_[i])
    end
    return arr
end


"""
    heavyatoms(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of non-hydrogen atoms
connected to 1 to ``n``th atoms of the given molecule.
"""
heavy_atoms(mol::SimpleMolGraph) = degree(mol) - explicit_hydrogens(mol)


"""
    hydrogenconnected(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total hydrogens
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `H` query.
"""
total_hydrogens(mol::SimpleMolGraph) = explicit_hydrogens(mol) + implicit_hydrogens(mol)


"""
    connectivity(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total atoms
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `X` query.
"""
connectivity(mol::SimpleMolGraph) = degree(mol) + implicit_hydrogens(mol)




# Hydrogen bond donor/acceptor

function is_hydrogen_acceptor(mol::SimpleMolGraph)
    ac = (sym, lp) -> lp === nothing ? false : sym in (:N, :O, :F) && lp > 0
    return ac.(atom_symbol(mol), lone_pair(mol))
end


"""
    hacceptorcount(mol::SimpleMolGraph) -> Int

Return the total number of hydrogen bond acceptors (N, O and F).
"""
hydrogen_acceptor_count(mol::SimpleMolGraph) = reduce(+, is_hydrogen_acceptor(mol); init=0)


function is_hydrogen_donor(mol::SimpleMolGraph)
    dc = (sym, h) -> sym in (:N, :O) && h > 0
    return dc.(atom_symbol(mol), total_hydrogens(mol))
end


"""
    hdonorcount(mol::SimpleMolGraph) -> Int

Return the total number of hydrogen bond donors (O and N attached to hydrogens).
"""
hydrogen_donor_count(mol::SimpleMolGraph) = reduce(+, is_hydrogen_donor(mol); init=0)




# Rotatable bonds

"""
    isrotatable(mol::SimpleMolGraph)

Return a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule are rotatable or not.
"""
function is_rotatable(mol::SimpleMolGraph)
    degree_ = degree(mol)
    edgeinring_ = is_edge_in_ring(mol)
    bondorder_ = bond_order(mol)
    arr = Vector{Bool}(undef, ne(mol))
    for (i, e) in enumerate(edges(mol))
        arr[i] = (!edgeinring_[i] && bondorder_[i] == 1
            && degree_[src(e)] != 1 && degree_[dst(e)] != 1)
    end
    return arr
end


"""
    rotatablecount(mol::SimpleMolGraph) -> Int

Return the total number of rotatable bonds.
"""
rotatable_count(mol::SimpleMolGraph) = reduce(+, is_rotatable(mol); init=0)




# Composition

"""
    atomcounter(mol::SimpleMolGraph) -> Dict{Symbol,Int}

Count the number of atoms and return symbol => count dict.
"""
function atom_counter(mol::SimpleMolGraph)
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
    heavyatomcount(mol::SimpleMolGraph) -> Int

Return the total number of non-hydrogen atoms.
"""
function heavy_atom_count(mol::SimpleMolGraph)
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
function molecular_formula(mol::SimpleMolGraph)
    counter = atom_counter(mol)
    return write_formula(counter)
end


"""
    empirical_formula(mol::MolGraph) -> String

Return the empirical formula in Hill system.
"""
function empirical_formula(mol::SimpleMolGraph)
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
    pi_electron(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the number of ``\\pi`` electrons
of 1 to ``n``th atoms of the given molecule.

The counting of ``\\pi`` electrons is based on the following rules.

- Any atom incident to a double bond -> +1
- Any atom incident to two double bond -> +2
- Any atom incident to a triple bond -> +2
- Any other uncharged N, O or S that are neighbor of multiple bonds -> +2

These rules are applied for only typical organic atoms. The values for inorganic atoms will be 0.
"""
function pi_electron(mol::SimpleMolGraph)
    atomsymbol_ = atom_symbol(mol)
    charge_ = charge(mol)
    pie_ = apparent_valence(mol) - degree(mol)
    arr = zeros(Int, nv(mol))
    for i in vertices(mol)
        arr[i] = pie_[i]
        atomsymbol_[i] in (:N, :O, :S) || continue
        for nbr in neighbors(mol, i)
            if pie_[i] == 0 && pie_[nbr] > 0 && charge_[i] == 0
                arr[i] = 2
                break
            end
        end
    end
    return arr
end


"""
    hybridization(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the orbital hybridization symbols
(`:sp3`, `:sp2`, `:sp` or `:none`) of 1 to ``n``th atoms of the given molecule.

The hybridization value in inorganic atoms and non-typical organic atoms will be `:none`
(e.g. s, sp3d and sp3d2 orbitals).
"""
function hybridization(mol::SimpleMolGraph)
    atomsymbol_ = atom_symbol(mol)
    pielectron_ = pi_electron(mol)
    connectivity_ = connectivity(mol)
    lonepair_ = lone_pair(mol)
    arr = Vector{Symbol}(undef, nv(mol))
    for i in vertices(mol)
        if lonepair_[i] === nothing || atomsymbol_[i] === :H
            arr[i] = :none
            continue
        end
        orbitals = connectivity_[i] + lonepair_[i]
        if orbitals == 4
            if (atomsymbol_[i] in [:N, :O] && pielectron_[i] == 2)
                arr[i] = :sp2  # adjacent to conjugated bonds
            else
                arr[i] = :sp3
            end
        elseif orbitals == 3
            arr[i] = :sp2
        elseif orbitals == 2
            arr[i] = :sp
        else
            arr[i] = :none
        end
    end
    return arr
end




# Aromaticity

function is_ring_aromatic(mol::SimpleMolGraph)
    has_descriptor(mol, :is_ring_aromatic) && return get_descriptor(mol, :is_ring_aromatic)
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
is_ring_aromatic!(mol::MolGraph) = set_descriptor!(mol, :is_ring_aromatic, is_ring_aromatic(mol))

"""
    is_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th atoms
of the given molecule belong to an aromatic ring or not.

Some kind of aromaticity resulting from long conjugated chains and charge
delocalization may be unrecognizable. Also, non-classical aromaticity
such as Moebius aromaticity is not considered.
"""
function is_aromatic(mol::SimpleMolGraph)
    arr = falses(nv(mol))
    for ring in sssr(mol)[findall(is_ring_aromatic(mol))]
        arr[ring] .= true
    end
    return arr
end


"""
    is_edge_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule belong to an aromatic ring or not.

See [`isaromatic`](@ref).
"""
function is_edge_aromatic(mol::SimpleMolGraph)
    arr = falses(ne(mol))
    for ring in sssr(mol)[findall(is_ring_aromatic(mol))]
        for i in 1:(length(ring) - 1)
            arr[edge_rank(mol, ring[i], ring[i + 1])] = true
        end
        arr[edge_rank(mol, ring[1], ring[end])] = true
    end
    return arr
end



"""
    precalculate!(mol::MolGraph)

Convenient method to pre-calculate and cache performance bottleneck descriptors.
"""
function precalculate!(mol::MolGraph)
    sssr!(mol)
    lone_pair!(mol)
    apparent_valence!(mol)
    valence!(mol)
    is_ring_aromatic!(mol)
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