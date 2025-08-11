#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


const SPn_ATOMS = (:B, :C, :N, :O, :Si, :P, :S)
# TODO: :P, :Se, :Te?
const SP2_CONJUGATING_HETEROATOMS = (:O, :N, :S)


"""
    hybridization(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the orbital hybridization symbols
(`:sp3`, `:sp2`, `:sp` or `:none`) of 1 to ``n``th atoms of the given molecule.

The hybridization value in inorganic atoms and non-typical organic atoms will be `:none`
(e.g. s, sp3d and sp3d2 orbitals). Note that this is a simplified geometry descriptor
for substructure matching and does not reflect actual molecular orbital hybridization.
"""
hybridization(mol::SimpleMolGraph) = hybridization(
    mol.graph, atom_symbol(mol), valence(mol), connectivity(mol), lone_pair(mol))

function hybridization(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return hybridization(
        mol.graph, atom_symbol(mol), valence(mol), connectivity(mol), lone_pair(mol))
end

function hybridization(
        g::SimpleGraph, symbol_arr::Vector{Symbol}, valence_arr::Vector{Int},
        connectivity_arr::Vector{Int}, lone_pair_arr::Vector{Int})
    arr = fill(:none, length(lone_pair_arr))
    hybmap = Dict(4 => :sp3, 3 => :sp2, 2 => :sp)
    for i in 1:length(arr)
        symbol_arr[i] in SPn_ATOMS || continue
        cnt = connectivity_arr[i] + lone_pair_arr[i]
        haskey(hybmap, cnt)  || continue
        arr[i] = hybmap[cnt]
    end
    # Hybridization of heteroatoms next to conjugated system
    for i in 1:length(arr)
        symbol_arr[i] in SP2_CONJUGATING_HETEROATOMS || continue
        (arr[i] == :sp3 && valence_arr[i] < 4 && lone_pair_arr[i] > 0) || continue
        for nbr in neighbors(g, i)
            arr[nbr] in (:sp, :sp2) || continue
            arr[i] = :sp2
            break
        end
    end
    return arr
end


"""
    pi_electron(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the number of ``\\pi`` electrons
of 1 to ``n``th atoms of the given molecule.

The number of ``\\pi`` electrons is calculated as `valence` - `connectivity`.
Typically, each atom connected to double bonds adds one pi electron for each, and each
atom connected to a triple bond adds two pi electrons.
"""
pi_electron(mol::SimpleMolGraph) = pi_electron(
    valence(mol), connectivity(mol), lone_pair(mol), hybridization(mol))

function pi_electron(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return pi_electron(valence(mol), connectivity(mol), lone_pair(mol), hybridization(mol))
end

function pi_electron(
        valence_arr::Vector{Int}, connectivity_arr::Vector{Int},
        lone_pair_arr::Vector{Int}, hyb_arr::Vector{Symbol})
    arr = fill(zero(Int), length(valence_arr))
    for i in 1:length(arr)
        pie = valence_arr[i] - connectivity_arr[i]
        if pie == 0 && lone_pair_arr[i] > 0 && hyb_arr[i] === :sp2
            arr[i] = 2
        else
            arr[i] = pie
        end
    end
    return arr
end


"""
    is_ring_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether first to ``n``-th rings
of a given molecule are aromatic or not.

This is a binary descriptor based on a chemoinformatic algorithm and may not reflect
actual molecular orbitals. Atypical aromaticities such as Moebius aromaticity are not considered.
"""
is_ring_aromatic(mol::SimpleMolGraph) = is_ring_aromatic(
    mol.graph, sssr(mol), edge_which_ring(mol),
    atom_symbol(mol), bond_order(mol),
    hybridization(mol), pi_electron(mol)
)

function is_ring_aromatic(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    if has_descriptor(mol, :is_ring_aromatic)
        return get_descriptor(mol, :is_ring_aromatic)
    end
    return is_ring_aromatic(
        mol.graph, sssr(mol), edge_which_ring(mol),
        atom_symbol(mol), bond_order(mol),
        hybridization(mol), pi_electron(mol)
    )
end

is_ring_aromatic!(mol::SimpleMolGraph) = setproperty!(
    mol.gprops.descriptors, :is_ring_aromatic,
    is_ring_aromatic(
        mol.graph, sssr(mol), edge_which_ring(mol),
        atom_symbol(mol), bond_order(mol),
        hybridization(mol), pi_electron(mol)
    )
)

function is_ring_aromatic(
        g::SimpleGraph, sssr_::Vector{Vector{Int}}, which_ring_arr::Vector{Vector{Int}},
        symbol_arr::Vector{Symbol}, order_arr::Vector{Int}, hyb_arr::Vector{Symbol},
        pi_arr::Vector{Int})
    # 1. evaluate each rings
    confirmed_ring = Int[]  # marked as aromatic
    not_aromatic = Int[]  # Unlikely to be aromatic (sp2 conjugation break)
    vs_confirmed = falses(nv(g))  # vertices belong to `confirmed_ring`
    vs_declined = falses(nv(g))  # vertices belong to `not_aromatic`
    rings_suspended = Vector{Int}[]  # depends on adjacent rings
    huckel_arr = copy(pi_arr)  # Huckel rule electron count
    er = Dict(e => i for (i, e) in enumerate(edges(g)))  # edge rank
    for (i, ring) in enumerate(sssr_)
        ring_sus = Int[]
        suspended = false
        broken = false
        for (i, r) in enumerate(ring)
            if hyb_arr[r] !== :sp2
                broken = true  # sp2 conjugation break
                break
            end
            # Check if double bonds are along with the ring or not
            outnbrs = setdiff(neighbors(g, r), ring)
            length(outnbrs) == 1 || continue
            outnbr = only(outnbrs)
            e = er[u_edge(g, r, outnbr)]
            order_arr[e] == 2 || continue
            # Process outgoing double bonds
            if length(which_ring_arr[e]) != 0
                suspended = true  # depends on adjacent rings
                push!(ring_sus, r)
            elseif symbol_arr[outnbr] === :O  # carbonyl
                huckel_arr[r] = 0
            else
                broken = true  # can not be aromatic
                break
            end
        end
        push!(rings_suspended, ring_sus)
        if broken
            push!(not_aromatic, i)
            vs_declined[sssr_[i]] .= true
        elseif !suspended && sum(huckel_arr[ring]) % 4 == 2
            push!(confirmed_ring, i)
            vs_confirmed[sssr_[i]] .= true
        end
    end
    # 2. Expand adjacent possible rings from confirmed ring
    found = true
    while found
        found = false
        for r in setdiff(1:length(sssr_), confirmed_ring, not_aromatic)
            sus = setdiff(rings_suspended[r], findall(vs_confirmed))
            isempty(sus) || continue
            rest = setdiff(sssr_[r], findall(vs_confirmed))
            cfcnt = length(sssr_[r]) - length(rest)
            if (sum(huckel_arr[rest]) + cfcnt) % 4 == 2
                push!(confirmed_ring, r)
                vs_confirmed[sssr_[r]] .= true
                found = true
            end
        end
    end
    # 3. Find indivisible fused aromatic rings
    ring_conn = SimpleGraph{Int}(length(sssr_))
    visited = vcat(confirmed_ring, not_aromatic)
    while length(visited) != length(sssr_)
        root = setdiff(1:length(sssr_), visited)[1]
        stack = [root]
        push!(visited, root)
        while !isempty(stack)
            n = popfirst!(stack)  # BFS
            for e in induced_subgraph_edges(g, sssr_[n])
                order_arr[er[e]] == 1 || continue
                rings = which_ring_arr[er[e]]
                length(rings) == 2 || continue
                nr = only(setdiff(rings, n))
                nr in visited && continue
                if pi_arr[e.src] == 2 || pi_arr[e.dst] == 2
                    # e.g. coelenterazine
                    add_edge!(ring_conn, n, nr)
                    push!(stack, nr)
                    push!(visited, nr)
                    continue
                end
                uv, vv = edge_neighbors(g, e)
                ud = []
                for u in uv
                    ur = er[u_edge(g, e.src, u)]
                    if order_arr[ur] == 2
                        push!(ud, ur)
                    end
                end
                vd = []
                for v in vv
                    vr = er[u_edge(g, e.dst, v)]
                    if order_arr[vr] == 2
                        push!(vd, vr)
                    end
                end
                (length(ud) != 1 || length(vd) != 1) && continue
                if which_ring_arr[ud[1]] != which_ring_arr[vd[1]]
                    add_edge!(ring_conn, n, nr)
                    push!(stack, nr)
                    push!(visited, nr)
                end
            end
        end
    end
    for conn in connected_components(ring_conn)
        length(conn) == 1 && continue
        ns = union([sssr_[n] for n in conn]...)
        any(vs_declined[ns]) && continue
        if sum(huckel_arr[ns]) % 4 == 2
            push!(confirmed_ring, conn...)
        end
    end
    res = falses(length(sssr_))
    for c in confirmed_ring
        res[c] = true
    end
    return res
end


"""
    is_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th atoms
of the given molecule belong to an aromatic ring or not.

See [`is_ring_aromatic`](@ref).
"""
is_aromatic(mol::SimpleMolGraph) = is_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol))

function is_aromatic(
        g::SimpleGraph, sssr_::Vector{Vector{Int}}, is_ring_arom::Vector{Bool})
    arr = falses(nv(g))
    for ring in sssr_[findall(is_ring_arom)]
        arr[ring] .= true
    end
    return arr
end



"""
    is_edge_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule belong to an aromatic ring or not.

See [`is_ring_aromatic`](@ref).
"""
is_edge_aromatic(mol::SimpleMolGraph) = is_edge_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol))

function is_edge_aromatic(
        g::SimpleGraph, sssr_::Vector{Vector{Int}}, is_ring_arom::Vector{Bool})
    arr = falses(ne(g))
    er = Dict(e => i for (i, e) in enumerate(edges(g)))
    for ring in sssr_[findall(is_ring_arom)]
        for i in 1:(length(ring) - 1)
            arr[er[u_edge(g, ring[i], ring[i + 1])]] = true
        end
        arr[er[u_edge(g, ring[1], ring[end])]] = true
    end
    return arr
end
