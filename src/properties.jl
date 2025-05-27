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
    pi_electron, pi_delocalized, hybridization, hybridization_delocalized,
    is_ring_aromatic, is_ring_aromatic!, is_aromatic, is_edge_aromatic


const LONEPAIR_COUNT = Dict(
    :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
    :Si => 0, :P => 1, :S => 2, :Cl => 3,
    :As => 1, :Se => 2, :Br => 3, :I => 3
)

# Atoms that can be hydrogen acceptors/donors
const HYDROGEN_ACCEPTOR_ATOMS = (:N, :O, :F)
const HYDROGEN_DONOR_ATOMS = (:N, :O)

# for pi electron count, hybridization and aromaticity calculation
# TODO: :P, :Se, :Te? 
const SP2_CONJUGATING_HETEROATOMS = (:O, :N, :S)


# Molecular graph topology descriptors

"""
    degree(mol::SimpleMolGraph{T}) -> Vector{T}

Return a vector of size ``n`` representing the node degree of the molecular graph
of 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `D` query.
"""
function Graphs.degree(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_degree) && return get_cache(mol, :v_degree)
    return degree(mol.graph)
end

degree!(mol::SimpleMolGraph) = set_cache!(mol, :v_degree, degree(mol.graph))


"""
    sssr(mol::SimpleMolGraph{T}) -> Vector{Vector{T}}

Return vectors of ring nodes representing small set of smallest rings (SSSR).

See [`mincyclebasis`](@ref).
"""
function sssr(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :sssr) && return get_cache(mol, :sssr)
    return mincyclebasis(mol.graph)
end

sssr!(mol::SimpleMolGraph) = set_cache!(mol, :sssr, mincyclebasis(mol.graph))


"""
    which_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th nodes of the given graph.

SSSR membership is represented as a vector of SSSR indices assigned to each rings.
This means nodes that have the same SSSR index belong to the same SSSR.
"""
function which_ring(mol::SimpleMolGraph{T,V,E}) where {T,V,E}
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_which_ring) && return get_cache(mol, :v_which_ring)
    arr = [Int[] for _ in vertices(mol)]
    for (i, cyc) in enumerate(sssr(mol))
        for n in cyc
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    edge_which_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th bonds of the given molecule.

SSSR membership is represented as a set of SSSR indices assigned to each rings.
This means bonds that have the same SSSR index belong to the same SSSR.
"""
function edge_which_ring(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_which_ring) && return get_cache(mol, :e_which_ring)
    arr = [Int[] for _ in 1:ne(mol)]
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
function fused_rings(g::SimpleGraph)
    cobr = setdiff(Set(edges(g)), bridges(g))
    subg, vmap = induced_subgraph(g, collect(cobr))
    return  [vmap[c] for c in connected_components(subg)]
end

function fused_rings(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :fused_rings) && return get_cache(mol, :fused_rings)
    return fused_rings(mol.graph)
end


"""
    which_fused_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`fused_rings`](@ref) membership of
1 to ``n``th atoms of the given molecule.

Fused ring membership is represented as a set of fused ring indices assigned to each fused rings.
This means atoms that have the same fused ring index belong to the same fused ring.
"""
function which_fused_ring(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_which_fused_ring) && return get_cache(mol, :v_which_fused_ring)
    arr = [Int[] for _ in vertices(mol)]
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
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_smallest_ring) && return get_cache(mol, :v_smallest_ring)
    sssr_ = sssr(mol)
    whichring_ = which_ring(mol)
    arr = zeros(Int, nv(mol))
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
function ring_count(mol::MolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_ring_count) && return get_cache(mol, :v_ring_count)
    return length.(which_ring(mol))
end


"""
    is_in_ring(mol::MolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th atoms of
the given molecule belong to a ring or not.
"""
function is_in_ring(mol::MolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_is_in_ring) && return get_cache(mol, :v_is_in_ring)
    return .!isempty.(which_ring(mol))
end


"""
    is_edge_in_ring(mol::MolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th bonds of
the given molecule belong to a ring or not.
"""
function is_edge_in_ring(mol::MolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_is_in_ring) && return get_cache(mol, :e_is_in_ring)
    return .!isempty.(edge_which_ring(mol))
end



"""
    bmscaffold(mol::MolGraph) -> Vector{Int}

Return Bemis-Murcko scaffold as a vector of scaffold vertices.
"""
function bmscaffold(g::SimpleGraph)
    (length(connected_components(g)) == 1
        || error("bmscaffold for a disconnected graph cannot be determined"))
    degree_ = degree(g)
    while true
        terms = findall(degree_ .== 1)
        isempty(terms) && break
        for t in terms
            degree_[t] = 0
            for nbr in neighbors(g, t)
                degree_[nbr] == 0 && continue
                degree_[nbr] -= 1
            end
        end
    end
    return findall(degree_ .!= 0)
end
bmscaffold(mol::SimpleMolGraph) = bmscaffold(mol.graph)


"""
    scaffold_fragments(g::SimpleGraph) -> Vector{Int}

Return scaffold fragments. Intended for a MMP derivative.
"""
function scaffold_fragments(g::SimpleGraph)
    bms = bmscaffold(g)
    isempty(bms) && return scaffolds  # no scaffold
    queue = [bms]  # BFS
    scaffolds = Vector{Int}[bms]
    dag = SimpleDiGraph{Int}(1)  # fragmentation relationship (parent -> child)
    while !isempty(queue)
        # generate new fragments
        scaffold = pop!(queue)
        scf, scfvmap = induced_subgraph(g, scaffold)
        parent =findfirst(x -> x == scaffold, scaffolds)
        for fr in fused_rings(scf)
            # enumerate connections to the fused ring
            nbrs = setdiff(union([neighbors(scf, v) for v in fr]...), fr)
            length(nbrs) == 1 || continue
            # remove a fused ring
            subg, vmap = induced_subgraph(g, setdiff(scaffold, scfvmap[fr]))
            fragment = vmap[bmscaffold(subg)]
            if fragment in scaffolds  # duplicate, skip the branch
                dst = findfirst(x -> x == fragment, scaffolds)
                add_edge!(dag, parent, dst)
                continue
            end
            pushfirst!(queue, fragment)
            push!(scaffolds, fragment)
            add_vertex!(dag)
            add_edge!(dag, parent, length(scaffolds))
            if length(scaffolds) == 400
                @info "max fragment count reached"
                break
            end
        end
    end
    return scaffolds, dag
end



# Primary properties

"""
    atom_symbol(mol::MolGraph) -> Vector{Symbol}

Return a vector of size ``n`` representing atom symbols of 1 to ``n``th atoms of
the given molecule.
"""
function atom_symbol(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_symbol) && return get_cache(mol, :v_symbol)
    return [atom_symbol(props(mol, i)) for i in vertices(mol)]
end

atom_symbol!(mol::SimpleMolGraph
    ) = set_cache!(mol, :v_symbol, [atom_symbol(props(mol, i)) for i in vertices(mol)])


"""
    atom_number(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom numbers of 1 to ``n``th atoms of
the given molecule.
"""
function atom_number(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_number) && return get_cache(mol, :v_number)
    return [atom_number(props(mol, i)) for i in vertices(mol)]
end

atom_number!(mol::SimpleMolGraph
    ) = set_cache!(mol, :v_number, [atom_number(props(mol, i)) for i in vertices(mol)])


"""
    charge(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom charges of 1 to ``n``th atoms of
the given molecule.
"""
function charge(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_charge) && return get_cache(mol, :v_charge)
    return [charge(props(mol, i)) for i in vertices(mol)]
end

# kekulize! or charge standardization (e.g. polarize!) would be reset
charge!(mol::SimpleMolGraph
    ) = set_cache!(mol, :v_charge, [charge(props(mol, i)) for i in vertices(mol)])


"""
    multiplicity(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom multiplicities of 1 to ``n``th atoms of
the given molecule (1: non-radical, 2: radical, 3: biradical).
"""
function multiplicity(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_multiplicity) && return get_cache(mol, :v_multiplicity)
    return [multiplicity(props(mol, i)) for i in vertices(mol)]
end


"""
    bond_order(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing bond order of 1 to ``n``th bonds of
the given molecule.
"""
function bond_order(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_order) && return get_cache(mol, :e_order)
    return [bond_order(props(mol, e)) for e in edges(mol)]
end

# mass -> src/mass.jl
# coords -> src/coords.jl


# Secondary properties (descriptors)

# Valence

"""
    apparent_valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total bond order
incident to 1 to ``n``th atoms of the given molecule.
"""
function apparent_valence(g, order_arr)
    arr = fill(zero(Int), nv(g))
    er = Dict(e => i for (i, e) in enumerate(edges(g)))
    for e in edges(g)
        arr[src(e)] += order_arr[er[e]]
        arr[dst(e)] += order_arr[er[e]]
    end
    return arr
end

function apparent_valence(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_apparent_valence) && return get_cache(mol, :v_apparent_valence)
    return apparent_valence(mol.graph, bond_order(mol))
end

apparent_valence!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_apparent_valence, apparent_valence(mol.graph, bond_order(mol)))


"""
    valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the intrinsic valence of
1 to ``n``th atoms of the given molecule.

The number of implicit hydrogens would be calculated based on the valence.
The valence of a hypervalent atom or a non-organic atom is the same as
its `apparent_valence`. This property corresponds to SMARTS `v` query.
"""
function valence(symbol_arr, charge_arr, apparent_valence_arr)
    arr = fill(zero(Int), length(symbol_arr))
    for i in 1:length(symbol_arr)
        if haskey(LONEPAIR_COUNT, symbol_arr[i])
            val = 4 - abs(LONEPAIR_COUNT[symbol_arr[i]] - charge_arr[i])
            arr[i] = max(val, apparent_valence_arr[i])
        else
            arr[i] = apparent_valence_arr[i]
        end
    end
    return arr
end

function valence(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_valence) && return get_cache(mol, :v_valence)
    return valence(atom_symbol(mol), charge(mol), apparent_valence(mol))
end

valence!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_valence,
    valence(atom_symbol(mol), charge(mol), apparent_valence(mol))
)


"""
    explicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of explicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Explicit" means hydrogens are explicitly represented as graph nodes.
"""
function explicit_hydrogens(g, symbol_arr)
    arr = fill(zero(Int), nv(g))
    for i in vertices(g)
        symbol_arr[i] === :H || continue
        for nbr in neighbors(g, i)
            arr[nbr] += 1
        end
    end
    return arr
end

function explicit_hydrogens(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_explicit_hydrogens) && return get_cache(mol, :v_explicit_hydrogens)
    return explicit_hydrogens(mol.graph, atom_symbol(mol))
end

explicit_hydrogens!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_explicit_hydrogens, explicit_hydrogens(mol.graph, atom_symbol(mol)))


"""
    implicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of implicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Implicit" means hydrogens are not represented as graph nodes,
but it can be infered from the intrinsic valence of typical organic atoms.
"""
function implicit_hydrogens(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_implicit_hydrogens) && return get_cache(mol, :v_implicit_hydrogens)
    return valence(mol) - apparent_valence(mol)
end

implicit_hydrogens!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_implicit_hydrogens, valence(mol) - apparent_valence(mol))


"""
    heavy_atoms(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of non-hydrogen atoms
connected to 1 to ``n``th atoms of the given molecule.
"""
function heavy_atoms(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_heavy_atoms) && return get_cache(mol, :v_heavy_atoms)
    return degree(mol) - explicit_hydrogens(mol)
end

heavy_atoms!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_heavy_atoms, degree(mol) - explicit_hydrogens(mol))



"""
    total_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total hydrogens
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `H` query.
"""
function total_hydrogens(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_total_hydrogens) && return get_cache(mol, :v_total_hydrogens)
    return explicit_hydrogens(mol) + implicit_hydrogens(mol)
end

total_hydrogens!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_total_hydrogens, explicit_hydrogens(mol) + implicit_hydrogens(mol))


"""
    connectivity(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total atoms
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `X` query.
"""
function connectivity(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_connectivity) && return get_cache(mol, :v_connectivity)
    return degree(mol) + implicit_hydrogens(mol)
end

connectivity!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_connectivity, degree(mol) + implicit_hydrogens(mol))


"""
    lone_pair(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of lone pairs of
1 to ``n``th atoms of the given molecule.
"""
function lone_pair(symbol_arr, charge_arr, connectivity_arr)
    arr = fill(zero(Int), length(symbol_arr))
    for i in 1:length(symbol_arr)
        if haskey(LONEPAIR_COUNT, symbol_arr[i])
            excess = 4 - connectivity_arr[i]
            intrinsic = LONEPAIR_COUNT[symbol_arr[i]] - charge_arr[i]
            arr[i] = max(min(excess, intrinsic), 0)
        else
            arr[i] = 0
        end
    end
    return arr
end

function lone_pair(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_lone_pair) && return get_cache(mol, :v_lone_pair)
    return lone_pair(atom_symbol(mol), charge(mol), connectivity(mol))
end

lone_pair!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_lone_pair, lone_pair(atom_symbol(mol), charge(mol), connectivity(mol)))




# Hydrogen bond donor/acceptor

function is_hydrogen_acceptor(symbol_arr, lone_pair_arr)
    arr = falses(length(symbol_arr))
    for i in 1:length(symbol_arr)
        if symbol_arr[i] in HYDROGEN_ACCEPTOR_ATOMS && lone_pair_arr[i] > 0
            arr[i] = true
        end
    end
    return arr
end

function is_hydrogen_acceptor(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_is_hydrogen_acceptor) && return get_cache(mol, :v_is_hydrogen_acceptor)
    return is_hydrogen_acceptor(atom_symbol(mol), lone_pair(mol))
end

is_hydrogen_acceptor!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_is_hydrogen_acceptor, is_hydrogen_acceptor(atom_symbol(mol), lone_pair(mol)))


"""
    hydrogen_acceptor_count(mol::SimpleMolGraph) -> Int

Return the total number of hydrogen bond acceptors (N, O and F).
"""
hydrogen_acceptor_count(mol::SimpleMolGraph) = reduce(+, is_hydrogen_acceptor(mol); init=0)


function is_hydrogen_donor(symbol_arr, total_hydrogens_arr)
    arr = falses(length(symbol_arr))
    for i in 1:length(symbol_arr)
        if symbol_arr[i] in HYDROGEN_DONOR_ATOMS && total_hydrogens_arr[i] > 0
            arr[i] = true
        end
    end
    return arr
end

function is_hydrogen_donor(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_is_hydrogen_donor) && return get_cache(mol, :v_is_hydrogen_donor)
    return is_hydrogen_donor(atom_symbol(mol), total_hydrogens(mol))
end

is_hydrogen_donor!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_is_hydrogen_donor, is_hydrogen_donor(atom_symbol(mol), total_hydrogens(mol)))


"""
    hydrogen_donor_count(mol::SimpleMolGraph) -> Int

Return the total number of hydrogen bond donors (O and N attached to hydrogens).
"""
hydrogen_donor_count(mol::SimpleMolGraph) = reduce(+, is_hydrogen_donor(mol); init=0)




# Rotatable bonds

"""
    is_rotatable(mol::SimpleMolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule are rotatable or not.
"""
function is_rotatable(edge_list, degree_arr, edge_in_ring_arr, order_arr)
    arr = falses(length(edge_list))
    for (i, e) in enumerate(edge_list)
        if (!edge_in_ring_arr[i] && order_arr[i] == 1
                && degree_arr[src(e)] != 1 && degree_arr[dst(e)] != 1)
            arr[i] = true
        end
    end
    return arr
end

function is_rotatable(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_is_rotatable) && return get_cache(mol, :e_is_rotatable)
    return is_rotatable(edges(mol), degree(mol), is_edge_in_ring(mol), bond_order(mol))
end

is_rotatable!(mol::SimpleMolGraph) = set_cache!(
    mol, :e_is_rotatable,
    is_rotatable(edges(mol), degree(mol), is_edge_in_ring(mol), bond_order(mol)))


"""
    rotatable_count(mol::SimpleMolGraph) -> Int

Return the total number of rotatable bonds.
"""
rotatable_count(mol::SimpleMolGraph) = reduce(+, is_rotatable(mol); init=0)




# Composition

"""
    atom_counter(mol::SimpleMolGraph) -> Dict{Symbol,Int}

Count the number of atoms and return symbol => count dict.
"""
function atom_counter(mol::SimpleMolGraph)
    counter = Dict{Symbol,Int}()
    for sym in atom_symbol(mol)
        if !haskey(counter, sym)
            counter[sym] = 0
        end
        counter[sym] += 1
    end
    hcnt = reduce(+, implicit_hydrogens(mol); init=0)
    if hcnt > 0
        if !haskey(counter, :H)
            counter[:H] = 0
        end
        counter[:H] += hcnt
    end
    return counter
end


"""
    heavy_atom_count(mol::SimpleMolGraph) -> Int

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




# Hybridization and aromaticity

"""
    hybridization(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the orbital hybridization symbols
(`:sp3`, `:sp2`, `:sp` or `:none`) of 1 to ``n``th atoms of the given molecule.

The hybridization value in inorganic atoms and non-typical organic atoms will be `:none`
(e.g. s, sp3d and sp3d2 orbitals). Note that this is a simplified geometry descriptor
for substructure matching and does not reflect actual molecular orbital hybridization.
"""
function hybridization(g, symbol_arr, valence_arr, connectivity_arr, lone_pair_arr)
    arr = fill(:none, length(lone_pair_arr))
    hybmap = Dict(4 => :sp3, 3 => :sp2, 2 => :sp)
    for i in 1:length(arr)
        cnt = connectivity_arr[i] + lone_pair_arr[i]
        arr[i] = get(hybmap, cnt, :none)
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

function hybridization(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_hybridization) && return get_cache(mol, :v_hybridization)
    return hybridization(
        mol.graph, atom_symbol(mol), valence(mol), connectivity(mol), lone_pair(mol))
end

hybridization!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_hybridization,
    hybridization(mol.graph, atom_symbol(mol), valence(mol), connectivity(mol), lone_pair(mol)))

hybridization_delocalized = hybridization  # TODO: for backward compatibility. to be removed.


"""
    pi_electron(mol::SimpleMolGraph) -> Vector{Int}

Returns a vector of size ``n`` representing the number of ``\\pi`` electrons
of 1 to ``n``th atoms of the given molecule.

The number of ``\\pi`` electrons is calculated as `valence` - `connectivity`.
Typically, each atom connected to double bonds adds one pi electron for each, and each
atom connected to a triple bond adds two pi electrons.
"""
function pi_electron(valence_arr, connectivity_arr, lone_pair_arr, hyb_arr)
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

function pi_electron(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_pi_electron) && return get_cache(mol, :v_pi_electron)
    return pi_electron(valence(mol), connectivity(mol), lone_pair(mol), hybridization(mol))
end

pi_electron!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_pi_electron,
    pi_electron(valence(mol), connectivity(mol)), lone_pair(mol), hybridization(mol))

pi_delocalized = pi_electron  # TODO: for backward compatibility. to be removed.


"""
    is_ring_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether first to ``n``-th rings
of a given molecule are aromatic or not.

This is a binary descriptor based on a chemoinformatic algorithm and may not reflect
actual molecular orbitals. Atypical aromaticities such as Moebius aromaticity are not considered.
"""
function is_ring_aromatic(g, sssr_, which_ring_arr, symbol_arr, order_arr, hyb_arr, pi_arr)
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

function is_ring_aromatic(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :is_ring_aromatic) && return get_cache(mol, :is_ring_aromatic)
    return is_ring_aromatic(
        mol.graph, sssr(mol), edge_which_ring(mol), atom_symbol(mol), bond_order(mol),
        hybridization(mol), pi_electron(mol))
end

function is_ring_aromatic!(mol::SimpleMolGraph)
    set_cache!(mol, :is_ring_aromatic, is_ring_aromatic(
        mol.graph, sssr(mol), edge_which_ring(mol), atom_symbol(mol), bond_order(mol),
        hybridization(mol), pi_electron(mol)))
end


"""
    is_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th atoms
of the given molecule belong to an aromatic ring or not.

See [`is_ring_aromatic`](@ref).
"""
function is_aromatic(g, sssr_, is_ring_arom)
    arr = falses(nv(g))
    for ring in sssr_[findall(is_ring_arom)]
        arr[ring] .= true
    end
    return arr
end


function is_aromatic(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_is_aromatic) && return get_cache(mol, :v_is_aromatic)
    return is_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol))
end

is_aromatic!(mol::SimpleMolGraph) = set_cache!(
    mol, :v_is_aromatic, is_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol)))


"""
    is_edge_aromatic(mol::SimpleMolGraph) -> Vector{Bool}

Returns a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule belong to an aromatic ring or not.

See [`is_ring_aromatic`](@ref).
"""
function is_edge_aromatic(g, sssr_, is_ring_arom)
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

function is_edge_aromatic(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :e_is_aromatic) && return get_cache(mol, :e_is_aromatic)
    return is_edge_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol))
end

is_edge_aromatic!(mol::SimpleMolGraph) = set_cache!(
    mol, :e_is_aromatic, is_edge_aromatic(mol.graph, sssr(mol), is_ring_aromatic(mol)))
