#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

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


# Primary properties

"""
    atom_symbol(mol::MolGraph) -> Vector{Symbol}

Return a vector of size ``n`` representing atom symbols of 1 to ``n``th atoms of
the given molecule.
"""
function atom_symbol(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return [atom_symbol(props(mol, i)) for i in vertices(mol)]
end


"""
    atom_number(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom numbers of 1 to ``n``th atoms of
the given molecule.
"""
function atom_number(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return [atom_number(props(mol, i)) for i in vertices(mol)]
end


"""
    atom_charge(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom charges of 1 to ``n``th atoms of
the given molecule.
"""
function atom_charge(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return mol.gprops.descriptors.atom_charge
end

default_atom_charge!(mol::SimpleMolGraph) = setproperty!(
    mol.gprops.descriptors, :atom_charge,
    [atom_charge(props(mol, i)) for i in vertices(mol)]
)


"""
    multiplicity(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom multiplicities of 1 to ``n``th atoms of
the given molecule (1: non-radical, 2: radical, 3: biradical).
"""
function multiplicity(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return [multiplicity(props(mol, i)) for i in vertices(mol)]
end


"""
    bond_order(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing bond order of 1 to ``n``th bonds of
the given molecule.
"""
function bond_order(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return mol.gprops.descriptors.bond_order
end

default_bond_order!(mol::SimpleMolGraph) = setproperty!(
    mol.gprops.descriptors, :bond_order,
    [bond_order(props(mol, e)) for e in edges(mol)]
)

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
    dispatch_update!(mol)
    return mol.gprops.descriptors.apparent_valence
end

apparent_valence!(mol::SimpleMolGraph) = setproperty!(
    mol.gprops.descriptors, :apparent_valence,
    apparent_valence(mol.graph, bond_order(mol))
)


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
    dispatch_update!(mol)
    return mol.gprops.descriptors.valence
end

valence!(mol::SimpleMolGraph) = setproperty!(
    mol.gprops.descriptors, :valence,
    valence(atom_symbol(mol), atom_charge(mol), apparent_valence(mol))
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
    dispatch_update!(mol)
    return explicit_hydrogens(mol.graph, atom_symbol(mol))
end


"""
    implicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of implicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Implicit" means hydrogens are not represented as graph nodes,
but it can be infered from the intrinsic valence of typical organic atoms.
"""
function implicit_hydrogens(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return valence(mol) - apparent_valence(mol)
end


"""
    heavy_atoms(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of non-hydrogen atoms
connected to 1 to ``n``th atoms of the given molecule.
"""
function heavy_atoms(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return degree(mol) - explicit_hydrogens(mol)
end


"""
    total_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total hydrogens
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `H` query.
"""
function total_hydrogens(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return explicit_hydrogens(mol) + implicit_hydrogens(mol)
end


"""
    connectivity(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total atoms
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `X` query.
"""
function connectivity(mol::SimpleMolGraph)
    dispatch_update!(mol)
    return degree(mol) + implicit_hydrogens(mol)
end


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
    dispatch_update!(mol)
    return mol.gprops.descriptors.lone_pair
end

lone_pair!(mol::SimpleMolGraph) = setproperty!(
    mol.gprops.descriptors, :lone_pair,
    lone_pair(atom_symbol(mol), atom_charge(mol), connectivity(mol))
)




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
    dispatch_update!(mol)
    return is_hydrogen_acceptor(atom_symbol(mol), lone_pair(mol))
end


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
    dispatch_update!(mol)
    return is_hydrogen_donor(atom_symbol(mol), total_hydrogens(mol))
end


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
    dispatch_update!(mol)
    return is_rotatable(edges(mol), degree(mol), is_edge_in_ring(mol), bond_order(mol))
end


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
