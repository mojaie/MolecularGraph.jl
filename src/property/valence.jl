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


# Primary properties

"""
    atom_symbol(mol::MolGraph) -> Vector{Symbol}

Return a vector of size ``n`` representing atom symbols of 1 to ``n``th atoms of
the given molecule.
"""
atom_symbol(mol::SimpleMolGraph) = Symbol[atom_symbol(props(mol, i)) for i in vertices(mol)]

function atom_symbol(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return Symbol[atom_symbol(props(mol, i)) for i in vertices(mol)]
end


"""
    atom_number(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom numbers of 1 to ``n``th atoms of
the given molecule.
"""
atom_number(mol::SimpleMolGraph) = Int[atom_number(props(mol, i)) for i in vertices(mol)]

function atom_number(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return Int[atom_number(props(mol, i)) for i in vertices(mol)]
end


"""
    atom_charge(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom charges of 1 to ``n``th atoms of
the given molecule.
"""
atom_charge(mol::SimpleMolGraph) = Int[atom_charge(props(mol, i)) for i in vertices(mol)]

function atom_charge(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    if has_descriptor(mol, :atom_charge)
        return get_descriptor(mol, :atom_charge)
    end
    return Int[atom_charge(props(mol, i)) for i in vertices(mol)]
end

default_atom_charge!(mol::ReactiveMolGraph) = setproperty!(
    mol.gprops.descriptors, :atom_charge,
    Int[atom_charge(props(mol, i)) for i in vertices(mol)]
)


"""
    multiplicity(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom multiplicities of 1 to ``n``th atoms of
the given molecule (1: non-radical, 2: radical, 3: biradical).
"""
multiplicity(mol::SimpleMolGraph) = Int[multiplicity(props(mol, i)) for i in vertices(mol)]

function multiplicity(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return Int[multiplicity(props(mol, i)) for i in vertices(mol)]
end


"""
    isotope(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing isotope numbers (or 0 if not specified) of 1 to ``n``th atoms of
the given molecule.
"""
isotope(mol::SimpleMolGraph) = Int[isotope(props(mol, i)) for i in vertices(mol)]

function isotope(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return Int[isotope(props(mol, i)) for i in vertices(mol)]
end


"""
    bond_order(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing bond order of 1 to ``n``th bonds of
the given molecule.
"""
bond_order(mol::SimpleMolGraph) = Int[bond_order(props(mol, e)) for e in edges(mol)]

function bond_order(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    if has_descriptor(mol, :bond_order)
        return get_descriptor(mol, :bond_order)
    end
    return Int[bond_order(props(mol, e)) for e in edges(mol)]
end

default_bond_order!(mol::ReactiveMolGraph) = setproperty!(
    mol.gprops.descriptors, :bond_order,
    Int[bond_order(props(mol, e)) for e in edges(mol)]
)

# coords -> src/coords.jl


# Secondary properties (descriptors)

# Valence

"""
    apparent_valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total bond order
incident to 1 to ``n``th atoms of the given molecule.
"""
apparent_valence(mol::SimpleMolGraph) = apparent_valence(mol.graph, bond_order(mol))

function apparent_valence(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    if has_descriptor(mol, :apparent_valence)
        return get_descriptor(mol, :apparent_valence)
    end
    return apparent_valence(mol.graph, bond_order(mol))
end

apparent_valence!(mol::ReactiveMolGraph) = setproperty!(
    mol.gprops.descriptors, :apparent_valence,
    apparent_valence(mol.graph, bond_order(mol))
)

function apparent_valence(g::SimpleGraph, order_arr::Vector{Int})
    arr = fill(zero(Int), nv(g))
    er = Dict(e => i for (i, e) in enumerate(edges(g)))
    for e in edges(g)
        arr[src(e)] += order_arr[er[e]]
        arr[dst(e)] += order_arr[er[e]]
    end
    return arr
end


"""
    valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the intrinsic valence of
1 to ``n``th atoms of the given molecule.

The number of implicit hydrogens would be calculated based on the valence.
The valence of a hypervalent atom or a non-organic atom is the same as
its `apparent_valence`. This property corresponds to SMARTS `v` query.
"""
valence(mol::SimpleMolGraph
    ) = valence(atom_symbol(mol), atom_charge(mol), apparent_valence(mol))

function valence(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    if has_descriptor(mol, :valence)
        return get_descriptor(mol, :valence)
    end
    return valence(atom_symbol(mol), atom_charge(mol), apparent_valence(mol))
end

valence!(mol::ReactiveMolGraph) = setproperty!(
    mol.gprops.descriptors, :valence,
    valence(atom_symbol(mol), atom_charge(mol), apparent_valence(mol))
)

function valence(
        symbol_arr::Vector{Symbol}, charge_arr::Vector{Int},
        apparent_valence_arr::Vector{Int})
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


"""
    explicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of explicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Explicit" means hydrogens are explicitly represented as graph nodes.
"""
explicit_hydrogens(mol::SimpleMolGraph) = explicit_hydrogens(mol.graph, atom_symbol(mol))

function explicit_hydrogens(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return explicit_hydrogens(mol.graph, atom_symbol(mol))
end

function explicit_hydrogens(g::SimpleGraph, symbol_arr::Vector{Symbol})
    arr = fill(zero(Int), nv(g))
    for i in vertices(g)
        symbol_arr[i] === :H || continue
        for nbr in neighbors(g, i)
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
    return valence(mol) - apparent_valence(mol)
end


"""
    total_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total hydrogens
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `H` query.
"""
function total_hydrogens(mol::SimpleMolGraph)
    return explicit_hydrogens(mol) + implicit_hydrogens(mol)
end


"""
    connectivity(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total atoms
(implicit and explicit) connected to 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `X` query.
"""
function connectivity(mol::SimpleMolGraph)
    return degree(mol) + implicit_hydrogens(mol)
end


"""
    lone_pair(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of lone pairs of
1 to ``n``th atoms of the given molecule.
"""
lone_pair(mol::SimpleMolGraph
    ) = lone_pair(atom_symbol(mol), atom_charge(mol), connectivity(mol))

function lone_pair(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    if has_descriptor(mol, :lone_pair)
        return get_descriptor(mol, :lone_pair)
    end
    return lone_pair(atom_symbol(mol), atom_charge(mol), connectivity(mol))
end

lone_pair!(mol::ReactiveMolGraph) = setproperty!(
    mol.gprops.descriptors, :lone_pair,
    lone_pair(atom_symbol(mol), atom_charge(mol), connectivity(mol))
)

function lone_pair(
        symbol_arr::Vector{Symbol}, charge_arr::Vector{Int}, connectivity_arr::Vector{Int})
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



# Hydrogen bond donor/acceptor

"""
    is_hydrogen_acceptor(mol::SimpleMolGraph) -> Vector{Bool}
    is_hydrogen_acceptor(
        atom_symbol::Symbol, lone_pair_count::Int) -> Bool

Return whether the atom is a hydrogen bond acceptor (N, O and F).
"""
is_hydrogen_acceptor(mol::SimpleMolGraph
    ) = is_hydrogen_acceptor.(atom_symbol(mol), lone_pair(mol))

is_hydrogen_acceptor(symbol::Symbol, lone_pairs::Int
    ) = symbol in HYDROGEN_ACCEPTOR_ATOMS && lone_pairs > 0


"""
    hydrogen_acceptor_count(mol::SimpleMolGraph) -> Int

Return the total number of hydrogen bond acceptors (N, O and F).
"""
hydrogen_acceptor_count(mol::SimpleMolGraph
    ) = hydrogen_acceptor_count(mol, vproptype(mol), eproptype(mol))

hydrogen_acceptor_count(
    mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond}
) = sum(is_hydrogen_acceptor(mol))

function hydrogen_acceptor_count(
        mol::SimpleMolGraph, ::Type{<:AbstractAtom}, ::Type{<:AbstractBond})
    cnt = sum(is_hydrogen_acceptor(mol))
    vcnt = sum(hydrogen_acceptor_count(props(mol, i)) for i in vertices(mol))
    return cnt + vcnt
end

hydrogen_acceptor_count(atom, ::Val{true}) = hydrogen_acceptor_count(group.mol)
hydrogen_acceptor_count(atom, ::Val{false}) = 0
hydrogen_acceptor_count(atom::AbstractAtom
    ) = hydrogen_acceptor_count(atom, Val(has_mol(typeof(atom))))


"""
    is_hydrogen_donor(mol::SimpleMolGraph) -> Vector{Bool}
    is_hydrogen_donor(
        atom_symbol::Symbol, hydrogen_count::Int) -> Bool

Return whether the atom is a hydrogen bond donor (O and N attached to hydrogens).
"""
is_hydrogen_donor(mol::SimpleMolGraph
    ) = is_hydrogen_donor.(atom_symbol(mol), total_hydrogens(mol))

is_hydrogen_donor(symbol::Symbol, total_hydrogens::Int
    ) = symbol in HYDROGEN_DONOR_ATOMS && total_hydrogens > 0


"""
    hydrogen_donor_count(mol::SimpleMolGraph) -> Int

Return the total number of hydrogen bond donors (O and N attached to hydrogens).
"""
hydrogen_donor_count(mol::SimpleMolGraph
    ) = hydrogen_donor_count(mol, vproptype(mol), eproptype(mol))

hydrogen_donor_count(
    mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond}
) = sum(is_hydrogen_donor(mol))

function hydrogen_donor_count(
        mol::SimpleMolGraph, ::Type{<:AbstractAtom}, ::Type{<:AbstractBond})
    cnt = sum(is_hydrogen_donor(mol))
    vcnt = sum(hydrogen_donor_count(props(mol, i)) for i in vertices(mol))
    return cnt + vcnt
end

hydrogen_donor_count(atom, ::Val{true}) = hydrogen_donor_count(group.mol)
hydrogen_donor_count(atom, ::Val{false}) = 0
hydrogen_donor_count(atom::AbstractAtom
    ) = hydrogen_donor_count(atom, Val(has_mol(typeof(atom))))



# Rotatable bonds
# TODO: options
# consider conjugation - e.g. amide and esters

"""
    is_rotatable(mol::SimpleMolGraph) -> Vector{Bool}

Return a vector of size ``n`` representing whether 1 to ``n``th bonds
of the given molecule are rotatable or not.
"""
is_rotatable(mol::SimpleMolGraph) = is_rotatable(mol, vproptype(mol), eproptype(mol))

function is_rotatable(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return is_rotatable(mol, vproptype(mol), eproptype(mol))
end

function is_rotatable(mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond})
    srcdeg = [degree(mol, e.src) for e in edges(mol)]
    dstdeg = [degree(mol, e.dst) for e in edges(mol)]
    return is_rotatable(srcdeg, dstdeg, is_edge_in_ring(mol), bond_order(mol))
end

function is_rotatable(mol::SimpleMolGraph, ::Type{<:AbstractAtom}, ::Type{<:AbstractBond})
    srcdeg = Int[]
    dstdeg = Int[]
    for e in edges(mol)
        ep = props(mol, e)
        if has_submap(typeof(ep))
            src = props(mol, e.src)
            if has_mol(typeof(src))
                snbr = only(neighbors(src.mol, ep.src[2]))
                push!(srcdeg, degree(src.mol, snbr))
            else
                push!(srcdeg, degree(mol, e.src))
            end
            dst = props(mol, e.dst)
            if has_mol(typeof(dst))
                dnbr = only(neighbors(dst.mol, ep.dst[2]))
                push!(dstdeg, degree(dst.mol, dnbr))
            else
                push!(dstdeg, degree(mol, e.dst))
            end
        else
            push!(srcdeg, degree(mol, e.src))
            push!(dstdeg, degree(mol, e.dst))
        end
    end
    return is_rotatable(srcdeg, dstdeg, is_edge_in_ring(mol), bond_order(mol))
end

is_rotatable(
    srcdeg::Int, dstdeg::Int, edge_in_ring::Bool, order::Int
) = ~edge_in_ring & (order == 1) & (srcdeg != 1) & (dstdeg != 1)

is_rotatable(
    srcdeg_arr::Vector{Int}, dstdeg_arr::Vector{Int},
    edge_in_ring_arr::BitVector, order_arr::Vector{Int}
) = is_rotatable.(srcdeg_arr, dstdeg_arr, edge_in_ring_arr, order_arr)


"""
    rotatable_count(mol::SimpleMolGraph) -> Int

Return the total number of rotatable bonds.
"""
rotatable_count(mol::SimpleMolGraph) = rotatable_count(mol, vproptype(mol), eproptype(mol))

rotatable_count(
    mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond}
) = sum(is_rotatable(mol))

function rotatable_count(
        mol::SimpleMolGraph, ::Type{<:AbstractAtom}, ::Type{<:AbstractBond})
    cnt = sum(is_rotatable(mol))
    vcnt = sum(rotatable_count(props(mol, i)) for i in vertices(mol))
    return cnt + vcnt
end

rotatable_count(atom::AbstractAtom
    ) = rotatable_count(atom, Val(has_mol(typeof(atom))))
rotatable_count(atom, ::Val{false}) = 0
rotatable_count(atom, ::Val{true}) = rotatable_count(atom.mol)

# Composition

"""
    atom_counter(mol::SimpleMolGraph) -> Dict{Symbol,Int}
    atom_counter(
        symbol_arr::Vector{Symbol}, implh_arr::Vector{Int}) -> Dict{Symbol,Int}
    atom_counter(atom::AbstractAtom) -> Dict{Symbol,Int}

Count the number of atoms and return symbol => count dict.
"""
atom_counter(mol::SimpleMolGraph) = atom_counter(mol, vproptype(mol), eproptype(mol))

atom_counter(
    mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond}
) = atom_counter(atom_symbol(mol), implicit_hydrogens(mol))

function _inc!(counter::Dict{Symbol,Int}, sym::Symbol, cnt=1)
    if !haskey(counter, sym)
        counter[sym] = 0
    end
    counter[sym] += cnt
end

function atom_counter(mol::SimpleMolGraph, ::Type{<:AbstractAtom}, ::Type{<:AbstractBond})
    counter = atom_counter(atom_symbol(mol), implicit_hydrogens(mol))
    # count groups
    for i in vertices(mol)
        for (sym, cnt) in atom_counter(props(mol, i))
            _inc!(counter, sym, cnt)
        end
    end
    return counter
end

function atom_counter(symbol_arr::Vector{Symbol}, implh_arr::Vector{Int})
    counter = Dict{Symbol,Int}()
    for sym in symbol_arr
        haskey(ATOMSYMBOLMAP, sym) || continue
        _inc!(counter, sym)
    end
    _inc!(counter, :H, sum(implh_arr))
    return counter
end

atom_counter(atom::AbstractAtom) = atom_counter(
    atom, 
    Val(has_mol(typeof(atom))),
    Val(has_formula(typeof(atom))),
    Val(has_hydrogens(typeof(atom))),
    Val(has_label(typeof(atom)))
)
atom_counter(atom, ::Val{false}, ::Val{true}, ::Val, ::Val) = atom.formula
atom_counter(atom, ::Val{false}, ::Val{false}, ::Val{true}, ::Val
    ) = Dict(atom_symbol(atom.center) => 1, :H => atom.hydrogens)
atom_counter(atom, ::Val{false}, ::Val{false}, ::Val{false}, ::Val{true}
    ) = Dict(atom_symbol(atom) => 1)
atom_counter(atom, ::Val{false}, ::Val{false}, ::Val{false}, ::Val{false}
    ) = Dict{Symbol,Int}()

function atom_counter(atom, ::Val{true}, ::Val, ::Val, ::Val)
    cnt = atom_counter(atom.mol)
    sym = atom_symbol(atom.mol)
    implh = implicit_hydrogens(atom.mol)
    for (_, v) in atom.term  # remove terminal atoms
        cnt[sym[v]] -= 1
        cnt[:H] -= implh[v]
    end
    return cnt
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



"""
    chargesign(charge::Int) -> String

Get a charge sign.
"""
function chargesign(charge::Int)
    charge == 0 && return ""
    sign = charge > 0 ? "+" : "–" # en dash, not hyphen-minus
    num = abs(charge)
    return num > 1 ? string(num, sign) : sign
end

function sanitize_html(s::AbstractString)
    s = replace(s, "&" => "&amp;")
    s = replace(s, "<" => "&lt;")
    s = replace(s, ">" => "&gt;")
    s = replace(s, "\"" => "&quot;")
    s = replace(s, "'" => "&#39;")
    return s
end

function sanitize_markup(markup::Vector{Vector{Tuple{Symbol,String}}})
    san = Vector{Tuple{Symbol,String}}[]
    for gp in markup
        push!(san, [(sym, sanitize_html(unsafe_str)) for (sym, unsafe_str) in gp])
    end
    return san
end


"""
    atom_markup(mol::SimpleMolGraph) -> Vector{Vector{Tuple{Symbol,String}}}
    atom_markup(
        symbol_arr::Vector{Symbol}, charge_arr::Vector{Int}, implh_arr::Vector{Int}
    ) -> Vector{Vector{Tuple{Symbol,String}}}
    atom_markup(atom::AbstractAtom) -> Vector{Vector{Tuple{Symbol,String}}}

Return atom labels with style annotations which is used for drawing 2D structure and writing formula
"""
atom_markup(mol::SimpleMolGraph) = atom_markup(mol, vproptype(mol), eproptype(mol))

function atom_markup(mol::SimpleMolGraph, ::Type{<:AbstractAtom}, ::Type{<:AbstractBond})
    markup = atom_markup(atom_symbol(mol), isotope(mol), implicit_hydrogens(mol))
    for i in vertices(mol)
        markup[i] = something(atom_markup(props(mol, i)), markup[i])
    end
    return markup
end

atom_markup(
    mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond}
) = atom_markup(atom_symbol(mol), isotope(mol), implicit_hydrogens(mol))

function atom_markup(sym::Symbol, isotope::Int, implicith::Int)
    contents = [[(:default, string(sym))]]
    if isotope > 0
        pushfirst!(contents[1], (:sup, string(isotope)))
    end
    if implicith == 1
        push!(contents, [(:default, "H")])
    elseif implicith > 1
        push!(contents, [(:default, "H"), (:sub, string(implicith))])
    end
    return contents
end

atom_markup(
    symbol_arr::Vector{Symbol}, isotope_arr::Vector{Int}, implh_arr::Vector{Int}
) = atom_markup.(symbol_arr, isotope_arr, implh_arr)

atom_markup(atom::AbstractAtom) = atom_markup(
    atom,
    Val(has_formula(typeof(atom))),
    Val(has_hydrogens(typeof(atom))),
    Val(has_label(typeof(atom)))
)
atom_markup(atom, ::Val{true}, ::Val{false}, ::Val{true}
    ) = isempty(atom.label) ? markup_formula(atom.formula) : atom.label
atom_markup(atom, ::Val{false}, ::Val{true}, ::Val{false}
    ) = [[(:default, string(atom_symbol(atom.center))), (:sub, string(atom.hydrogens))]]
atom_markup(atom, ::Val{false}, ::Val{false}, ::Val{true}) = atom.label
atom_markup(atom, ::Val{false}, ::Val{false}, ::Val{false}) = nothing



"""
    markup_formula(
        counter::Dict{Symbol,Int}; charge=0) -> Vector{Vector{Tuple{Symbol,String}}}

Generate markup formula array
"""
function markup_formula(counter_::Dict{Symbol,Int})
    contents = Vector{Tuple{Symbol,String}}[]
    counter = copy(counter_)
    if haskey(counter, :C)
        elems = [(:default, "C")]
        c = pop!(counter, :C)
        c > 1 && push!(elems, (:sub, string(c)))
        push!(contents, elems)
        if haskey(counter, :H)
            elems = [(:default, "H")]
            h = pop!(counter, :H)
            h > 1 && push!(elems, (:sub, string(h)))
            push!(contents, elems)
        end
    end
    for sym in sort(collect(keys(counter)))
        elems = [(:default, string(sym))]
        cnt = counter[sym]
        cnt > 1 && push!(elems, (:sub, string(cnt)))
        push!(contents, elems)
    end
    return contents
end


write_formula(markup::Vector{Vector{Tuple{Symbol,String}}}) = join(e[2] for e in vcat(markup...))
write_formula(counter::Dict{Symbol,Int}) = write_formula(markup_formula(counter))


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
