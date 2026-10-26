#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# Number of lone pairs of the neutral atom with the typical valence (octet).
# The number of valence electrons is `4 + LONEPAIR_COUNT[sym]`.
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
atom_symbol(mol::SimpleMolGraph) = Symbol[atom_symbol(mol[i]) for i in vertices(mol)]


"""
    atom_number(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom numbers of 1 to ``n``th atoms of
the given molecule.
"""
atom_number(mol::SimpleMolGraph) = Int[atom_number(mol[i]) for i in vertices(mol)]


"""
    atom_charge(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom charges of 1 to ``n``th atoms of
the given molecule.
"""
atom_charge(mol::SimpleMolGraph) = Int[atom_charge(mol[i]) for i in vertices(mol)]
function atom_charge(mol::ReactiveMolGraph)
    if has_descriptor(mol, Val(:atom_charge))
        return get_descriptor(mol, Val(:atom_charge))
    end
    return Int[atom_charge(mol[i]) for i in vertices(mol)]
end

"""
    default_atom_charge!(mol::SimpleMolGraph) -> Nothing

Set atom charges of the atom objects to the `atom_charge` descriptor (see [`atom_charge`](@ref)).
Called before preprocessing that modifies charges (e.g. in [`smiles_on_update!`](@ref)).
"""
default_atom_charge!(mol::SimpleMolGraph) = set_descriptor!(
    mol, Val(:atom_charge), Int[atom_charge(mol[i]) for i in vertices(mol)])


"""
    atom_valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the valences specified to 1 to ``n``th atoms of
the given molecule (`-1` if not specified). See [`valence`](@ref) for the actual valences.
"""
atom_valence(mol::SimpleMolGraph) = Int[atom_valence(mol[i]) for i in vertices(mol)]


"""
    multiplicity(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing atom multiplicities of 1 to ``n``th atoms of
the given molecule (1: non-radical, 2: radical, 3: biradical).
"""
multiplicity(mol::SimpleMolGraph) = Int[multiplicity(mol[i]) for i in vertices(mol)]


"""
    isotope(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing isotope numbers (or 0 if not specified) of 1 to ``n``th atoms of
the given molecule.
"""
isotope(mol::SimpleMolGraph) = Int[isotope(mol[i]) for i in vertices(mol)]


"""
    bond_order(mol::MolGraph) -> Vector{Int}

Return a vector of size ``n`` representing bond order of 1 to ``n``th bonds of
the given molecule.
"""
bond_order(mol::SimpleMolGraph) = Int[bond_order(mol[e]) for e in edges(mol)]
function bond_order(mol::ReactiveMolGraph)
    if has_descriptor(mol, Val(:bond_order))
        return get_descriptor(mol, Val(:bond_order))
    end
    return Int[bond_order(mol[e]) for e in edges(mol)]
end

"""
    default_bond_order!(mol::SimpleMolGraph) -> Nothing

Set bond orders of the bond objects to the `bond_order` descriptor (see [`bond_order`](@ref)).
Called before preprocessing that modifies bond orders (e.g. [`kekulize!`](@ref) in
[`smiles_on_update!`](@ref)).
"""
default_bond_order!(mol::SimpleMolGraph) = set_descriptor!(
    mol, Val(:bond_order), Int[bond_order(mol[e]) for e in edges(mol)])


# coords -> src/coords.jl


# Secondary properties (descriptors)

# Valence

"""
    apparent_valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of total bond order
incident to 1 to ``n``th atoms of the given molecule.

Before kekulization, an aromatic atom (`isaromatic` of the atom type, e.g. SMILES lowercase
atoms) counts its aromatic bonds as single bonds plus one ``\\pi`` bond, e.g. 3 for `c` and
`n` of pyridine (this overestimates `[nH]` and `o`, which have no ``\\pi`` bond).
[`kekulize!`](@ref) overwrites it with the total order of the kekulized bonds.
"""
function apparent_valence(mol::SimpleMolGraph)
    if has_descriptor(mol, Val(:apparent_valence))
        return get_descriptor(mol, Val(:apparent_valence))
    end
    return apparent_valence(mol.graph, bond_order(mol), atom_isaromatic(mol))
end

apparent_valence(g::SimpleGraph, order_arr::Vector{Int}
    ) = bond_order_sum(nv(g), edges(g), order_arr)

apparent_valence(g::SimpleGraph, order_arr::Vector{Int}, arom::BitVector
    ) = apparent_valence(g, order_arr) + arom

"""
    apparent_valence!(mol::SimpleMolGraph)

Calculate [`apparent_valence`](@ref) and store the result in the descriptors of `mol`.
"""
apparent_valence!(mol::SimpleMolGraph) = set_descriptor!(
    mol, Val(:apparent_valence),
    apparent_valence(mol.graph, bond_order(mol), atom_isaromatic(mol)))

# `isaromatic` of the atoms (false if the atom type does not have it)
atom_isaromatic(mol::SimpleMolGraph) = BitVector(
    has_isaromatic(typeof(mol[i])) && isaromatic(mol[i]) for i in vertices(mol))

# Edge list versions (also used by the parsers before the graph is constructed)

# Total bond order of each atom
function bond_order_sum(n::Int, edges, order_arr::Vector{Int})
    arr = zeros(Int, n)
    for (e, o) in zip(edges, order_arr)
        arr[src(e)] += o
        arr[dst(e)] += o
    end
    return arr
end


# Typical valence of the element (octet rule), 0 if the element is not supported
formal_valence(sym::Symbol, charge::Int
    ) = haskey(LONEPAIR_COUNT, sym) ? 4 - abs(LONEPAIR_COUNT[sym] - charge) : 0

formal_valence(symbol_arr::Vector{Symbol}, charge_arr::Vector{Int}
    ) = formal_valence.(symbol_arr, charge_arr)

# Valence of the atom with `hcount` hydrogens that are not graph nodes (e.g. CDXML
# NumHydrogens) and the total bond order `bondsum`. Returns -1 (not specified) if it is the
# same as the valence inferred from the element. Elements other than typical organic atoms
# (except hydrogen) always have the specified valence (e.g. 0 for `[Sn]`), because their
# valence is not inferred from the element (other software may add hydrogens).
function hydrogen_count_valence(sym::Symbol, charge::Int, bondsum::Int, hcount::Int)
    v = bondsum + hcount
    (haskey(LONEPAIR_COUNT, sym) || sym === :H) || return v
    return v == max(formal_valence(sym, charge), bondsum) ? -1 : v
end

# Specified valences of atoms (-1: not specified) from the number of hydrogens that are not
# graph nodes (`hcount_arr`, -1: not specified, e.g. 0 for SMILES bracket atoms) and the bonds
# (the same as `apparent_valence` before kekulization).
function hydrogen_count_valence(
        edges, order_arr::Vector{Int}, arom::BitVector,
        symbol_arr::Vector{Symbol}, charge_arr::Vector{Int}, hcount_arr::Vector{Int})
    n = length(arom)
    bondsum = bond_order_sum(n, edges, order_arr) + arom
    return Int[hcount_arr[i] == -1 ? -1 : hydrogen_count_valence(
        symbol_arr[i], charge_arr[i], bondsum[i], hcount_arr[i]) for i in 1:n]
end

# Specified valence if available, otherwise the typical valence
expected_valence(symbol_arr::Vector{Symbol}, charge_arr::Vector{Int}, atomval_arr::Vector{Int}
    ) = [v >= 0 ? v : formal_valence(sym, chg) for (sym, chg, v) in zip(symbol_arr, charge_arr, atomval_arr)]


"""
    valence(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the intrinsic valence of
1 to ``n``th atoms of the given molecule.

The number of implicit hydrogens would be calculated based on the valence.
If the valence is specified to the atom (see [`atom_valence`](@ref), e.g. carbenes and
pentavalent nitrogens), the specified valence is used. Otherwise, the typical valence of
organic atoms is used. The valence of a hypervalent atom or a non-organic atom is the same as
its `apparent_valence`. This property corresponds to SMARTS `v` query.
"""
function valence(mol::SimpleMolGraph)
    if has_descriptor(mol, Val(:valence))
        return get_descriptor(mol, Val(:valence))
    end
    return valence(atom_symbol(mol), atom_charge(mol), atom_valence(mol), apparent_valence(mol))
end

valence(
    symbol_arr::Vector{Symbol}, charge_arr::Vector{Int}, apparent_valence_arr::Vector{Int}
) = max.(formal_valence(symbol_arr, charge_arr), apparent_valence_arr)

valence(
    symbol_arr::Vector{Symbol}, charge_arr::Vector{Int},
    atomval_arr::Vector{Int}, apparent_valence_arr::Vector{Int}
) = max.(expected_valence(symbol_arr, charge_arr, atomval_arr), apparent_valence_arr)

"""
    valence!(mol::SimpleMolGraph)

Calculate [`valence`](@ref) and store the result in the descriptors of `mol`.
"""
valence!(mol::SimpleMolGraph) = set_descriptor!(
    mol, Val(:valence),
    valence(atom_symbol(mol), atom_charge(mol), atom_valence(mol), apparent_valence(mol)))


function check_valence!(mol::SimpleMolGraph)
    av = apparent_valence(mol)
    fv = expected_valence(atom_symbol(mol), atom_charge(mol), atom_valence(mol))
    comments = String[]
    for i in vertices(mol)[av .> fv]
         push!(comments, "unusual valence at #$(i) (actual $(av[i]), expected $(fv[i]))")
    end
    mol[:logs]["warning_valence"] = join(comments, "; ")
end


"""
    explicit_hydrogens(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of explicit hydrogens
connected to 1 to ``n``th atoms of the given molecule.

"Explicit" means hydrogens are explicitly represented as graph nodes.
"""
explicit_hydrogens(mol::SimpleMolGraph) = explicit_hydrogens(mol.graph, atom_symbol(mol))

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
    lone_pair(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of lone pairs of
1 to ``n``th atoms of the given molecule.

The number of lone pairs is calculated from the number of valence electrons, the charge,
the [`valence`](@ref) and the unpaired electrons ([`multiplicity`](@ref) - 1), e.g. 1 for
singlet carbenes, 0 for triplet carbenes and pentavalent nitrogens. Elements other than
typical organic atoms have no lone pairs.
"""
function lone_pair(mol::SimpleMolGraph)
    if has_descriptor(mol, Val(:lone_pair))
        return get_descriptor(mol, Val(:lone_pair))
    end
    return lone_pair(atom_symbol(mol), atom_charge(mol), multiplicity(mol), valence(mol))
end

function lone_pair(sym::Symbol, charge::Int, multiplicity::Int, valence::Int)
    haskey(LONEPAIR_COUNT, sym) || return 0
    nonbonding = 4 + LONEPAIR_COUNT[sym] - charge - valence - (multiplicity - 1)
    return max(fld(nonbonding, 2), 0)
end

lone_pair(
    symbol_arr::Vector{Symbol}, charge_arr::Vector{Int},
    multiplicity_arr::Vector{Int}, valence_arr::Vector{Int}
) = lone_pair.(symbol_arr, charge_arr, multiplicity_arr, valence_arr)

"""
    lone_pair!(mol::SimpleMolGraph)

Calculate [`lone_pair`](@ref) and store the result in the descriptors of `mol`.
"""
lone_pair!(mol::SimpleMolGraph) = set_descriptor!(
    mol, Val(:lone_pair),
    lone_pair(atom_symbol(mol), atom_charge(mol), multiplicity(mol), valence(mol)))




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
    vcnt = sum(hydrogen_acceptor_count(mol[i]) for i in vertices(mol))
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
    vcnt = sum(hydrogen_donor_count(mol[i]) for i in vertices(mol))
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

function is_rotatable(mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond})
    srcdeg = [degree(mol, e.src) for e in edges(mol)]
    dstdeg = [degree(mol, e.dst) for e in edges(mol)]
    return is_rotatable(srcdeg, dstdeg, is_edge_in_ring(mol), bond_order(mol))
end

function is_rotatable(mol::SimpleMolGraph, ::Type{<:AbstractAtom}, ::Type{<:AbstractBond})
    srcdeg = Int[]
    dstdeg = Int[]
    for e in edges(mol)
        ep = mol[e]
        if has_submap(typeof(ep))
            src = mol[e.src]
            if has_mol(typeof(src))
                snbr = only(neighbors(src.mol, ep.src[2]))
                push!(srcdeg, degree(src.mol, snbr))
            else
                push!(srcdeg, degree(mol, e.src))
            end
            dst = mol[e.dst]
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
    vcnt = sum(rotatable_count(mol[i]) for i in vertices(mol))
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
        for (sym, cnt) in atom_counter(mol[i])
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
    implh = sum(implh_arr; init=0)
    implh > 0 && _inc!(counter, :H, implh)  # no H entry (count 0) for e.g. [Na+].[Cl-]
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
        implh[v] == 0 || (cnt[:H] -= implh[v])
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
    atom_markup(mol::SimpleMolGraph) -> Vector{Vector{Vector{Tuple{Symbol,String}}}}
    atom_markup(
        symbol_arr::Vector{Symbol}, isotope_arr::Vector{Int}, implh_arr::Vector{Int}
    ) -> Vector{Vector{Vector{Tuple{Symbol,String}}}}
    atom_markup(atom::AbstractAtom) -> Union{Vector{Vector{Tuple{Symbol,String}}},Nothing}

Return atom labels with style annotations which is used for drawing 2D structure and writing formula
(one markup for each atom). The single atom form returns `nothing` if the atom type
does not have its own label.
"""
atom_markup(mol::SimpleMolGraph) = atom_markup(mol, vproptype(mol), eproptype(mol))

function atom_markup(mol::SimpleMolGraph, ::Type{<:AbstractAtom}, ::Type{<:AbstractBond})
    markup = atom_markup(atom_symbol(mol), isotope(mol), implicit_hydrogens(mol))
    for i in vertices(mol)
        markup[i] = something(atom_markup(mol[i]), markup[i])
    end
    return markup
end

atom_markup(
    mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond}
) = atom_markup(atom_symbol(mol), isotope(mol), implicit_hydrogens(mol))

function atom_markup(sym::Symbol, isotope::Int, implicith::Int)
    contents = [[(:default, string(sym))]]
    # D and T already imply their mass number
    if isotope > 0 && isotope != hydrogen_isotope(sym)
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
    markup_formula(counter::Dict{Symbol,Int}) -> Vector{Vector{Tuple{Symbol,String}}}

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


write_formula(markup::Vector{Vector{Tuple{Symbol,String}}}
    ) = join(e[2] for e in Iterators.flatten(markup))
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
