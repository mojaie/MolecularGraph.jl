#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const STEREOCENTER_STATE = Dict(
    (1, 2, 3) => true, (1, 2, 4) => false,
    (1, 3, 2) => false, (1, 3, 4) => true,
    (1, 4, 2) => true, (1, 4, 3) => false,
    (2, 1, 3) => false, (2, 1, 4) => true,
    (2, 3, 1) => true, (2, 3, 4) => false,
    (2, 4, 1) => false, (2, 4, 3) => true,
    (3, 1, 2) => true, (3, 1, 4) => false,
    (3, 2, 1) => false, (3, 2, 4) => true,
    (3, 4, 1) => true, (3, 4, 2) => false,
    (4, 1, 2) => false, (4, 1, 3) => true,
    (4, 2, 1) => true, (4, 2, 3) => false,
    (4, 3, 1) => false, (4, 3, 2) => true
)


stereoneighbors(c::Stereocenter) = [c.lookingFrom, c.first, c.second]
isclockwise(c::Stereocenter) = c.isclockwise


function isclockwise(c::Stereocenter{T}, f::T, s::T, t::T) where T <: Integer
    fp = something(findfirst(==(f), stereoneighbors(c)), 4)
    sp = something(findfirst(==(s), stereoneighbors(c)), 4)
    tp = something(findfirst(==(t), stereoneighbors(c)), 4)
    return isclockwise(c) === STEREOCENTER_STATE[(fp, sp, tp)]
end



"""
    set_stereocenter!(mol::SimpleMolGraph, center, looking_from, v1, v2, is_clockwise) -> Nothing

Set stereocenter information to graph properties.
"""
function set_stereocenter!(
        mol::SimpleMolGraph{T}, center::T,
        looking_from::T, v1::T, v2::T, is_clockwise::Bool) where T
    mol[:stereocenter][center] = Stereocenter{T}(looking_from, v1, v2, is_clockwise)
    return
end


"""
    set_stereobond!(mol::SimpleMolGraph, bond::Edge, v1, v2, is_cis::Bool) -> Nothing

Set stereobond (cis-trans) information to graph properties.
"""
function set_stereobond!(
        mol::SimpleMolGraph{T}, bond::Edge{T}, v1::T, v2::T, is_cis::Bool) where T
    mol[:stereobond][bond] = Stereobond{T}(v1, v2, is_cis)
    return
end


"""
    stereo_hydrogen(mol::SimpleMolGraph, v::Integer) -> Union{Int,Nothing}

Return the explicit hydrogen node adjacent to the stereocenter `v`, or `nothing` if there is no
hydrogen neighbor (e.g. quaternary center or the hydrogen is already removed).
"""
function stereo_hydrogen(mol::SimpleMolGraph, v::Integer)
    nbrs = neighbors(mol, v)
    hpos = findfirst(x -> atom_symbol(mol[x]) === :H, nbrs)
    isnothing(hpos) && return  # 4° center or already removed
    return nbrs[hpos]
end

function stereo_hydrogen(mol::ReactiveMolGraph{T,V,E}, v::Integer
        ) where {T,V<:QueryTree,E<:QueryTree}
    nbrs = neighbors(mol, v)
    hpos = findfirst(x -> mol[x].vprops[1] == qeq(:symbol, "H"), nbrs)
    isnothing(hpos) && return  # 4° center or already removed
    return nbrs[hpos]
end


"""
    safe_stereo_hydrogen!(mol::SimpleMolGraph{T}, center::T) -> Union{T,Nothing}

Rearrange stereocenter properties to safely remove stereo hydrogen nodes
and return the hydrogen node index (or `nothing` if there is no stereo hydrogen).

This function is called inside `rem_vertex!` and `rem_vertices!` functions
to safely remove hydrogen nodes while preserving stereocenter information.
"""
function safe_stereo_hydrogen!(mol::SimpleMolGraph{T}, center::T) where T
    h = stereo_hydrogen(mol, center)
    isnothing(h) && return # 4° center or already removed
    stereo = mol[:stereocenter][center]
    spos = findfirst(==(h), stereoneighbors(stereo))
    isnothing(spos) && return h  # hydrogen at the lowest priority can be removed safely
    nonh = setdiff(ordered_neighbors(mol, center), h)
    f, s, t = nonh[1], nonh[2], nonh[3]  # three non-hydrogen neighbors
    set_stereocenter!(mol, center, f, s, t, isclockwise(stereo, f, s, t))
    return h
end



function angeval(u::Point2d, v::Point2d)
    # 0deg -> 1, 90deg -> 0, 180deg -> -1, 270deg-> -2, 360deg -> -3
    uv = dot(u, v) / (norm(u) * norm(v))
    return cross(u, v) >= 0 ? uv : -2 - uv
end


function anglesort(
        coords::Coords2d, center::T, ref::T,
        vertices::Vector{T}) where T <: Integer
    # return vertices order by clockwise direction
    c = coords[center]
    r = coords[ref]
    ps = [coords[v] for v in vertices]
    vs = [p - c for p in ps]
    return sortperm([angeval(r - c, v) for v in vs])
end


# Chirality of a three-coordinate center (with an implicit hydrogen) from the 2D coordinates
# and wedges (directions) of the neighbors. Each wedge gives the neighbor a z coordinate
# (:up => 1, :down => -1), and the implicit hydrogen is placed opposite to the neighbors.
# Returns whether `nbrs[2]` and `nbrs[3]` are clockwise looking from `nbrs[1]`, or `nothing`
# if the geometry is degenerate. This is consistent with the interpretation of InChI also in
# the case that the neighbors are in a half-plane.
function threecoord_isclockwise(
        coords::Coords2d, center::T, nbrs::AbstractVector{T}, drs::Vector{Symbol}) where T
    c = coords[center]
    vs = NTuple{3,Float64}[]
    for (n, dr) in zip(nbrs, drs)
        x, y = coords[n][1] - c[1], coords[n][2] - c[2]
        z = dr === :up ? 1.0 : (dr === :down ? -1.0 : 0.0)
        r = sqrt(x^2 + y^2 + z^2)
        r == 0 && return nothing
        push!(vs, (x / r, y / r, z / r))
    end
    h = (-(vs[1][1] + vs[2][1] + vs[3][1]),
        -(vs[1][2] + vs[2][2] + vs[3][2]), -(vs[1][3] + vs[2][3] + vs[3][3]))
    hr = sqrt(h[1]^2 + h[2]^2 + h[3]^2)
    hr < 1e-9 && return nothing
    h = h ./ hr
    a, b, d = vs[1] .- h, vs[2] .- h, vs[3] .- h
    triple = (a[1] * (b[2] * d[3] - b[3] * d[2]) - a[2] * (b[1] * d[3] - b[3] * d[1])
        + a[3] * (b[1] * d[2] - b[2] * d[1]))
    abs(triple) < 1e-9 && return nothing
    return triple < 0
end


"""
    stereocenter_from_sdf2d(
        g::SimpleGraph{T}, v_symbol::Vector{Symbol}, e_order::Vector{Int}, e_notation::Vector{Int},
        e_isordered::Vector{Bool}, v_coords2d::Coords2d) where T -> Tuple{StereocenterMap{T},Vector{String}}

Return stereocenter information obtained from 2D SDFile and log messages.
"""
function stereocenter_from_sdf2d(
        g::SimpleGraph{T}, v_symbol::Vector{Symbol}, e_order::Vector{Int},
        e_notation::Vector{Int}, e_isordered::Vector{Bool}, v_coords2d::Coords2d) where T
    centers = StereocenterMap{T}()
    ernk = edge_rank(g)
    comments = String[]
    for i in vertices(g)
        degree(g, i) in (3, 4) || continue
        nbrs = ordered_neighbors(g, i)
        drs = Symbol[]  # lookingFrom, atom1, atom2, (atom3)
        for nbr in nbrs
            nbe = edge_rank(ernk, i, nbr)
            e_order[nbe] == 1 || break
            if e_isordered[nbe] == (i < nbr)  # only outgoing wedges are considered
                if e_notation[nbe] == 1
                    push!(drs, :up)
                    continue
                elseif e_notation[nbe] == 6
                    push!(drs, :down)
                    continue
                end
            end
            push!(drs, :unspecified)
        end
        length(drs) == degree(g, i) || (@debug "$(i): multiple bond"; continue)
        degree(g, i) - count(x -> v_symbol[x] === :H, nbrs) > 2 || (@debug "$(i): a stereocenter cannot be connected to multiple hydrogens"; continue)
        upcnt = count(x -> x === :up, drs)
        dwcnt = count(x -> x === :down, drs)
        (upcnt == 0 && dwcnt == 0) && (@debug "$(i): unspecified"; continue)
        sortorder = pushfirst!(anglesort(v_coords2d, i, nbrs[1], nbrs[2:end]) .+ 1, 1)
        ons = nbrs[sortorder]  # node indices of neighbors in angular order
        ods = drs[sortorder]  # direction symbols of neighbors in angular order
        if length(nbrs) == 3
            if upcnt == 1 && dwcnt == 1 || upcnt + dwcnt == 3
                @debug "ignored stereocenter #$(i) (too many wedges)"
                push!(comments, "ignored stereocenter #$(i) (too many wedges)")
                continue
            end
            # The implicit hydrogen may be in the plane (e.g. the neighbors are in a
            # half-plane), so the angular order of the neighbors is not sufficient.
            cw = threecoord_isclockwise(v_coords2d, i, nbrs, drs)
            if isnothing(cw)
                @debug "ignored stereocenter #$(i) (degenerate geometry)"
                push!(comments, "ignored stereocenter #$(i) (degenerate geometry)")
            else
                centers[i] = Stereocenter{T}(nbrs[1], nbrs[2], nbrs[3], cw)
            end
            continue
        end
        if upcnt == 4 || dwcnt == 4
            @debug "ignored stereocenter #$(i) (too many wedges)"
            push!(comments, "ignored stereocenter #$(i) (too many wedges)")
        elseif upcnt == 3
            centers[i] = Stereocenter{T}(ons[1], ons[2], ons[3], ods[1] === :up && ods[3] === :up)
        elseif dwcnt == 3
            centers[i] = Stereocenter{T}(ons[1], ons[2], ons[3], ods[2] === :down && ods[4] === :down)
        elseif upcnt == 2
            if (ods[1] === :up && ods[3] === :up) || (ods[2] === :up && ods[4] === :up)
                centers[i] = Stereocenter{T}(ons[1], ons[2], ons[3], ods[1] === :up)
            else
                @debug "ignored stereocenter #$(i) (adjacent up wedges)"
                push!(comments, "ignored stereocenter #$(i) (adjacent up wedges)")
            end
        elseif dwcnt == 2
            if (ods[1] === :down && ods[3] === :down) || (ods[2] === :down && ods[4] === :down)
                centers[i] = Stereocenter{T}(ons[1], ons[2], ons[3], ods[1] !== :down)
            else
                @debug "ignored stereocenter #$(i) (adjacent down wedges)"
                push!(comments, "ignored stereocenter #$(i) (adjacent down wedges)")
            end
        elseif upcnt == 1 && dwcnt == 1
            if ((ods[1] === :unspecified && ods[3] === :unspecified)
                    || (ods[2] === :unspecified && ods[4] === :unspecified))
                @debug "ignored stereocenter #$(i) (up and down wedges facing each other)"
                push!(comments, "ignored stereocenter #$(i) (up and down wedges facing each other)")
            else
                centers[i] = Stereocenter{T}(ons[1], ons[2], ons[3], ods[1] === :up || ods[3] === :up)
            end
        elseif upcnt == 1
            centers[i] = Stereocenter{T}(ons[1], ons[2], ons[3], ods[1] === :up || ods[3] === :up)
        else  # dwcnt == 1
            centers[i] = Stereocenter{T}(ons[1], ons[2], ons[3], ods[2] === :down || ods[4] === :down)
        end
    end
    return centers, comments
end

stereocenter_from_sdf2d(mol::SimpleMolGraph) = stereocenter_from_sdf2d(
    mol.graph,
    Symbol[mol[i].symbol for i in vertices(mol)],
    Int[mol[e].order for e in edges(mol)],
    Int[mol[e].notation for e in edges(mol)],
    Bool[mol[e].isordered for e in edges(mol)],
    get_descriptor(mol, Val(:coords2d))[1]
)

"""
    stereocenter_from_draw2d(
        g::SimpleGraph{T}, v_symbol::Vector{Symbol}, e_order::Vector{Int},
        styles::Draw2dBondStyle, v_coords2d::Coords2d) where T -> StereocenterMap{T}

Return stereocenter information represented by the 2D bond styles (wedges) and coordinates.
"""
function stereocenter_from_draw2d(
        g::SimpleGraph{T}, v_symbol::Vector{Symbol}, e_order::Vector{Int},
        styles::Draw2dBondStyle, v_coords2d::Coords2d) where T
    # :up/:down start at the source vertex, :revup/:revdown start at the destination vertex
    notation = Int[s in (:up, :revup) ? 1 : s in (:down, :revdown) ? 6 : 0 for s in styles]
    isordered = Bool[s in (:up, :down) for s in styles]
    centers, _ = stereocenter_from_sdf2d(g, v_symbol, e_order, notation, isordered, v_coords2d)
    return centers
end


"""
    assign_stereocenter_wedges!(
        styles::Draw2dBondStyle, g::SimpleGraph{T}, v_symbol::Vector{Symbol},
        e_order::Vector{Int}, coords::Coords2d, centers::StereocenterMap{T}) where T -> Nothing
    assign_stereocenter_wedges!(mol::SimpleMolGraph) -> Nothing

Update 2D bond styles (wedges) so that they represent the stereocenters.

The stereocenters (`mol[:stereocenter]`) are the reference of the stereochemistry, and the
wedges are its 2D notation. For each stereocenter, if the wedges and the 2D coordinates do not
represent the chirality (e.g. the wedges were on removed hydrogens, or the 2D coordinates were
generated without wedges), the wedges starting at the center are replaced with a new wedge.
Wedges that represent the chirality are kept as they are.

The new wedge is put on a single bond, preferably on an acyclic bond to a non-stereocenter
terminal atom (e.g. hydrogen).

For `ReactiveMolGraph`, this is called automatically after the molecule is modified
(see [`dispatch_update!`](@ref)) and in [`coordgen!`](@ref).
"""
function assign_stereocenter_wedges!(
        styles::Draw2dBondStyle, g::SimpleGraph{T}, v_symbol::Vector{Symbol},
        e_order::Vector{Int}, coords::Coords2d, centers::StereocenterMap{T}) where T
    isempty(centers) && return
    ernk = edge_rank(g)
    outgoing(c, n, dr) = c < n ? dr : (dr === :up ? :revup : :revdown)
    is_represented(derived, c, stereo) = haskey(derived, c) && isclockwise(
        derived[c], stereo.lookingFrom, stereo.first, stereo.second) === stereo.isclockwise
    derived = stereocenter_from_draw2d(g, v_symbol, e_order, styles, coords)
    acyclic = Set{Edge{T}}()
    acyclic_ready = false  # calculated only if needed
    for (c, stereo) in centers
        # skip stereocenters inconsistent with the graph (e.g. a neighbor was disconnected)
        nbrs = neighbors(g, c)
        (degree(g, c) in (3, 4) && all(n -> n in nbrs, stereoneighbors(stereo))) || continue
        is_represented(derived, c, stereo) && continue
        if !acyclic_ready
            union!(acyclic, bridges(g))
            acyclic_ready = true
        end
        # discard wedges starting at the center
        for n in nbrs
            i = edge_rank(ernk, c, n)
            styles[i] in (outgoing(c, n, :up), outgoing(c, n, :down)) && (styles[i] = :none)
        end
        cands = sort(collect(nbrs), by=n -> (
            u_edge(T, c, n) ∉ acyclic, haskey(centers, n), degree(g, n)))
        done = false
        for n in cands
            i = edge_rank(ernk, c, n)
            (e_order[i] == 1 && styles[i] === :none) || continue
            for dr in (:up, :down)
                styles[i] = outgoing(c, n, dr)
                done = is_represented(
                    stereocenter_from_draw2d(g, v_symbol, e_order, styles, coords), c, stereo)
                done && break
            end
            done && break
            styles[i] = :none
        end
        done || @debug "$(c): failed to assign stereocenter wedges"
    end
    return
end

function assign_stereocenter_wedges!(mol::SimpleMolGraph)
    # gprops are accessed directly so as not to dispatch updates
    has_descriptor(mol, Val(:draw2d_bond_style)) || return
    hasfield(typeof(mol.gprops), :stereocenter) || return
    centers = mol.gprops.stereocenter
    isempty(centers) && return
    desc = mol.gprops.descriptors
    symbols = Symbol[atom_symbol(mol[i]) for i in vertices(mol)]
    # standardized (kekulized) bond orders if available
    orders = (has_descriptor(mol, Val(:bond_order)) && length(desc.bond_order) == ne(mol)
        ) ? desc.bond_order : Int[bond_order(mol[e]) for e in edges(mol)]
    for (coords, styles) in zip(desc.coords2d, desc.draw2d_bond_style)
        (length(coords) == nv(mol) && length(styles) == ne(mol)) || continue
        assign_stereocenter_wedges!(styles, mol.graph, symbols, orders, coords, centers)
    end
    return
end


"""
    stereocenter_from_sdf2d!(mol::MolGraph) -> Nothing

Set stereocenter information obtained from 2D SDFile.
"""
function stereocenter_from_sdf2d!(mol::SimpleMolGraph)
    centers, comments = stereocenter_from_sdf2d(mol)
    empty!(mol[:stereocenter])
    merge!(mol[:stereocenter], centers)
    for c in keys(mol[:stereocenter])
        safe_stereo_hydrogen!(mol, c)
    end
    if length(comments) > 0
        mol[:logs]["warning_stereocenter"] = join(comments, "; ")
    end
    return
end


"""
    stereocenter_from_smiles(
        g::SimpleGraph{T}, succ::Vector{Vector{T}}, v_stereo::Vector{Symbol}) where T -> StereocenterMap{T}

Return stereocenter information obtained from SMILES.
"""
function stereocenter_from_smiles(
        g::SimpleGraph{T}, succ::Vector{Vector{T}}, v_stereo::Vector{Symbol}) where T
    centers = StereocenterMap{T}()
    for i in vertices(g)
        degree(g, i) in (3, 4) || continue
        v_stereo[i] === :unspecified && continue
        # sort vertices by SMARTS lexicographic order
        nbrs = neighbors(g, i)
        prec = setdiff(nbrs, succ[i])
        if isempty(prec)
            # The first atom of the SMILES (or of a disconnected component) has no preceding
            # atom. The first neighbor in lexical order (the implicit hydrogen if exists) is
            # the looking-from atom (OpenSMILES).
            from, v1, v2 = succ[i][1], succ[i][2], succ[i][3]
        else
            from, v1, v2 = only(prec), succ[i][1], succ[i][2]
        end
        centers[i] = Stereocenter{T}(from, v1, v2, v_stereo[i] in (:clockwise, :clockwiseuns))
    end
    return centers
end

stereocenter_from_smiles(mol::SimpleMolGraph) = stereocenter_from_smiles(
    mol.graph, mol[:smarts_lexical_succ],
    Symbol[mol[i].stereo for i in vertices(mol)]
)

function stereocenter_from_smiles(mol::ReactiveMolGraph{T,V,E}) where {T,V<:QueryTree,E<:QueryTree}
    v_stereo = fill(:unspecified, nv(mol))
    for i in vertices(mol)
        for p in values(mol[i].vprops)
            if p.key == :stereo
                if v_stereo[i] !== :unspecified
                    error("$(i): stereochem logical operation is not supported")
                end
                v_stereo[i] = Symbol(p.value)
            end
        end
    end
    return stereocenter_from_smiles(mol.graph, mol[:smarts_lexical_succ], v_stereo)
end

"""
    stereocenter_from_smiles!(mol::MolGraph) -> Nothing

Set stereocenter information obtained from SMILES.
"""
function stereocenter_from_smiles!(mol::SimpleMolGraph)
    centers = stereocenter_from_smiles(mol)
    empty!(mol[:stereocenter])
    merge!(mol[:stereocenter], centers)
    for c in keys(mol[:stereocenter])
        safe_stereo_hydrogen!(mol, c)
    end
    return
end


"""
    stereobond_from_sdf2d(
        g::SimpleGraph{T}, sssr::Vector{Vector{T}}, e_order::Vector{Int}, e_notation::Vector{Int},
        v_coords2d::Coords2d) where T -> StereobondMap{T}

Return cis-trans diastereomerism information obtained from 2D SDFile.
"""
function stereobond_from_sdf2d(
        g::SimpleGraph{T}, sssr::Vector{Vector{T}}, e_order::Vector{Int}, e_notation::Vector{Int},
        v_coords2d::Coords2d) where T
    stereobonds = StereobondMap{T}()
    smallcycles = Edge{T}[]
    for cyc in edgemincyclebasis(g, sssr)
        length(cyc) < 8 && append!(smallcycles, cyc)
    end
    for (i, e) in enumerate(edges(g))
        e_order[i] == 2 || continue
        e_notation[i] == 3 && continue  # stereochem unspecified
        e in smallcycles && continue
        degree(g, src(e)) in (2, 3) || continue
        degree(g, dst(e)) in (2, 3) || continue
        snbrs, dnbrs = ordered_edge_neighbors(g, e)
        # Check coordinates
        d1, d2 = (v_coords2d[src(e)], v_coords2d[dst(e)])
        n1, n2 = (v_coords2d[snbrs[1]], v_coords2d[dnbrs[1]])
        cond(a, b) = (d1[1] - d2[1]) * (b - d1[2]) + (d1[2] - d2[2]) * (d1[1] - a)
        n1p = cond(n1[1], n1[2])
        n2p = cond(n2[1], n2[2])
        n1p * n2p == 0 && continue  # 180° bond angle
        is_cis = n1p * n2p > 0
        stereobonds[e] = Stereobond{T}(snbrs[1], dnbrs[1], is_cis)
    end
    return stereobonds
end

stereobond_from_sdf2d(mol::SimpleMolGraph) = stereobond_from_sdf2d(
    mol.graph,
    sssr(mol),
    Int[mol[e].order for e in edges(mol)],
    Int[mol[e].notation for e in edges(mol)],
    mol[:descriptors].coords2d[1]
)

"""
    stereobond_from_sdf2d!(mol::MolGraph) -> Nothing

Set cis-trans diastereomerism information obtained from 2D SDFile.
"""
function stereobond_from_sdf2d!(mol::SimpleMolGraph)
    empty!(mol[:stereobond])
    merge!(mol[:stereobond], stereobond_from_sdf2d(mol))
    return
end



function adj_direction!(
        v::T, issrc::Bool, enbrs::Vector{T}, dir::Vector{Symbol}, ernk::Dict{Edge{T},Int},
        comments::Vector{String}) where T
    directions = Tuple{T,Symbol}[] # (nbr, :up/:down) of bonds connected to the src(e) or dst(e)
    """
    Follows OpenSMILES specification http://opensmiles.org/opensmiles.html#chirality
        -> "up-ness" or "down-ness" of each single bond is relative to the carbon atom
    e.g.  C(\\F)=C/F -> trans
    """
    for n in enbrs
        d = dir[edge_rank(ernk, n, v)]
        d === :unspecified && continue
        pred = issrc ? n < v : v < n
        conv = Dict(:up => :down, :upuns => :down, :down => :up, :downuns => :up)
        push!(directions, (n, pred ? d : conv[d]))
    end
    # Conflicting cis/trans will be ignored (adopts OpenSMILES specs)
    if length(directions) == 2 && directions[1][2] == directions[2][2]
        v1, v2 = (directions[1][1], directions[2][1])
        @debug "ignored stereobond #$(v)-$(v1), #$(v)-$(v2) (conflicting up and down bonds)"
        push!(comments, "ignored stereobond #$(v)-$(v1), #$(v)-$(v2) (conflicting up and down bonds)")
        empty!(directions)
    end
    return directions
end

"""
    stereobond_from_smiles(
        g::SimpleGraph{T}, isdouble::BitVector, e_direction::Vector{Symbol}
    ) where T -> Tuple{StereobondMap{T},Vector{String}}

Return cis-trans diastereomerism information obtained from SMILES and log messages.
"""
function stereobond_from_smiles(
        g::SimpleGraph{T}, isdouble::BitVector, e_direction::Vector{Symbol}) where T
    stereobonds = StereobondMap{T}()
    comments = String[]
    ernk = edge_rank(g)
    for (i, e) in enumerate(edges(g))
        isdouble[i] || continue
        degree(g, src(e)) in (2, 3) || continue
        degree(g, dst(e)) in (2, 3) || continue
        snbrs, dnbrs = ordered_edge_neighbors(g, e)
        sds = adj_direction!(src(e), true, snbrs, e_direction, ernk, comments)
        isempty(sds) && continue
        dds = adj_direction!(dst(e), false, dnbrs, e_direction, ernk, comments)
        isempty(dds) && continue
        stereobonds[e] = Stereobond{T}(sds[1][1], dds[1][1], sds[1][2] !== dds[1][2])
    end
    return stereobonds, comments
end

stereobond_from_smiles(mol::SimpleMolGraph) = stereobond_from_smiles(
    mol.graph,
    BitVector(mol[e].order == 2 for e in edges(mol)),
    Symbol[mol[e].direction for e in edges(mol)]
)

function stereobond_from_smiles(mol::ReactiveMolGraph{T,V,E}) where {T,V<:QueryTree,E<:QueryTree}
    e_direction = fill(:unspecified, ne(mol))
    for (i, e) in enumerate(edges(mol))
        for p in values(mol[e].vprops)
            if p.key === :direction
                if e_direction[i] !== :unspecified
                    error("$(e): stereochem logical operation is not supported")
                end
                e_direction[i] = Symbol(p.value)
            end
        end
    end
    return stereobond_from_smiles(mol.graph, trues(ne(mol)), e_direction)
end


"""
    stereobond_from_smiles!(mol::MolGraph) -> Nothing

Set cis-trans diastereomerism information obtained from SMILES.
"""
function stereobond_from_smiles!(mol::SimpleMolGraph)
    bonds, comments = stereobond_from_smiles(mol)
    empty!(mol[:stereobond])
    merge!(mol[:stereobond], bonds)
    if length(comments) > 0
        mol[:logs]["warning_stereobond"] = join(comments, "; ")
    end
end

# TODO: axial chirality
# TODO: hypervalent chirality
