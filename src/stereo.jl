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


function isclockwise(stereo::Tuple{T,T,T,Bool}, f::T, s::T, t::T) where T
    fp = something(findfirst(==(f), stereo[1:3]), 4)
    sp = something(findfirst(==(s), stereo[1:3]), 4)
    tp = something(findfirst(==(t), stereo[1:3]), 4)
    return stereo[4] === STEREOCENTER_STATE[(fp, sp, tp)]
end


function remap!(
        ::Val{:stereocenter}, gprop::SimpleMolProperty{T}, vmap::Dict{T,T}
        ) where T <: Integer
    newmap = Dict{T,Tuple{T,T,T,Bool}}()
    for (k, v) in gprop.stereocenter
        isempty(setdiff([k, v[1:3]...], keys(vmap))) || continue
        newmap[vmap[k]] = (vmap[v[1]], vmap[v[2]], vmap[v[3]], v[4])
    end
    gprop.stereocenter = newmap
    return
end

function reconstruct(::Val{:stereocenter}, T::Type{<:SimpleMolProperty}, data)
    U = eltype(T)
    return Dict{U,Tuple{U,U,U,Bool}}(parse(U, i) => tuple(val...) for (i, val) in data)
end

function to_dict(::Val{:stereocenter}, ::Val{:default}, gprop::AbstractProperty)
    return Dict{String,Any}(string(i) => collect(val) for (i, val) in gprop.stereocenter)
end


function remap!(
        ::Val{:stereobond}, gprop::SimpleMolProperty{T}, vmap::Dict{T,T}
        ) where T <: Integer
    newmap = Dict{Edge{T},Tuple{T,T,Bool}}()
    for (k, v) in gprop.stereobond
        isempty(setdiff([src(k), dst(k), v[1:2]...], keys(vmap))) || continue
        newmap[u_edge(T, vmap[src(k)], vmap[dst(k)])] = (vmap[v[1]], vmap[v[2]], v[3])
    end
    gprop.stereobond = newmap
    return
end

function reconstruct(::Val{:stereobond}, T::Type{<:SimpleMolProperty}, data)
    U = eltype(T)
    return Dict{Edge{U},Tuple{U,U,Bool}}(
        Edge{U}(s, d) => tuple(val...) for (s, d, val) in data)
end

function to_dict(::Val{:stereobond}, ::Val{:default}, gprop::AbstractProperty)
    return [[src(e), dst(e), collect(val)] for (e, val) in gprop.stereobond]
end


"""
    set_stereocenter!(mol::ReactiveMolGraph, center, looking_from, v1, v2, is_clockwise) -> Nothing

Set stereocenter information to graph properties.
"""
function set_stereocenter!(
        mol::ReactiveMolGraph{T,V,E}, center,
        looking_from, v1, v2, is_clockwise) where {T,V,E}
    mol.gprops.stereocenter[center] = (looking_from, v1, v2, is_clockwise)
end


"""
    set_stereocenter!(mol::ReactiveMolGraph, bond, v1, v2, is_cis) -> Nothing

Set stereocenter information to graph properties.
"""
function set_stereobond!(
        mol::ReactiveMolGraph{T,V,E}, bond, v1, v2, is_cis) where {T,V,E}
    mol.gprops.stereobond[bond] = (v1, v2, is_cis)
end


function stereo_hydrogen(mol::SimpleMolGraph, v::Integer)
    nbrs = neighbors(mol, v)
    hpos = findfirst(x -> atom_symbol(props(mol, x)) === :H, nbrs)
    isnothing(hpos) && return  # 4° center or already removed
    return nbrs[hpos]
end


"""
    safe_stereo_hydrogen!(mol::SimpleMolGraph, v::Integer) -> Bool

Rearrange stereocenter properties to safely remove stereo hydrogen nodes
and return the hydrogen node index.

This function is called inside `rem_vertex!` and `rem_vertices!` functions
to safely remove hydrogen nodes while preserving stereocenter information.
"""
function safe_stereo_hydrogen!(mol::ReactiveMolGraph{T,V,E}, center::T) where {T,V,E}
    h = stereo_hydrogen(mol, center)
    isnothing(h) && return # 4° center or already removed
    stereo = get_prop(mol, :stereocenter)[center]
    spos = findfirst(==(h), stereo[1:3])
    isnothing(spos) && return h  # hydrogen at the lowest priority can be removed safely
    nonh = setdiff(ordered_neighbors(mol, center), h)
    set_stereocenter!(mol, center, nonh..., isclockwise(stereo, nonh...))
    return h
end



function angeval(u::Point2d, v::Point2d)
    # 0deg -> 1, 90deg -> 0, 180deg -> -1, 270deg-> -2, 360deg -> -3
    uv = dot(u, v) / (norm(u) * norm(v))
    return cross(u, v) >= 0 ? uv : -2 - uv
end


function anglesort(coords, center, ref, vertices)
    # return vertices order by clockwise direction
    c = coords[center]
    r = coords[ref]
    ps = [coords[v] for v in vertices]
    vs = [p - c for p in ps]
    return sortperm([angeval(r - c, v) for v in vs])
end


"""
    stereocenter_from_sdf2d(mol::MolGraph{T,V,E}) where {T,V,E} -> Dict{T,Tuple{T,T,T,Bool}}

Return stereocenter information obtained from 2D SDFile.
"""
function stereocenter_from_sdf2d(
        g::SimpleGraph{T}, v_symbol, e_order, e_notation, e_isordered, v_coords2d) where T
    centers = Dict{T,Tuple{T,T,T,Bool}}()
    edgerank = Dict(e => i for (i, e) in enumerate(edges(g)))
    comments = String[]
    for i in vertices(g)
        degree(g, i) in (3, 4) || continue
        nbrs = ordered_neighbors(g, i)
        drs = Symbol[]  # lookingFrom, atom1, atom2, (atom3)
        for nbr in nbrs
            e_order[edgerank[u_edge(g, i, nbr)]] == 1 || break
            if e_isordered[edgerank[u_edge(g, i, nbr)]] == (i < nbr)  # only outgoing wedges are considered
                if e_notation[edgerank[u_edge(g, i, nbr)]] == 1
                    push!(drs, :up)
                    continue
                elseif e_notation[edgerank[u_edge(g, i, nbr)]] == 6
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
        sortorder = [1, map(x -> x + 1, anglesort(v_coords2d, i, nbrs[1], nbrs[2:end]))...]
        ons = nbrs[sortorder]  # node indices of neighbors in angular order
        ods = drs[sortorder]  # direction symbols of neighbors in angular order
        if length(nbrs) == 3
            if upcnt == 1 && dwcnt == 1 || upcnt + dwcnt == 3
                @debug "$(i): a reference plane cannot be defined $(ons) $(ods)"
                push!(comments, "$(i): a reference plane cannot be defined")
            elseif upcnt == 2 || dwcnt == 2
                centers[i] = (ons[1], ons[2], ons[3], upcnt == 0)
            else
                centers[i] = (ons[1], ons[2], ons[3], dwcnt == 0)
            end
            continue
        end
        if upcnt == 4 || dwcnt == 4
            @debug "$(i): a reference plane cannot be defined $(ons) $(ods)"
            push!(comments, "$(i): a reference plane cannot be defined")
        elseif upcnt == 3
            centers[i] = (ons[1], ons[2], ons[3], ods[1] === :up && ods[3] === :up)
        elseif dwcnt == 3
            centers[i] = (ons[1], ons[2], ons[3], ods[2] === :down && ods[4] === :down)
        elseif upcnt == 2
            if (ods[1] === :up && ods[3] === :up) || (ods[2] === :up && ods[4] === :up)
                centers[i] = (ons[1], ons[2], ons[3], ods[1] === :up)
            else
                @debug "$(i): adjacent wedges lie on the same side of the reference plane $(ons) $(ods)"
                push!(comments, "$(i): adjacent wedges lie on the same side of the reference plane")
            end
        elseif dwcnt == 2
            if (ods[1] === :down && ods[3] === :down) || (ods[2] === :down && ods[4] === :down)
                centers[i] = (ons[1], ons[2], ons[3], ods[1] !== :down)
            else
                @debug "$(i): adjacent wedges lie on the same side of the reference plane $(ons) $(ods)"
                push!(comments, "$(i): adjacent wedges lie on the same side of the reference plane")
            end
        elseif upcnt == 1 && dwcnt == 1
            if ((ods[1] === :unspecified && ods[3] === :unspecified)
                    || (ods[2] === :unspecified && ods[4] === :unspecified))
                @debug "$(i): non-adjacent wedges lie on the opposite side of the reference plane $(ons) $(ods)"
                push!(comments, "$(i): non-adjacent wedges lie on the opposite side of the reference plane")
            else
                centers[i] = (ons[1], ons[2], ons[3], ods[1] === :up || ods[3] === :up)
            end
        elseif upcnt == 1
            centers[i] = (ons[1], ons[2], ons[3], ods[1] === :up || ods[3] === :up)
        else  # dwcnt == 1
            centers[i] = (ons[1], ons[2], ons[3], ods[2] === :down || ods[4] === :down)
        end
    end
    return centers, comments
end

stereocenter_from_sdf2d(mol::ReactiveMolGraph) = stereocenter_from_sdf2d(
    mol.graph,
    [get_prop(mol, i, :symbol) for i in vertices(mol)],
    [get_prop(mol, e, :order) for e in edges(mol)],
    [get_prop(mol, e, :notation) for e in edges(mol)],
    [get_prop(mol, e, :isordered) for e in edges(mol)],
    get_descriptor(mol, :coords2d)[1]
)

"""
    stereocenter_from_sdf2d!(mol::MolGraph) -> Nothing

Set stereocenter information obtained from 2D SDFile.
"""
function stereocenter_from_sdf2d!(mol::ReactiveMolGraph)
    centers, comments = stereocenter_from_sdf2d(mol)
    set_prop!(mol, :stereocenter, centers)
    for c in keys(get_prop(mol, :stereocenter))
        safe_stereo_hydrogen!(mol, c)
    end
    if length(comments) > 0
        get_prop(mol, :logs)["stereocenter_ignored"] = join(comments, "; ")
    end
    return
end


"""
    stereocenter_from_smiles(g::SimpleGraph{T}, v_stereo) where T -> Dict{T,Tuple{T,T,T,Bool}}
Return stereocenter information obtained from SMILES.
"""
function stereocenter_from_smiles(g::SimpleGraph{T}, succ, v_stereo) where T
    centers = Dict{T,Tuple{T,T,T,Bool}}()
    for i in vertices(g)
        degree(g, i) in (3, 4) || continue
        v_stereo[i] === :unspecified && continue
        # sort vertices by SMARTS lexicographic order
        nbrs = neighbors(g, i)
        first = only(setdiff(nbrs, succ[i]))
        centers[i] = (first, succ[i][1], succ[i][2], v_stereo[i] === :clockwise)
    end
    return centers
end

stereocenter_from_smiles(mol::SimpleMolGraph) = stereocenter_from_smiles(
    mol.graph, get_prop(mol, :smarts_lexical_succ),
    [get_prop(mol, i, :stereo) for i in vertices(mol)]
)

"""
    stereocenter_from_smiles!(mol::MolGraph) -> Nothing

Set stereocenter information obtained from SMILES.
"""
function stereocenter_from_smiles!(mol::ReactiveMolGraph)
    centers = stereocenter_from_smiles(mol)
    set_prop!(mol, :stereocenter, centers)
    for c in keys(get_prop(mol, :stereocenter))
        safe_stereo_hydrogen!(mol, c)
    end
    return
end


"""
    stereobond_from_sdf2d(g::SimpleGraph{T}, e_order, e_notation, v_coords2d) where T -> Dict{Edge{T},Tuple{T,T,Bool}}

Return cis-trans diastereomerism information obtained from 2D SDFile.
"""
function stereobond_from_sdf2d(g::SimpleGraph{T}, e_order, e_notation, v_coords2d) where T
    stereobonds = Dict{Edge{T},Tuple{T,T,Bool}}()
    smallcycles = Edge{T}[]
    for cyc in edgemincyclebasis(g)
        length(cyc) < 8 && push!(smallcycles, cyc...)
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
        stereobonds[e] = (snbrs[1], dnbrs[1], is_cis)
    end
    return stereobonds
end

stereobond_from_sdf2d(mol::ReactiveMolGraph) = stereobond_from_sdf2d(
    mol.graph,
    [get_prop(mol, e, :order) for e in edges(mol)],
    [get_prop(mol, e, :notation) for e in edges(mol)],
    mol.gprops.descriptors.coords2d[1]
)

"""
    stereobond_from_sdf2d!(mol::MolGraph) -> Nothing

Set cis-trans diastereomerism information obtained from 2D SDFile.
"""
stereobond_from_sdf2d!(mol::ReactiveMolGraph) = setproperty!(
    mol.gprops, :stereobond, stereobond_from_sdf2d(mol)
)


"""
    stereobond_from_smiles(g::SimpleGraph{T}, e_order, e_direction) where T -> Dict{Edge{T},Tuple{T,T,Bool}}

Return cis-trans diastereomerism information obtained from SMILES.
"""
function stereobond_from_smiles(g::SimpleGraph{T}, e_order, e_direction) where T
    stereobonds = Dict{Edge{T},Tuple{T,T,Bool}}()
    comments = String[]
    edgerank = Dict(e => i for (i, e) in enumerate(edges(g)))
    for (i, e) in enumerate(edges(g))
        e_order[i] == 2 || continue
        degree(g, src(e)) in (2, 3) || continue
        degree(g, dst(e)) in (2, 3) || continue
        snbrs, dnbrs = ordered_edge_neighbors(g, e)
        sds = []  # (nbr, :up/:down) of bonds connected to src(e)
        dds = []  # (nbr, :up/:down) of bonds connected to dst(e)
        """
        Follows OpenSMILES specification http://opensmiles.org/opensmiles.html#chirality
          -> "up-ness" or "down-ness" of each single bond is relative to the carbon atom
        e.g.  C(\\F)=C/F -> trans
        """
        for sn in snbrs
            sd = e_direction[edgerank[u_edge(g, sn, src(e))]]
            sd === :unspecified && continue
            sd = sn < src(e) ? sd : (sd === :up ? :down : :up)
            push!(sds, (sn, sd))
        end
        for dn in dnbrs
            dd = e_direction[edgerank[u_edge(g, dn, dst(e))]]
            dd === :unspecified && continue
            dd = dst(e) < dn ? dd : (dd === :up ? :down : :up)
            push!(dds, (dn, dd))
        end
        (isempty(sds) || isempty(dds)) && continue  # no :up or :down bonds
        # Conflicting cis/trans will be ignored (adopts OpenSMILES specs)
        if (length(sds) == 2 && sds[1][2] == sds[2][2]) || (
                length(dds) == 2 && dds[1][2] == dds[2][2])
            @debug "$(e): conflicting up and down bonds $(sds) $(dds)"
            push!(comments, "$(e): conflicting up and down bonds")
            continue
        end
        stereobonds[e] = (sds[1][1], dds[1][1], sds[1][2] !== dds[1][2])
    end
    return stereobonds, comments
end

stereobond_from_smiles(mol::ReactiveMolGraph) = stereobond_from_smiles(
    mol.graph,
    [get_prop(mol, e, :order) for e in edges(mol)],
    [get_prop(mol, e, :direction) for e in edges(mol)]
)

"""
    stereobond_from_smiles!(mol::MolGraph) -> Nothing

Set cis-trans diastereomerism information obtained from SMILES.
"""
function stereobond_from_smiles!(mol::ReactiveMolGraph)
    bonds, comments = stereobond_from_smiles(mol)
    setproperty!(mol.gprops, :stereobond, bonds)
    if length(comments) > 0
        mol.gprops.logs["stereobond_ignored"] = join(comments, "; ")
    end
end

# TODO: axial chirality
# TODO: hypervalent chirality
