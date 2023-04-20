#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Stereocenter, Stereobond,
    set_stereocenter!, set_stereobond!, remove_stereo_hydrogen!,
    stereocenter_from_smiles!, stereocenter_from_sdf2d!,
    stereobond_from_smiles!, stereobond_from_sdf2d!

struct Stereocenter{T} <: AbstractDict{T,Tuple{T,T,T,Bool}}
    mapping::Dict{T,Tuple{T,T,T,Bool}}
end

Stereocenter{T}(data::Vector=[]) where T = Stereocenter{T}(
    Dict{T,Tuple{T,T,T,Bool}}(i => tuple(val...) for (i, val) in data))

Base.iterate(stereo::Stereocenter) = iterate(stereo.mapping)
Base.iterate(stereo::Stereocenter, i) = iterate(stereo.mapping, i)
Base.length(stereo::Stereocenter) = length(stereo.mapping)
Base.get(stereo::Stereocenter, k, v) = get(stereo.mapping, k, v)
Base.setindex!(stereo::Stereocenter, v, k) = setindex!(stereo.mapping, v, k)
to_dict(stereo::Stereocenter) = [[i, val] for (i, val) in stereo.mapping]


struct Stereobond{T} <: AbstractDict{Edge{T},Tuple{T,T,Bool}}
    mapping::Dict{Edge{T},Tuple{T,T,Bool}}
end

Stereobond{T}(data::Vector=[]) where T = Stereobond{T}(
    Dict{Edge{T},Tuple{T,T,Bool}}(Edge{T}(s, d) => tuple(val...) for (s, d, val) in data))

Base.iterate(stereo::Stereobond) = iterate(stereo.mapping)
Base.iterate(stereo::Stereobond, i) = iterate(stereo.mapping, i)
Base.length(stereo::Stereobond) = length(stereo.mapping)
Base.get(stereo::Stereobond, k, v) = get(stereo.mapping, k, v)
Base.setindex!(stereo::Stereobond, v, k) = setindex!(stereo.mapping, v, k)
to_dict(stereo::Stereobond) = [[src(e), dst(e), val] for (e, val) in stereo.mapping]


"""
    set_stereocenter!(mol::SimpleMolGraph, center, looking_from, v1, v2, is_clockwise) -> Nothing

Set stereocenter information to graph properties.
"""
function set_stereocenter!(
        mol::SimpleMolGraph{T,V,E}, center,
        looking_from, v1, v2, is_clockwise) where {T,V,E}
    if !haskey(mol.gprops, :stereocenter)
        mol.gprops[:stereocenter] = Stereocenter{T}()
    end
    mol.gprops[:stereocenter][center] = (looking_from, v1, v2, is_clockwise)
end


"""
    set_stereocenter!(mol::SimpleMolGraph, bond, v1, v2, is_cis) -> Nothing

Set stereocenter information to graph properties.
"""
function set_stereobond!(
        mol::SimpleMolGraph{T,V,E}, bond, v1, v2, is_cis) where {T,V,E}
    if !haskey(mol.gprops, :stereobond)
        mol.gprops[:stereobond] = Stereobond{T}()
    end
    mol.gprops[:stereobond][bond] = (v1, v2, is_cis)
end


"""
    remove_stereo_hydrogen!(mol::SimpleMolGraph, v::Integer) -> Bool

Safely remove explicit hydrogens connected to stereocenter node `v`.
"""
function remove_stereo_hydrogen!(mol::SimpleMolGraph, v::Integer)
    """
    [C@@]([H])(C)(N)O -> C, N, O, (H), @
    ([H])[C@@](C)(N)O -> C, N, O, (H), @
    C[C@@]([H])(N)O -> C, N, O, (H), @@
    C[C@@](N)([H])O -> C, N, O, (H), @
    C[C@@](N)(O)[H] -> C, N, O, (H), @@
    """
    nbrs = neighbors(mol, v)
    hpos = findfirst(x -> get_prop(mol, x, :symbol) === :H, nbrs)
    hpos === nothing && return false  # no removable hydrogens
    h = nbrs[hpos]
    rem_vertex!(mol, h) || return false  # failed to remove hydrogen node
    stereo = get_prop(mol, :stereocenter)[v]
    vs = collect(stereo[1:3])
    spos = findfirst(x -> x == h, vs)
    spos === nothing && return true  # hydrogen at the lowest priority can be removed safely
    is_rev = spos in [1, 3]
    rest = setdiff(nbrs, vs)  # the lowest node index
    resti = isempty(rest) ? h : rest[1]  # rest is empty if the lowest node is the end node.
    popat!(vs, spos)
    push!(vs, resti)
    set_stereocenter!(mol, v, vs[1], vs[2], vs[3], xor(stereo[4], is_rev))
    return true
end



function angeval(u::Point2D, v::Point2D)
    # 0deg -> 1, 90deg -> 0, 180deg -> -1, 270deg-> -2, 360deg -> -3
    uv = dot(u, v) / (norm(u) * norm(v))
    return cross2d(u, v) >= 0 ? uv : -2 - uv
end


function anglesort(coords, center, ref, vertices)
    # return vertices order by clockwise direction
    c = Point2D(coords, center)
    r = Point2D(coords, ref)
    ps = [Point2D(coords, v) for v in vertices]
    vs = [p - c for p in ps]
    return sortperm([angeval(r, v) for v in vs])
end


"""
    stereocenter_from_sdf2d(mol::MolGraph{T,V,E}) where {T,V,E} -> Stereocenter{T}

Return stereocenter information obtained from 2D SDFile.
"""
function stereocenter_from_sdf2d(mol::MolGraph{T,V,E}) where {T,V,E}
    centers = Stereocenter{T}()
    crds = coords2d(mol)
    for i in vertices(mol)
        degree(mol.graph, i) in (3, 4) || continue
        nbrs = ordered_neighbors(mol, i)
        drs = Symbol[]  # lookingFrom, atom1, atom2, (atom3)
        for nbr in nbrs
            if get_prop(mol, i, nbr, :isordered) == i < nbr  # only outgoing wedges are considered
                if get_prop(mol, i, nbr, :notation) == 1
                    push!(drs, :up)
                    continue
                elseif get_prop(mol, i, nbr, :notation) == 6
                    push!(drs, :down)
                    continue
                end
            end
            push!(drs, :unspecified)
        end
        upcnt = count(x -> x === :up, drs)
        dwcnt = count(x -> x === :down, drs)
        (upcnt == 0 && dwcnt == 0) && continue  # unspecified
        sortorder = [1, map(x -> x + 1, anglesort(crds, i, nbrs[1], nbrs[2:end]))...]
        ons = nbrs[sortorder]
        ods = drs[sortorder]
        if length(nbrs) == 3
            if upcnt + dwcnt == 1
                # if there is implicit hydrogen, an up-wedge -> clockwise, a down-edge -> ccw
                centers[i] = (ons[1], ons[2], ons[3], dwcnt == 0)
            else
                # @debug "Ambiguous stereochemistry (maybe need explicit hydrogen)"
                # need warning?
                continue
            end
        elseif (ods[1] !== :unspecified && ods[3] !== :unspecified && ods[1] !== ods[3] ||
                ods[2] !== :unspecified && ods[4] !== :unspecified && ods[2] !== ods[4])
            # @debug "Ambiguous stereochemistry (opposite direction wedges at opposite side)"
            # need warning?
            continue
        elseif (ods[1] !== :unspecified && ods[1] === ods[2] ||
                ods[2] !== :unspecified && ods[2] === ods[3] ||
                ods[3] !== :unspecified && ods[3] === ods[4] ||
                ods[4] !== :unspecified && ods[4] === ods[1])
            # @debug "Ambiguous stereochemistry (same direction wedges at same side)"
            # need warning?
            continue
        elseif upcnt != 0
            centers[i] = (ons[1], ons[2], ons[3], ods[1] === :up || ods[3] === :up)
        else  # down wedges only
            centers[i] = (ons[1], ons[2], ons[3], ods[2] === :down || ods[4] === :down)
        end
    end
    return centers
end

"""
    stereocenter_from_sdf2d!(mol::MolGraph) -> Nothing

Set stereocenter information obtained from 2D SDFile.
"""
stereocenter_from_sdf2d!(mol::MolGraph
    ) = begin mol.gprops[:stereocenter] = stereocenter_from_sdf2d(mol) end


"""
    stereocenter_from_smiles(mol::MolGraph{T,V,E}) where {T,V,E} -> Stereocenter{T}
Return stereocenter information obtained from SMILES.
"""
function stereocenter_from_smiles(mol::MolGraph{T,V,E}) where {T,V,E}
    centers = Stereocenter{T}()
    for i in vertices(mol)
        degree(mol.graph, i) in (3, 4) || continue
        direction = get_prop(mol, i, :stereo)
        direction === :unspecified && continue
        nbrs = ordered_neighbors(mol, i)
        centers[i] = (nbrs[1], nbrs[2], nbrs[3], direction === :clockwise)
    end
    return centers
end

"""
    stereocenter_from_smiles!(mol::MolGraph) -> Nothing

Set stereocenter information obtained from SMILES.
"""
stereocenter_from_smiles!(mol::MolGraph
    ) = begin mol.gprops[:stereocenter] = stereocenter_from_smiles(mol) end


"""
    stereobond_from_sdf2d(mol::MolGraph{T,V,E}) where {T,V,E} -> Stereobond{T}

Return cis-trans diastereomerism information obtained from 2D SDFile.
"""
function stereobond_from_sdf2d(mol::MolGraph{T,V,E}) where {T,V,E}
    stereobonds = Stereobond{T}()
    crds = coords2d(mol)
    for e in edges(mol)
        get_prop(mol, e, :order) == 2 || continue
        get_prop(mol, e, :notation) == 3 && continue  # stereochem unspecified
        degree(mol.graph, src(e)) in (2, 3) || continue
        degree(mol.graph, dst(e)) in (2, 3) || continue
        snbrs, dnbrs = ordered_edge_neighbors(mol, e)
        # Check coordinates
        d1, d2 = (Point2D(crds, src(e)), Point2D(crds, dst(e)))
        n1, n2 = (Point2D(crds, snbrs[1]), Point2D(crds, dnbrs[1]))
        cond(a, b) = (d1.x - d2.x) * (b - d1.y) + (d1.y - d2.y) * (d1.x - a)
        n1p = cond(n1.x, n1.y)
        n2p = cond(n2.x, n2.y)
        n1p * n2p == 0 && continue  # 180° bond angle
        is_cis = n1p * n2p > 0
        stereobonds[e] = (snbrs[1], dnbrs[1], is_cis)
    end
    return stereobonds
end

"""
    stereobond_from_sdf2d!(mol::MolGraph) -> Nothing

Set cis-trans diastereomerism information obtained from 2D SDFile.
"""
stereobond_from_sdf2d!(mol::MolGraph
    ) = begin mol.gprops[:stereobond] = stereobond_from_sdf2d(mol) end


"""
    stereobond_from_smiles(mol::MolGraph{T,V,E}) where {T,V,E} -> Stereobond{T}

Return cis-trans diastereomerism information obtained from SMILES.
"""
function stereobond_from_smiles(mol::MolGraph{T,V,E}) where {T,V,E}
    stereobonds = Stereobond{T}()
    for e in edges(mol)
        get_prop(mol, e, :order) == 2 || continue
        degree(mol.graph, src(e)) in (2, 3) || continue
        degree(mol.graph, dst(e)) in (2, 3) || continue
        snbrs, dnbrs = ordered_edge_neighbors(mol, e)
        sds = []
        dds = []
        for sn in snbrs
            sd = get_prop(mol, sn, src(e), :direction)
            sd === :unspecified && continue
            push!(sds, (sn, sd))
        end
        for dn in dnbrs
            dd = get_prop(mol, dn, dst(e), :direction)
            dd === :unspecified && continue
            push!(dds, (dn, dd))
        end
        isempty(sds) && continue
        isempty(dds) && continue
        length(sds) == 2 && sds[1][2] == sds[2][2] && error("Invalid diastereomer representation")
        length(dds) == 2 && dds[1][2] == dds[2][2] && error("Invalid diastereomer representation")
        """
        follows OpenSMILES specification http://opensmiles.org/opensmiles.html#chirality
          -> "up-ness" or "down-ness" of each single bond is relative to the carbon atom
        e.g.  C(\\F)=C/F -> trans
        """
        is_cis = (sds[1][2] !== dds[1][2]) == (sds[1][1] < src(e))
        stereobonds[e] = (sds[1][1], dds[1][1], is_cis)
    end
    return stereobonds
end

"""
    stereobond_from_smiles!(mol::MolGraph) -> Nothing

Set cis-trans diastereomerism information obtained from SMILES.
"""
stereobond_from_smiles!(mol::MolGraph
    ) = begin mol.gprops[:stereobond] = stereobond_from_smiles(mol) end

# TODO: axial chirality
# TODO: hypervalent chirality
