#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    addstereohydrogens, removestereohydrogens,
    setdiastereo!, optimizewedges!, setstereocenter!


"""
    addstereohydrogens(mol::SMILES) -> SMILES

Return new molecule with explicit hydrogen nodes attached to stereocenters.
"""
function addstereohydrogens(mol::SMILES)
    """
    [C@@H](C)(C)C -> [C@@]([H])(C)(C)C
    C[C@@H](C)C -> C[C@@]([H])(C)C
    """
    atoms = nodeattrs(mol)
    hcount_ = hcount(mol)
    cpmol = graphmol(mol)
    mapper = Dict{Int,Int}()
    offset = 0
    for i in 1:nodecount(mol)
        mapper[i + offset] = i
        atoms[i].stereo === :unspecified && continue
        degree(mol, i) == 3 || continue
        hcount_[i] == 1 || continue
        n = addnode!(cpmol, nodeattrtype(mol)(:H))
        addedge!(cpmol, i, n, edgeattrtype(mol)())
        offset += 1
        mapper[i + offset] = n
    end
    return remapnodes(cpmol, mapper)
end


"""
    removestereohydrogens(mol::SMILES) -> SMILES

Return new molecule without explicit hydrogen nodes attached to the stereocenters.
"""
function removestereohydrogens(mol::SMILES)
    """
    [C@@]([H])(C)(N)O -> [C@@H](C)(N)O
    ([H])[C@@](C)(N)O -> [C@@H](C)(N)O
    C[C@@]([H])(N)O -> C[C@@H](N)O
    C[C@@](N)([H])O -> C[C@H](N)O (reverse direction)
    C[C@@](N)(O)[H] -> C[C@@H](N)O
    """
    cpmol = clone(mol)
    hcount_ = hcount(cpmol)
    atoms = nodeattrs(cpmol)
    to_be_removed = Int[]
    rev = Dict(:clockwise => :anticlockwise, :anticlockwise => :clockwise)
    for c in 1:nodecount(cpmol)
        atoms[c].stereo === :unspecified && continue
        degree(cpmol, c) == 4 || continue
        hcount_[c] == 1 || continue
        for (i, adj) in enumerate(sort(collect(adjacencies(cpmol, c))))
            if atoms[adj].symbol === :H
                if i == 3
                    a = setstereo(atoms[c], rev[atoms[c].stereo])
                    setnodeattr!(cpmol, c, a)
                end
                push!(to_be_removed, adj)
            end
        end
    end
    ns = setdiff(nodeset(cpmol), to_be_removed)
    return graphmol(nodesubgraph(cpmol, ns))
end


"""
    setdiastereo!(mol::Union{SMILES,SDFile}) -> nothing

Set diastereomerism flags to `Bond.stereo` fields of double bonds.

`:unspecified`, `:cis` or `:trans` will be assigned according to configurations of bonds labeled by lower indices. (ex. SMILES `C/C=C(C)/C` have 4 explicit bonds 1:`C/C`, 2:'C=C', 3:`C(C)` and 4:`C/C`. 1st bond is prior to the other implicit hydrogen attached to 2nd atom, and 3rd bond is prior to 4th bond in index label order. 1st and 3th atom are in ``cis`` position, so `Bond.stereo` of 2nd bond will be set to `:cis`.)
"""
function setdiastereo!(mol::SMILES)
    edges = edgeattrs(mol)
    hcount_ = hcount(mol)
    for i in 1:edgecount(mol)
        edges[i].order == 2 || continue
        d = Symbol[]
        for n in mol.edges[i]
            inc = incidences(mol, n)
            if length(inc) == 3
                f, s = sort([e for e in inc if e != i])
                if edges[f].direction !== :unspecified
                    push!(d, edges[f].direction)
                elseif edges[s].direction !== :unspecified
                    push!(d, edges[s].direction === :up ? :down : :up)
                end
            elseif length(inc) == 2 && hcount_[n] == 1
                f = [e for e in inc if e != i][1]
                if edges[f].direction !== :unspecified
                    push!(d, edges[f].direction)
                end
            end
            # TODO: axial chirality
        end
        if length(d) == 2
            stereo = d[1] === d[2] ? :trans : :cis
            bond = setstereo(edges[i], stereo)
            setedgeattr!(mol, i, bond)
        end
    end
end


function setdiastereo!(mol::SDFile)
    nodes = nodeattrs(mol)
    edges = edgeattrs(mol)
    hcount_ = hcount(mol)
    coords_ = coords2d(mol)
    for i in 1:edgecount(mol)
        edges[i].order == 2 || continue
        edges[i].notation == 3 && continue  # stereochem unspecified
        # Get lower indexed edges connected to each side of the bond
        ns = Int[]
        for n in mol.edges[i]
            incs = incidences(mol, n)
            length(incs) in (2, 3) || continue
            f = sort([e for e in incs if e != i])[1]
            push!(ns, neighbors(mol, n)[f])
        end
        length(ns) == 2 || continue
        # Check coordinates
        d1, d2 = [point(coords_, n) for n in mol.edges[i]]
        n1, n2 = [point(coords_, n) for n in ns]
        cond(a, b) = (x(d1)-x(d2)) * (b-y(d1)) + (y(d1)-y(d2)) * (x(d1)-a)
        n1p = cond(x(n1), y(n1))
        n2p = cond(x(n2), y(n2))
        if n1p * n2p < 0
            stereo = :trans
        elseif n1p * n2p > 0
            stereo = :cis
        else
            stereo = :unspecified
        end
        bond = setstereo(edges[i], stereo)
        setedgeattr!(mol, i, bond)
    end
end



function angeval(u::Point2D, v::Point2D)
    # 0deg -> 1, 90deg -> 0, 180deg -> -1, 270deg-> -2, 360deg -> -3
    uv = dot(u, v) / (norm(u) * norm(v))
    return cross(u, v) >= 0 ? uv : -2 - uv
end


function anglesort(coords, center, ref, vertices)
    c = point(coords, center)
    r = point(coords, ref)
    ps = [point(coords, v) for v in vertices]
    vs = [p - c for p in ps]
    per = sortperm([angeval(r, v) for v in vs])
    return vertices[per]
end


"""
    optimizewedges!(mol::SDFile)

Optimize dashes and wedges representations. Typical stereocenters can be drawn as four bonds including only a wedge and/or a dash, so if there are too many dashes and wedges, some of them will be converted to normal single bond without changing any stereochemistry.
"""
function optimizewedges!(mol::SDFile)
    edges = edgeattrs(mol)
    coords_ = coords2d(mol)
    for i in 1:nodecount(mol)
        incs = incidences(mol, i)
        length(incs) in (3, 4) || continue
        upincs = Int[]
        downincs = Int[]
        for inc in incs
            if edges[inc].notation == 1
                push!(mol.edges[inc][1] == i ? upincs : downincs, inc)
            elseif edges[inc].notation == 6
                push!(downincs, inc)
            end
        end
        (isempty(upincs) && isempty(downincs)) && continue  # unspecified
        newbonds = Tuple{Int,Int}[]
        if length(incs) == 4
            if length(upincs) == 3
                downb = (setdiff(incs, upincs)[1], 6)
                others = [(inc, 0) for inc in upincs]
                push!(newbonds, downb)
                append!(newbonds, others)
            elseif length(downincs) == 3
                upb = (setdiff(incs, downincs)[1], 1)
                others = [(inc, 0) for inc in downincs]
                push!(newbonds, upb)
                append!(newbonds, others)
            elseif length(upincs) == 2 && length(downincs) == 1
                append!(newbonds, [(inc, 0) for inc in upincs])
            elseif length(downincs) == 2 && length(upincs) == 1
                append!(newbonds, [(inc, 0) for inc in downincs])
            else
                nbrs = neighbors(mol, i)
                badjs = [nbrs[b] for b in setdiff(incs, upincs, downincs)]
                upadjs = [nbrs[b] for b in upincs]
                downadjs = [nbrs[b] for b in downincs]
                if length(upincs) == 2 && isempty(downincs)
                    ordered = anglesort(
                        coords_, i, upadjs[1], union(badjs, [upadjs[2]]))
                    if findfirst(isequal(upadjs[2]), ordered) == 2
                        push!(newbonds,(upincs[2], 0))
                    end
                elseif length(upincs) == 2 && length(downincs) == 2
                    ordered = anglesort(
                        coords_, i, upadjs[1], union(downadjs, [upadjs[2]]))
                    if findfirst(isequal(upadjs[2]), ordered) == 2
                        push!(newbonds, (upincs[2], 0))
                    end
                elseif isempty(upincs) && length(downincs) == 2
                    ordered = anglesort(
                        coords_, i, downadjs[1], union(badjs, [downadjs[2]]))
                    if findfirst(isequal(downadjs[2]), ordered) == 2
                        push!(newbonds, (downincs[2], 0))
                    end
                elseif length(upincs) == 1 && length(downincs) == 1
                    ordered = anglesort(
                        coords_, i, badjs[1],
                        [badjs[2], upadjs[1], downadjs[1]]
                    )
                    u = point(coords_, badjs[1])
                    v = point(coords_, badjs[2])
                    d = cross(u, v)
                    if findfirst(isequal(badjs[2]), ordered) == 2
                        # if d == 0, stereochem is ambigious
                        if d > 0
                            push!(newbonds, (ordered[1], 0))
                        elseif d < 0
                            push!(newbonds, (ordered[3], 0))
                        end
                    end
                end
            end
        elseif length(incs) == 3
            if length(upinc) == 3
                append!(newbonds, [(inc, 0) for inc in upincs[1:2]])
            elseif length(downincs) == 3
                append!(newbonds, [(inc, 0) for inc in downincs[1:2]])
            elseif length(upinc) == 2
                downb = (setdiff(incs, upincs)[1], 6)
                others = [(inc, 0) for inc in upincs]
                push!(newbonds, downb)
                append!(newbonds, others)
            elseif length(downinc) == 2
                upb = (setdiff(incs, downincs)[1], 1)
                others = [(inc, 0) for inc in downincs]
                push!(newbonds, upb)
                append!(newbonds, others)
            end
        end
        for (inc, notation) in newbonds
            b = setnotation(edges[inc], notation)
            setedgeattr!(mol, inc, b)
        end
    end
end



"""
    setstereocenter!(mol::SDFile)

Set stereocenter information to Atom.stereo (`:unspecified`, `:clockwise`, `:anticlockwise` or `:atypical`). Clockwise/anticlockwise means the configuration of 2-4th nodes in index label order when we see the chiral center from the node labeled by the lowest index. If there is an implicit hydrogen, its index label will be regarded as the same as the stereocenter atom.

Note that `setstereocenter!` will optimize dashes and wedges by using `optimizewedges!` inside it. If you do not want to change the bond notation, `clone` the molecule before running `setstereocenter!`.
"""
function setstereocenter!(mol::SDFile)
    optimizewedges!(mol)
    nodes = nodeattrs(mol)
    edges = edgeattrs(mol)
    coords_ = coords2d(mol)
    for i in 1:nodecount(mol)
        adjs = Int[]
        upadjs = Int[]
        downadjs = Int[]
        for (inc, adj) in neighbors(mol, i)
            push!(adjs, adj)
            if edges[inc].notation == 1
                push!(mol.edges[inc][1] == i ? upadjs : downadjs, adj)
            elseif edges[inc].notation == 6
                push!(downadjs, adj)
            end
        end
        (isempty(upadjs) && isempty(downadjs)) && continue  # unspecified
        length(adjs) < 3 && continue # atoms attached to the chiral center
        stereo = :unspecified
        if length(adjs) > 4  # Hypervalent
            stereo = :atypical
        elseif length(adjs) == 4
            if length(upadjs) > 1 || length(downadjs) > 1
                stereo = :atypical  # Not pyramidal
            else
                # Select a pivot
                pivot = isempty(upadjs) ? downadjs[1] : upadjs[1]
                tri = setdiff(adjs, [pivot])
                quad = adjs
                rev = isempty(upadjs) ? true : false
            end
        elseif length(adjs) == 3
            if length(upadjs) + length(downadjs) > 1
                stereo = :atypical  # Axial chirality?
            else
                # Implicit hydrogen is the pivot
                pivot = i
                tri = adjs
                quad = union(adjs, [i])
                rev = isempty(upadjs) ? true : false
            end
        end
        if stereo === :unspecified
            cw = isclockwise(cartesian2d(coords_, sort(tri)))
            if rev
                cw = !cw
            end
            # Arrange the configuration so that the lowest index node is the pivot.
            pividx = findfirst(isequal(pivot), sort(quad))
            if pividx in (2, 4)
                cw = !cw
            end
            stereo = cw ? :clockwise : :anticlockwise
        end
        a = setstereo(nodes[i], stereo)
        setnodeattr!(mol, i, a)
    end
end
