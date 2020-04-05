#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    addstereohydrogens, removestereohydrogens,
    setdiastereo!, optimizewedges!, setstereocenter!


"""
    addstereohydrogens(mol::GraphMol) -> GraphMol

Return new molecule with explicit hydrogen nodes attached to stereocenters.
"""
function addstereohydrogens(mol::GraphMol)
    """
    [C@@H](C)(C)C -> [C@@]([H])(C)(C)C
    C[C@@H](C)C -> C[C@@]([H])(C)C
    """
    atoms = nodeattrs(mol)
    hcount_ = hcount(mol)
    newmol = clone(mol)
    mapper = Dict{Int,Int}() 
    offset = 0
    for i in 1:nodecount(mol)
        mapper[i] = i + offset
        atoms[i].stereo === :unspecified && continue
        degree(mol, i) == 3 || continue
        hcount_[i] == 1 || continue
        n = addnode!(newmol, nodeattrtype(mol)(:H))
        addedge!(newmol, i, n, edgeattrtype(mol)())
        offset += 1
        mapper[n] = i + offset
    end
    return remapnodes(newmol, mapper)
end


"""
    removestereohydrogens(mol::GraphMol) -> GraphMol

Return new molecule without explicit hydrogen nodes attached to stereocenters.
"""
function removestereohydrogens(mol::GraphMol)
    """
    [C@@]([H])(C)(N)O -> [C@@H](C)(N)O
    ([H])[C@@](C)(N)O -> [C@@H](C)(N)O
    C[C@@]([H])(N)O -> C[C@@H](N)O
    C[C@@](N)([H])O -> C[C@H](N)O (reverse direction)
    C[C@@](N)(O)[H] -> C[C@@H](N)O
    """
    atoms = nodeattrs(mol)
    nodes_to_rm = Int[]
    nodes_to_rev = Int[]
    for c in 1:nodecount(mol)
        atoms[c].stereo === :unspecified && continue
        adjs = adjacencies(mol, c)
        length(adjs) == 4 || continue
        hs = Int[]
        rev = false
        for (i, adj) in enumerate(sort(collect(adjs)))
            if atoms[adj].symbol === :H
                if i == 3
                    rev = true
                end
                push!(hs, adj)
            end
        end
        length(hs) == 1 || continue
        append!(nodes_to_rm, hs)
        rev && push!(nodes_to_rev, c)
    end
    newmol = clone(mol)
    for c in nodes_to_rev
        a = setstereo(atoms[c],
            atoms[c].stereo === :clockwise ? :anticlockwise : :clockwise)
        setnodeattr!(newmol, c, a)
    end
    ns = setdiff(nodeset(mol), nodes_to_rm)
    return graphmol(nodesubgraph(newmol, ns))
end


"""
    setdiastereo!(mol::Union{SMILES,SDFile}) -> nothing

Set diastereomerism flags to `Bond.stereo` fields of double bonds.

`:unspecified`, `:cis` or `:trans` will be assigned according to configurations of bonds labeled by lower indices. (ex. SMILES `C/C=C(C)/C` have 4 explicit bonds 1:`C/C`, 2:`C=C`, 3:`C(C)` and 4:`C/C`. 1st bond is prior to the other implicit hydrogen attached to 2nd atom, and 3rd bond is prior to 4th bond in index label order. 1st and 3th atom are in ``cis`` position, so `Bond.stereo` of 2nd bond will be set to `:cis`.)
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
        d1, d2 = [Point2D(coords_, n) for n in mol.edges[i]]
        n1, n2 = [Point2D(coords_, n) for n in ns]
        cond(a, b) = (d1.x - d2.x) * (b - d1.y) + (d1.y - d2.y) * (d1.x - a)
        n1p = cond(n1.x, n1.y)
        n2p = cond(n2.x, n2.y)
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
    return cross2d(u, v) >= 0 ? uv : -2 - uv
end


function anglesort(coords, center, ref, vertices)
    c = Point2D(coords, center)
    r = Point2D(coords, ref)
    ps = [Point2D(coords, v) for v in vertices]
    vs = [p - c for p in ps]
    per = sortperm([angeval(r, v) for v in vs])
    return vertices[per]
end


"""
    optimizewedges!(mol::SDFile)

Optimize dashes and wedges representations. Typical stereocenters can be drawn as four bonds including only a wedge and/or a dash, so if there are too many dashes and wedges, some of them will be converted to normal single bond without changing any stereochemistry.
"""
function optimizewedges!(mol::SDFile)
    # TODO: no longer used
    edges = edgeattrs(mol)
    coords_ = coords2d(mol)
    for i in 1:nodecount(mol)
        incs = sort(collect(incidences(mol, i)))
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
                downb = (pop!(setdiff(incs, upincs)), 6)
                others = [(inc, 0) for inc in upincs]
                push!(newbonds, downb)
                append!(newbonds, others)
            elseif length(downincs) == 3
                upb = (pop!(setdiff(incs, downincs)), 1)
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
                    u = Point2D(coords_, badjs[1])
                    v = Point2D(coords_, badjs[2])
                    d = cross2d(u, v)
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
            if length(upincs) == 3
                append!(newbonds, [(inc, 0) for inc in upincs[1:2]])
            elseif length(downincs) == 3
                append!(newbonds, [(inc, 0) for inc in downincs[1:2]])
            elseif length(upincs) == 2
                downb = (pop!(setdiff(incs, upincs)), 6)
                others = [(inc, 0) for inc in upincs]
                push!(newbonds, downb)
                append!(newbonds, others)
            elseif length(downincs) == 2
                upb = (pop!(setdiff(incs, downincs)), 1)
                others = [(inc, 0) for inc in downincs]
                push!(newbonds, upb)
                append!(newbonds, others)
            end
        end
        for (inc, notation) in newbonds
            if mol.edges[inc][2] == i
                mol.edges[inc] = (i, mol.edges[inc][1])  # Reverse edge
            end
            b = setnotation(edges[inc], notation)
            setedgeattr!(mol, inc, b)
        end
    end
end



"""
    setstereocenter!(mol::SDFile)

Set stereocenter information to Atom.stereo (`:unspecified`, `:clockwise`, `:anticlockwise` or `:atypical`). Clockwise/anticlockwise means the configuration of 2-4th nodes in index label order when we see the chiral center from the node labeled by the lowest index. If there is an implicit hydrogen, its index label will be regarded as the same as the stereocenter atom.
"""
function setstereocenter!(mol::SDFile)
    nodes = nodeattrs(mol)
    edges = edgeattrs(mol)
    coords_ = coords2d(mol)
    for i in 1:nodecount(mol)
        nbrs = neighbors(mol, i)
        length(nbrs) < 3 && continue # atoms attached to the chiral center
        adjs = Int[]
        upadjs = Int[]
        downadjs = Int[]
        for (inc, adj) in nbrs
            push!(adjs, adj)
            if edges[inc].notation == 1 && mol.edges[inc][1] == i
                push!(upadjs, adj)
            elseif edges[inc].notation == 6 && mol.edges[inc][1] == i
                push!(downadjs, adj)
            end
        end
        (isempty(upadjs) && isempty(downadjs)) && continue  # unspecified
        if (length(adjs) > 4 ||
            length(adjs) - length(upadjs) - length(downadjs) < 2)
            # Hypervalent, not sp3, axial chirality or wrong wedges
            a = setstereo(nodes[i], :atypical)
            setnodeattr!(mol, i, a)
            continue
        end
        if length(adjs) == 4
            # Select a pivot
            pivot = isempty(upadjs) ? downadjs[1] : upadjs[1]
            tri = setdiff(adjs, [pivot])
            quad = adjs
            rev = isempty(upadjs) ? true : false
            # Check wedges
            if length(downadjs) == 2 || length(upadjs) == 2
                ordered = anglesort(coords_, i, pivot, tri)
                sec = isempty(upadjs) ? downadjs[2] : upadjs[2]
                if findfirst(isequal(sec), ordered) != 2
                    a = setstereo(nodes[i], :atypical)
                    setnodeattr!(mol, i, a)
                    continue
                end
            elseif length(downadjs) == 1 && length(upadjs) == 1
                ordered = anglesort(coords_, i, pivot, tri)
                if findfirst(isequal(downadjs[1]), ordered) == 2
                    a = setstereo(nodes[i], :atypical)
                    setnodeattr!(mol, i, a)
                    continue
                end
            end
        elseif length(adjs) == 3
            # Implicit hydrogen is the pivot
            pivot = i
            tri = adjs
            quad = union(adjs, [i])
            rev = isempty(upadjs) ? false : true
        end
        cw = isclockwise(toarray(coords_, sort(tri)))
        if rev
            cw = !cw
        end
        # Arrange the configuration so that the lowest index node is the pivot.
        pividx = findfirst(isequal(pivot), sort(quad))
        if pividx in (2, 4)
            cw = !cw
        end
        stereo = cw ? :clockwise : :anticlockwise
        a = setstereo(nodes[i], stereo)
        setnodeattr!(mol, i, a)
    end
end
