#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    addchiralhydrogens, removechiralhydrogens,
    chiralcenter, setdiastereo!, setchiralcenter!


"""
    addchiralhydrogens(mol::SMILES) -> SMILES

Return new molecule with explicit chiral hydrogen nodes.
"""
function addchiralhydrogens(mol::SMILES)
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
    removechiralhydrogens(mol::SMILES) -> SMILES

Return new molecule without explicit chiral hydrogen nodes.
"""
function removechiralhydrogens(mol::SMILES)
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
    chiralcenter(mol::SMILES) -> SMILES

Atom mapping for chirality information (chiral center => tuple of four nodes). `:unspecified`: no chirality information, `:clockwise` and `:anticlockwise`: If we see the chiral center from the 1st node, 2-4th nodes are placed clockwise/anticlockwise. If there is an implicit hydrogen, it will be shown as the pseudo node index of -1.
"""
@cache function chiralcenter(mol::SMILES)
    # TODO: unnecessary?
    dict = Dict{Int,NTuple{4,Int}}()
    atoms = nodeattrs(mol)
    implicith_ = implicithcount(mol)
    for i in 1:nodecount(mol)
        atoms[i].stereo === :unspecified && continue
        adj = collect(adjacencies(mol, i))
        if implicith_[i] == 1
            a1, a2, a3 = sort(adj)
            if atoms[i].stereo === :anticlockwise
                a2, a3 = (a3, a2)
            end
            dict[i] = i == 1 ? (-1, a1, a2, a3) : (a1, -1, a2, a3)
        else
            a1, a2, a3, a4 = sort(adj)
            if atoms[i].stereo === :anticlockwise
                a3, a4 = (a4, a3)
            end
            dict[i] = (a1, a2, a3, a4)
        end
    end
    return dict
end


"""
    setdiastereo!(mol::SMILES) -> SMILES

Set diastereomer infomation to double bonds.

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


function setchiralcenter!(mol::SDFile)
    nodes = nodeattrs(mol)
    edges = edgeattrs(mol)
    coords_ = coords2d(mol)
    for i in 1:nodecount(mol)
        adjs = Int[]
        upadjs = Int[]
        downadjs = Int[]
        for (inc, adj) in neighbors(mol, i)
            push!(adjs, adj)
            mol.edges[inc][1] == i || continue
            if edges[inc].notation == 1
                push!(upadjs, adj)
            elseif edges[inc].notation == 6
                push!(downadjs, adj)
            end
        end
        (isempty(upadjs) && isempty(downadjs)) && continue  # unspecified
        if length(adjs) != 4
            # TODO: hypervalent, syn/anti or axial chirality
            a = setstereo(nodes[i], :atypical)
            setnodeattr!(mol, i, a)
            continue
        end
        # Select a pivot
        rev = false
        if !isempty(upadjs)
            pivot = upadjs[1]
        elseif length(downadjs) == 1
            pivot = downadjs[1]
            rev = true
        else
            pivot = setdiff(adjs, downadjs)[1]
        end
        # Remove pivot and check direction
        sadjs = sort(adjs, by=p->p[1])
        triplet = setdiff(sadjs, [pivot])
        cw = isclockwise(cartesian2d(coords_, triplet))
        if rev
            cw = !cw
        end
        # Arrange the configuration so that the lowest index node is the pivot.
        pividx = findfirst(isequal(pivot), sadjs)
        if pividx in (2, 4)
            cw = !cw
        end
        a = setstereo(nodes[i], cw ? :clockwise : :anticlockwise)
        setnodeattr!(mol, i, a)
    end
end
