#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    addchiralhydrogens, removechiralhydrogens,
    chiralcenter, diastereobond


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
    diastereobond(mol::SMILES) -> SMILES

Bond mapping for diastereomeric bond information (diastereo bond => tuple of two edges and a symbol (`:unspecified`, `:cis` or `:trans`).
"""
@cache function diastereobond(mol::SMILES)
    dict = Dict{Int,Tuple{Int,Int,Symbol}}()
    edges = edgeattrs(mol)
    for i in 1:edgecount(mol)
        edges[i].order == 2 || continue
        u, v = mol.edges[i]
        ui = -1
        us = :unspecified
        for i in incidences(mol, u)
            edges[i].stereo === :unspecified && continue
            ui = i
            us = edges[i].stereo
            break
        end
        us === :unspecified && continue
        vi = -1
        vs = :unspecified
        for i in incidences(mol, v)
            edges[i].stereo === :unspecified && continue
            vi = i
            vs = edges[i].stereo
            break
        end
        vs === :unspecified && continue
        cistrans = us === vs ? :trans : :cis
        dict[i] = (ui, vi, cistrans)
    end
    return dict
end
