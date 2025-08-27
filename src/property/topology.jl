
# Molecular graph topology descriptors

"""
    degree(mol::SimpleMolGraph{T}) -> Vector{T}

Return a vector of size ``n`` representing the node degree of the molecular graph
of 1 to ``n``th atoms of the given molecule.

This property corresponds to SMARTS `D` query.
"""
Graphs.degree


"""
    sssr(mol::SimpleMolGraph{T}) -> Vector{Vector{T}}

Return vectors of ring nodes representing small set of smallest rings (SSSR).

See [`mincyclebasis`](@ref).
"""
function sssr(mol::SimpleMolGraph)
    if has_descriptor(mol, :sssr)
        return get_descriptor(mol, :sssr)
    end
    return mincyclebasis(mol.graph)
end

function sssr!(mol::SimpleMolGraph)
    mol.state.has_new_edges || return  # skip if new edges
    set_descriptor!(mol, :sssr, mincyclebasis(mol.graph))
    return
end

function remap!(
        ::Val{:sssr}, desc::MolDescriptor{T}, vmap::Vector{T},
        edges::Vector{Edge{T}}) where T <: Integer
    # Just remove incomplete cycles after `rem_vertex!`
    # to avoid expensive recalculation of myncycles
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    cont = Vector{T}[]
    oldmems = keys(revv)
    for cyc in desc.sssr
        isempty(setdiff(cyc, oldmems)) || continue
        push!(cont, [revv[v] for v in cyc])
    end
    empty!(desc.sssr)
    append!(desc.sssr, cont)
    return
end


"""
    which_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th nodes of the given graph.

SSSR membership is represented as a vector of SSSR indices assigned to each rings.
This means nodes that have the same SSSR index belong to the same SSSR.
"""
function which_ring(mol::SimpleMolGraph)
    arr = [Int[] for _ in vertices(mol)]
    for (i, cyc) in enumerate(sssr(mol))
        for n in cyc
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    edge_which_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`sssr`](@ref) membership of
1 to ``n``th bonds of the given molecule.

SSSR membership is represented as a set of SSSR indices assigned to each rings.
This means bonds that have the same SSSR index belong to the same SSSR.
"""
function edge_which_ring(mol::SimpleMolGraph{T}) where T
    arr = [Int[] for _ in 1:ne(mol)]
    ernk = edge_rank(mol)
    for (i, cyc) in enumerate(sssr(mol))
        for j in 1:(length(cyc) - 1)
            push!(arr[edge_rank(ernk, cyc[j], cyc[j + 1])], i)
        end
        push!(arr[edge_rank(ernk, cyc[1], cyc[end])], i)
    end
    return arr
end


"""
    fused_rings(mol::SimpleMolGraph{T}) -> Vector{Vector{T}}

Return vectors of fused ring node sets.

A fused ring is defined as a 2-edge connected components in terms of graph theory.
Spirocyclic structures are considered to be part of a fused ring.
"""
function fused_rings(g::SimpleGraph)
    cobr = setdiff(Set(edges(g)), bridges(g))
    subg, vmap = induced_subgraph(g, collect(cobr))
    return  [vmap[c] for c in connected_components(subg)]
end

fused_rings(mol::SimpleMolGraph) = fused_rings(mol.graph)


"""
    which_fused_ring(mol::SimpleMolGraph) -> Vector{Vector{Int}}

Return a vector of size ``n`` representing [`fused_rings`](@ref) membership of
1 to ``n``th atoms of the given molecule.

Fused ring membership is represented as a set of fused ring indices assigned to each fused rings.
This means atoms that have the same fused ring index belong to the same fused ring.
"""
function which_fused_ring(mol::SimpleMolGraph)
    arr = [Int[] for _ in vertices(mol)]
    for (i, conn) in enumerate(fused_rings(mol))
        for n in conn
            push!(arr[n], i)
        end
    end
    return arr
end


"""
    smallest_ring(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the size of the smallest [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to.

If the node is not in a ring, the value would be 0.
This property corresponds to SMARTS `r` query.
"""
function smallest_ring(mol::SimpleMolGraph)
    sssr_ = sssr(mol)
    whichring_ = which_ring(mol)
    arr = zeros(Int, nv(mol))
    for i in vertices(mol)
        rs = whichring_[i]
        isempty(rs) && continue
        arr[i] = minimum(length, sssr_[rs])
    end
    return arr
end


"""
    ring_count(mol::SimpleMolGraph) -> Vector{Int}

Return a vector of size ``n`` representing the number of [`sssr`](@ref)
that 1 to ``n``th atoms of the given molecule belong to.

This property corresponds to SMARTS `R` query.
"""
ring_count(mol::SimpleMolGraph) = length.(which_ring(mol))


"""
    is_in_ring(mol::SimpleMolGraph) -> BitVector

Return a vector of size ``n`` representing whether 1 to ``n``th atoms of
the given molecule belong to a ring or not.
"""
is_in_ring(mol::SimpleMolGraph) = .!isempty.(which_ring(mol))


"""
    is_edge_in_ring(mol::SimpleMolGraph) -> BitVector

Return a vector of size ``n`` representing whether 1 to ``n``th bonds of
the given molecule belong to a ring or not.
"""
is_edge_in_ring(mol::SimpleMolGraph) = .!isempty.(edge_which_ring(mol))



"""
    bmscaffold(mol::SimpleMolGraph) -> Vector{Int}

Return Bemis-Murcko scaffold as a vector of scaffold vertices.
"""
function bmscaffold(g::SimpleGraph)
    (length(connected_components(g)) == 1
        || error("bmscaffold for a disconnected graph cannot be determined"))
    degree_ = degree(g)
    while true
        terms = findall(degree_ .== 1)
        isempty(terms) && break
        for t in terms
            degree_[t] = 0
            for nbr in neighbors(g, t)
                degree_[nbr] == 0 && continue
                degree_[nbr] -= 1
            end
        end
    end
    return findall(degree_ .!= 0)
end
bmscaffold(mol::SimpleMolGraph) = bmscaffold(mol.graph)


"""
    scaffold_fragments(g::SimpleGraph) -> Vector{Int}

Return scaffold fragments. Intended for a MMP derivative.
"""
function scaffold_fragments(g::SimpleGraph)
    bms = bmscaffold(g)
    isempty(bms) && return scaffolds  # no scaffold
    queue = [bms]  # BFS
    scaffolds = Vector{Int}[bms]
    dag = SimpleDiGraph{Int}(1)  # fragmentation relationship (parent -> child)
    while !isempty(queue)
        # generate new fragments
        scaffold = pop!(queue)
        scf, scfvmap = induced_subgraph(g, scaffold)
        parent =findfirst(x -> x == scaffold, scaffolds)
        for fr in fused_rings(scf)
            # enumerate connections to the fused ring
            nbrs = setdiff(union([neighbors(scf, v) for v in fr]...), fr)
            length(nbrs) == 1 || continue
            # remove a fused ring
            subg, vmap = induced_subgraph(g, setdiff(scaffold, scfvmap[fr]))
            fragment = vmap[bmscaffold(subg)]
            if fragment in scaffolds  # duplicate, skip the branch
                dst = findfirst(x -> x == fragment, scaffolds)
                add_edge!(dag, parent, dst)
                continue
            end
            pushfirst!(queue, fragment)
            push!(scaffolds, fragment)
            add_vertex!(dag)
            add_edge!(dag, parent, length(scaffolds))
            if length(scaffolds) == 400
                @info "max fragment count reached"
                break
            end
        end
    end
    return scaffolds, dag
end
