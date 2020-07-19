#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    mincycles, edgemincycles, mincyclemembership, edgemincyclemembership,
    canonicalcycle


"""
    mincycles(graph::UndirectedGraph) -> Vector{Vector{Int}}

Calculate minimum cycle basis (also known as Smallest Set of Smallest Rings
in the context of molecular graph theory).
"""
@cache function edgemincycles(graph::UndirectedGraph)
    mincycs = Vector{Int}[]
    # TODO: twoedgeconnectedcomponents may not improve performance significantlly
    for biconn in twoedgeconnectedcomponents(graph)
        subg = nodesubgraph(graph, biconn)
        for cyc in mincyclebasis(subg)
            push!(mincycs, cyc)
        end
    end
    return mincycs
end

edgemincycles(view::SubgraphView) = edgemincycles(view.graph)


@cache function mincycles(graph::UndirectedGraph)
    mincycs = Vector{Int}[]
    for cy in edgemincycles(graph)
        cycy = vcat(cy, cy)
        nodes = Int[]
        for i in 1:length(cy)
            a = getedge(graph, cycy[i])
            b = getedge(graph, cycy[i + 1])
            push!(nodes, iterate(intersect(a, b))[1])
        end
        push!(mincycs, nodes)
    end
    return mincycs
end

mincycles(view::SubgraphView) = mincycles(view.graph)


function mincyclemembership(graph::OrderedGraph)
    nodes = [Set{Int}() for n in 1:nodecount(graph)]
    for (i, cyc) in enumerate(mincycles(graph))
        for n in cyc
            push!(nodes[n], i)
        end
    end
    return nodes
end

mincyclemembership(view::SubgraphView) = mincyclemembership(view.graph)


function edgemincyclemembership(graph::OrderedGraph)
    edges = [Set{Int}() for n in 1:edgecount(graph)]
    for (i, cyc) in enumerate(edgemincycles(graph))
        for e in cyc
            push!(edges[e], i)
        end
    end
    return edges
end

edgemincyclemembership(view::SubgraphView) = edgemincyclemembership(view.graph)


function mincyclebasis(graph::UndirectedGraph)
    # de Pina algorithm re-interpreted by Kavitha et al.
    cycles = Vector{Int}[]
    root = pop!(nodeset(graph))
    S = [Set(e) for e in cotree_edges(graph, root)]
    N = length(S)  # N: circuit rank
    for k in 1:N
        minedges = findmincycle(graph, S[k])
        push!(cycles, minedges)
        for i in (k + 1):N
            if length(intersect(S[i], minedges)) % 2 == 1
                S[i] = symdiff(S[i], S[k])
            end
        end
    end
    return cycles
end


function findmincycle(graph, S)
    G = edgesubgraph(graph, setdiff(edgeset(graph), S))
    U = disjointunion(G, G)
    for s in S
        (u, v) = getedge(graph, s)
        u1 = getunionnode(U, 1, u)
        v1 = getunionnode(U, 1, v)
        u2 = getunionnode(U, 2, u)
        v2 = getunionnode(U, 2, v)
        addedge!(U, u1, v2, DisjointUnionEdge(1, s))
        addedge!(U, u2, v1, DisjointUnionEdge(1, s))
    end
    paths = []
    for n in nodeset(graph)
        n1 = getunionnode(U, 1, n)
        n2 = getunionnode(U, 2, n)
        sp = shortestpathedges(U, n1, n2)
        push!(paths, sp)
    end
    minpath = sortstablemin(paths, by=length)
    return [getsourceedge(U, e).sourcekey for e in minpath]
end


function canonicalcycle(nodes::Vector{Int})
    """Align cycle indices to start from lowest index and following one of
    neighbors that have the lower index
    """
    # TODO: not used
    (fst, fstidx) = findmin(nodes)
    succidx = fstidx == lastindex(nodes) ? 1 : fstidx + 1
    succ = nodes[succidx]
    predidx = fstidx == 1 ? lastindex(nodes) : fstidx - 1
    pred = nodes[predidx]
    cp = succ < pred ? copy(nodes) : reverse(copy(nodes))
    while cp[1] != fst
        pushfirst!(cp, pop!(cp))
    end
    return cp
end
