#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    pathgraph, cyclegraph,
    bipartitegraph, completegraph,
    laddergraph, circularladder, moebiusladder


"""
    pathgraph(n::Int) -> VectorGraph{Node,Edge}

Generate path graph ``P_n``.
"""
function pathgraph(n::Int)
    n >= 2 || throw(DomainError(n, "n should be 2 or more"))
    return vectorgraph(n, (i, i + 1) for i in 1:(n-1))
end


"""
    cyclegraph(length::Int) -> VectorGraph{Node,Edge}

Generate cycle graph ``C_n``.
"""
function cyclegraph(n::Int)
    n >= 3 || throw(DomainError(n, "n should be 3 or more"))
    edges = [(i, i + 1) for i in 1:(n-1)]
    push!(edges, (n, 1))
    return vectorgraph(n, edges)
end


"""
    bipartitegraph(m::Int,n::Int) -> VectorGraph{Node,Edge}

Generate bipartite graph ``K_{m,n}``.
"""
function bipartitegraph(m::Int,n::Int)
    m >= 1 || throw(DomainError(m, "m should not be 1 or more"))
    n >= 1 || throw(DomainError(n, "n should not be 1 or more"))
    edges = Tuple{Int,Int}[]
    for i in 1:n
        for j in n+1:n+m
            push!(edges, (i, j))
        end
    end
    return vectorgraph(n + m, edges)
end


"""
    completegraph(length::Int) -> VectorGraph{Node,Edge}

Generate complete graph ``K_n``.
"""
function completegraph(n::Int)
    n >= 0 || throw(DomainError(n, "n should not be negative"))
    n == 0 && return mapgraph(Node,Edge)
    n == 1 && return mapgraph([1], [])
    return vectorgraph(n, ((u, v) for (u, v) in combinations(1:n)))
end


"""
    laddergraph(n::Int) -> VectorGraph{Node,Edge}

Generate ladder graph ``L_n``.
"""
function laddergraph(n::Int)
    n >= 1 || throw(DomainError(n, "n should be 1 or more"))
    edges = [(1, 2)]
    for i in 1:n-1
        push!(edges, (2i+1, 2i+2))
        push!(edges, (2i-1, 2i+1))
        push!(edges, (2i, 2i+2))
    end
    return vectorgraph(n * 2, edges)
end


"""
    circularladder(n::Int) -> VectorGraph{Node,Edge}

Generate circular ladder graph ``CL_n``.
"""
function circularladder(n::Int)
    n >= 3 || throw(DomainError(n, "n should be 3 or more"))
    graph = laddergraph(n)
    edgecount = 3n - 2
    updateedge!(graph, edgetype(graph)(1, 2n-1), edgecount + 1)
    updateedge!(graph, edgetype(graph)(2, 2n), edgecount + 2)
    return graph
end


"""
    moebiusladder(n::Int) -> VectorGraph{Node,Edge}

Generate Möbius ladder graph ``ML_n``.
"""
function moebiusladder(n::Int)
    n >= 3 || throw(DomainError(n, "n should be 3 or more"))
    graph = laddergraph(n)
    edgecount = 3n - 2
    updateedge!(graph, edgetype(graph)(1, 2n), edgecount + 1)
    updateedge!(graph, edgetype(graph)(2, 2n-1), edgecount + 2)
    return graph
end
