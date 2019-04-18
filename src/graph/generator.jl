#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    pathgraph, cyclegraph,
    completebipartite, completegraph,
    laddergraph, circularladder, moebiusladder,
    squaregrid,
    generalizedpetersen, petersengraph, dodecahedralgraph


"""
    pathgraph(n::Int; mutable=false) -> PlainGraph

Generate path graph ``P_n``.
"""
function pathgraph(n::Int; mutable=false)
    n >= 2 || throw(DomainError(n, "n should be 2 or more"))
    f = mutable ? plaingraph : immutableplaingraph
    return f(n, (i, i + 1) for i in 1:(n-1))
end


"""
    cyclegraph(length::Int; mutable=false) -> PlainGraph

Generate cycle graph ``C_n``.
"""
function cyclegraph(n::Int; mutable=false)
    n >= 3 || throw(DomainError(n, "n should be 3 or more"))
    edges = [(i, i + 1) for i in 1:(n-1)]
    push!(edges, (n, 1))
    f = mutable ? plaingraph : immutableplaingraph
    return f(n, edges)
end


"""
    completebipartite(m::Int,n::Int; mutable=false) -> PlainGraph

Generate bipartite graph ``K_{m,n}``.
"""
function completebipartite(m::Int, n::Int; mutable=false)
    m >= 1 || throw(DomainError(m, "m should not be 1 or more"))
    n >= 1 || throw(DomainError(n, "n should not be 1 or more"))
    edges = Tuple{Int,Int}[]
    for i in 1:n
        for j in n+1:n+m
            push!(edges, (i, j))
        end
    end
    f = mutable ? plaingraph : immutableplaingraph
    return f(n + m, edges)
end


"""
    completegraph(length::Int; mutable=false) -> PlainGraph

Generate complete graph ``K_n``.
"""
function completegraph(n::Int; mutable=false)
    n >= 0 || throw(DomainError(n, "n should not be negative"))
    f = mutable ? plaingraph : immutableplaingraph
    n == 0 && return f()
    n == 1 && return f(1, [])
    return f(n, ((u, v) for (u, v) in combinations(1:n)))
end


"""
    laddergraph(n::Int; mutable=false) -> PlainGraph

Generate ladder graph ``L_n``.
"""
function laddergraph(n::Int; mutable=false)
    n >= 1 || throw(DomainError(n, "n should be 1 or more"))
    edges = [(1, 2)]
    for i in 1:n-1
        push!(edges, (2i+1, 2i+2))
        push!(edges, (2i-1, 2i+1))
        push!(edges, (2i, 2i+2))
    end
    f = mutable ? plaingraph : immutableplaingraph
    return f(n * 2, edges)
end


"""
    circularladder(n::Int; mutable=false) -> PlainGraph

Generate circular ladder graph ``CL_n``.
"""
function circularladder(n::Int; mutable=false)
    n >= 3 || throw(DomainError(n, "n should be 3 or more"))
    graph = laddergraph(n, mutable=true)
    addedge!(graph, 1, 2n-1)
    addedge!(graph, 2, 2n)
    if !mutable
        graph = immutableplaingraph(graph)
    end
    return graph
end


"""
    moebiusladder(n::Int; mutable=false) -> PlainGraph

Generate Möbius ladder graph ``ML_n``.
"""
function moebiusladder(n::Int; mutable=false)
    n >= 3 || throw(DomainError(n, "n should be 3 or more"))
    graph = laddergraph(n, mutable=true)
    addedge!(graph, 1, 2n)
    addedge!(graph, 2, 2n-1)
    if !mutable
        graph = immutableplaingraph(graph)
    end
    return graph
end


"""
    squaregrid(m::Int,n::Int; mutable=false) -> PlainGraph

Generate ``m \\times n`` square grid graph.

Use `cartesianproduct` for higher dimensional grid graphs.
"""
function squaregrid(m::Int, n::Int; mutable=false)
    m >= 1 || throw(DomainError(m, "m should not be 1 or more"))
    n >= 1 || throw(DomainError(n, "n should not be 1 or more"))
    edges = Tuple{Int,Int}[]
    for i in 1:(m*n-1)
        mod(i, n) == 0 && continue
        push!(edges, (i, i + 1))
    end
    for i in 1:((m-1)*n)
        push!(edges, (i, i + n))
    end
    f = mutable ? plaingraph : immutableplaingraph
    return f(n * m, edges)
end


"""
    generalizedpetersen(n::Int, k::Int; mutable=false) -> PlainGraph

Generate generalized petersen graph ``G(n,k)``.
"""
function generalizedpetersen(n::Int, k::Int; mutable=false)
    n >= 3 || throw(DomainError(n, "n should be 3 or more"))
    k < n / 2 || throw(DomainError(n, "k should be less than n/2"))
    edges = Tuple{Int,Int}[]
    for i in 0:(n-1)
        push!(edges, (i + 1, mod(i + 1, n) + 1))
        push!(edges, (i + 1, i + n + 1))
        push!(edges, (i + n + 1, mod(i + n + k, n) + n + 1))
    end
    f = mutable ? plaingraph : immutableplaingraph
    return f(2n, edges)
end

petersengraph(; kwargs...) = generalizedpetersen(5, 2; kwargs...)
dodecahedralgraph(; kwargs...) = generalizedpetersen(10, 2; kwargs...)
