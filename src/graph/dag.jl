#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    DAG,
    ancestors,
    descendants,
    topologicalsort


struct DAG{N<:AbstractNode,E<:AbstractDirectedEdge} <: MapDGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    successors::Dict{Int,Dict{Int,Int}}
    predecessors::Dict{Int,Dict{Int,Int}}

    function DAG{N,E}() where {N<:AbstractNode,E<:AbstractDirectedEdge}
        new(Dict(), Dict(), Dict(), Dict())
    end
end

function DAG(nodes::AbstractArray{Int},
             edges::AbstractArray{Tuple{Int,Int}})
    graph = DAG{Node,Arrow}()
    for node in nodes
        updatenode!(graph, Node(), node)
    end
    for (i, edge) in enumerate(edges)
        updateedge!(graph, Arrow(edge...), i)
    end
    graph
end


function ancestors(graph::DAG, idx::Int)
    rev = ReverseGraph(graph)
    return reachablenodes(rev, idx)
end


function descendants(graph::DAG, idx::Int)
    return reachablenodes(graph, idx)
end


function topologicalsort(graph::DAG)
    result = Int[]
    stack = [i for i in nodekeys(graph) if indegree(graph, i) == 0]
    used = Set{Int}()
    while !isempty(stack)
        n = pop!(stack)
        push!(result, n)
        for (s, edge) in successors(graph, n)
            if edge in used
                continue
            end
            push!(used, edge)
            if isempty(setdiff(inedgekeys(graph, s), used))
                push!(stack, s)
            end
        end
    end
    if edgecount(graph) != length(used)
        throw(OperationError("Cycle found"))
    end
    return result
end
