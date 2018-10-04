

mutable struct UndirectedGraph
    adjacency::Dict
    nodes::Dict
    edges::Dict
    function UndirectedGraph(nodes, edges)
        initialize!(new(), nodes, edges)
    end
end


UndirectedGraph() = UndirectedGraph([], [])


function initialize!(graph::UndirectedGraph,
                     nodes::Array{Union{Any, Any}},
                     edges::Array{Union{Any, Any, Any}})
    for node in nodes
        updatenode!(graph, node)
    end
    for edge in edges
        updateedge!(graph, edge)
    end
    graph
end


function getnode(graph::UndirectedGraph, idx)
    graph.nodes[idx]
end


function getedge(graph::UndirectedGraph, u, v)
    graph.adjacency[u][v]
end


function updatenode!(graph::UndirectedGraph, node::Tuple{Any, Any})
    updatenode!(graph, node[0], node[1])
end


function updatenode!(graph::UndirectedGraph, idx, attr)
    graph.nodes[idx] = attr
    return
end


function updateedge!(graph::UndirectedGraph, edge::Tuple{Any, Any, Any})
    updateedge!(graph, edge[0], edge[1], edge[2])
end


function updateedge!(graph::UndirectedGraph, u, v, attr)
    graph.adjacency[u][v] = attr
    graph.adjacency[v][u] = attr
    return
end
