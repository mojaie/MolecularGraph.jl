
abstract type Graph end


struct UndirectedGraph <: Graph
    adjacency::Dict
    nodes::Array
    function UndirectedGraph(nodes::Array, edges::Array)
        new(adjacency, nodes)
    end
end


function node(graph::Graph, idx)
    graph.nodes[idx]
end


function edge(graph::Graph, u, v)
    graph.adjacency[u][v]
end


function updatenode!(graph::Graph, idx::Int32, attr::Dict)
    graph.nodes[idx] = attr
end


function updateedge!(graph::Graph, u::Int32, v::Int32, attr::Dict)
    graph.adjacency[u][v] = attr
    graph.adjacency[v][u] = attr
end
