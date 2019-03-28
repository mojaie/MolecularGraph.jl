#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


struct SubstructureView{T<:VectorGraph} <: GeneralMol
    graph::SubgraphView{T}
    attribute::Dict{Symbol,Any}
    cache::Dict{Symbol,Any}
end


function MolecularGraphModel.nodesubgraph(mol::VectorMol, nodes)
    subg = nodesubgraph(mol.graph, nodes)
    return SubstructureView(subg, mol.attribute, mol.cache)
end


function MolecularGraphModel.edgesubgraph(mol::VectorMol, edges)
    subg = edgesubgraph(mol.graph, edges)
    return SubstructureView(subg, mol.attribute, mol.cache)
end
