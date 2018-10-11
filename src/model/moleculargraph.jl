#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolecularGraph,
    getatom,
    atompos,
    getbond,
    neighbors,
    atomvector,
    bondvector,
    adjvector,
    newatom!,
    updatebond!,
    required_descriptor

import ..GraphModel: neighbors


mutable struct MolecularGraph
    graph::UndirectedGraph
    descriptors::Set
    rings
    scaffolds
    isolated
    data
    function MolecularGraph()
        mol = new()
        mol.graph = UndirectedGraph{UInt16}()
        mol.descriptors = Set()
        mol
    end
end


function getatom(mol::MolecularGraph, idx)
    getnode(mol.graph, idx)
end

function atompos(mol::MolecularGraph, idx)
    nodepos(mol.graph, idx)
end

function getbond(mol::MolecularGraph, u, v)
    getedge(mol.graph, u, v)
end


function neighbors(mol::MolecularGraph, idx)
    neighbors(mol.graph, idx)
end


function atomvector(mol::MolecularGraph)
    mol.graph.nodes
end


function bondvector(mol::MolecularGraph)
    mol.graph.edges
end


function adjvector(mol::MolecularGraph)
    mol.graph.adjacency
end


function newatom!(mol::MolecularGraph, atom::Atom)
    newnode!(mol.graph, atom)
end

function newatom!(mol::MolecularGraph, idx::Integer, atom::Atom)
    atom.index = idx
    newatom!(mol, atom)
end


function updatebond!(mol::MolecularGraph, bond::Bond)
    updateedge!(mol.graph, bond)
end

function updatebond!(mol::MolecularGraph, u::Integer, v::Integer, bond::Bond)
    bond.u = u
    bond.v = v
    updatebond!(mol, bond)
end


function required_descriptor(mol::MolecularGraph, desc::AbstractString)
    if desc ∉ mol.descriptors
        throw(ErrorException("$(desc) is not assigned"))
    end
end
