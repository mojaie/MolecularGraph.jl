#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    getindex, haskey, setindex!,
    getatom, getnode, getbond, getedge,
    neighbors, neighborcount, degree,
    nodecount, atomcount, edgecount, bondcount,
    updatenode!, updateatom!, updatebond!, updatebond!,
    unlinknode!, unlinkatom!, unlinkedge!, unlinkbond!


import Base: getindex, haskey, setindex!


Base.getindex(mol::VectorMol, k) = mol.vector[k]
Base.haskey(mol::VectorMol, k) = haskey(mol.vector, k)

function Base.setindex!(mol::VectorMol, v, k)
    mol.vector[k] = v
end


# Aliases

getatom = getnode
getbond = getedge
atomcount = nodecount
bondcount = edgecount
updateatom! = updatenode!
updatebond! = updateedge!
unlinkatom! = unlinknode!
unlinkbond! = unlinkedge!
