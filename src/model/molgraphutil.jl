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


Base.getindex(mol::VectorMol, k::Symbol) = mol.vector[k]
Base.getindex(
    mol::VectorMol, k1::Symbol, k2::Symbol, K::Symbol...
) = AnnotationArray(
    hcat([mol.vector[k] for k in [k1, k2, K...]]...), [k1, k2, K...])

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
