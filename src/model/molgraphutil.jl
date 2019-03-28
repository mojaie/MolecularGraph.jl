#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    getindex,
    getatom, getbond, atomcount, bondcount,
    updateatom!, updatebond!, unlinkatom!, unlinkbond!


Base.getindex(mol::GeneralMol, sym::Symbol) = eval(Expr(:call, sym, mol))
Base.getindex(
    mol::GeneralMol, k1::Symbol, k2::Symbol, K::Symbol...
) = hcat(eval(Expr(:call, sym, mol)) for k in [k1, k2, K...])


# Aliases

getatom = getnode
getbond = getedge
atomcount = nodecount
bondcount = edgecount
updateatom! = updatenode!
updatebond! = updateedge!
unlinkatom! = unlinknode!
unlinkbond! = unlinkedge!
