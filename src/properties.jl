#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


export
    molweight,
    hydrogen_acceptor_count,
    hydrogen_donor_count,
    wildman_crippen_logp,
    rotatable_count


"""
    molweight(mol::VectorMol; digits=2) -> Float64

Return standard molecular weight.
"""
function molweight(mol::VectorMol; digits=2)
    elemental!(mol)
    return round(reduce(+, mol[:MolWeight]; init=0), digits=digits)
end


"""
    hydrogen_acceptor_count(mol::VectorMol) -> Int

Return the number of hydrogen bond acceptors (N, O and F).
"""
function hydrogen_acceptor_count(mol::VectorMol)
    elemental!(mol)
    return reduce(+, mol[:H_Acceptor]; init=0)
end


"""
    hydrogen_donor_count(mol::VectorMol) -> Int

Return the number of hydrogen bond donors (O and N attached to hydrogens).
"""
function hydrogen_donor_count(mol::VectorMol)
    elemental!(mol)
    return reduce(+, mol[:H_Donor]; init=0)
end


"""
    wildman_crippen_logp(mol::VectorMol) -> Float64

Return predicted logP value calculated by using Wildman and Crippen method.

# Reference

1. Wildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical
   Parameters by Atomic Contributions. Journal of Chemical Information and
   Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l
"""
function wildman_crippen_logp(mol::VectorMol; digits=2)
    wclogpcalc!(mol)
    return round(reduce(+, mol[:WCLogPContrib]; init=0), digits=digits)
end


"""
    rotatable_count(mol::VectorMol) -> Int

Return the number of rotatable bonds.
"""
function rotatable_count(mol::VectorMol)
    rotatable!(mol)
    return reduce(+, mol[:Rotatable]; init=0)
end
