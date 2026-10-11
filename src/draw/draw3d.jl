#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function atom_radius(mol::SimpleMolGraph; mapping=ATOM_VANDERWAALS_RADII)
    isa(mapping, Real) && return fill(mapping, nv(mol))
    desc = zeros(Float64, nv(mol))
    for i in vertices(mol)
        an = atom_number(mol[i])
        desc[i] = mapping[an]
        mapping === ATOM_COVALENT_RADII || continue
        isa(r, Real) && continue
        # Carbon and a few metals have multiple values
        if an == 6
            d = degree(mol, i)
            key = d == 3 ? "Csp3" :
                d == 2 ? "Csp2" : "Csp"
            desc[i] = r[key]
        else
            # For metals, it's safest to choose the smallest radius
            desc[i] = minimum(values(r))
        end
    end
    return desc
end
