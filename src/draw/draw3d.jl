#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const DEFAULT_BALL_DIAMETER = float(0.4)
const DEFAULT_STICK_DIAMETER = float(0.33)
const DEFAULT_WIRE_DIAMETER = float(0.1)

const Z_DIR = [0,0,1]

function atom_radius(mol::SimpleMolGraph; mapping=ATOM_VANDERWAALS_RADII)
    isa(mapping, Real) && return fill(mapping, nv(mol))
    desc = zeros(Float64, nv(mol))
    for i in vertices(mol)
        an = atom_number(props(mol, i))
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


const RADII_TYPE = Dict(
    "van der Waals" => ATOM_VANDERWAALS_RADII,
    "covalent" => ATOM_COVALENT_RADII
)

"""
    spacefilling(mol::MolGraph; radii="van der Waals")

Represent `mol` as a space-filling (Calotte) model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. (3D SDF files can be downloaded from sites such as PubChem.) The two supported options for `radii` are
`"van der Waals"` and `"covalent"`; the former are available only for main-group elements, and the latter are available for
all.

This function requires that you load one of the backends of the Makie/GLMakie/CairoMakie family.
"""
spacefilling(
    args...;
    radii="van der Waals",
    kwargs...
) = moldisplay(
    args...;
    radii=get(RADII_TYPE, radii, radii),
    showbonds=false,
    kwargs...
)

spacefilling!(
    args...;
    radii="van der Waals",
    kwargs...
) = moldisplay!(
    args...;
    radii=get(RADII_TYPE, radii, radii),
    showbonds=false,
    kwargs...
)



"""
    ballstick(mol::UndirectedGraph; radii=0.3, bonddiameter=0.1)

Represent `mol` as a ball-and-stick model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. 3D SDF files can be downloaded from sites such as PubChem.

`radii` optionally specifies the radii of the balls, in Angstroms.
`bonddiameter` optionally specifies the radii of the sticks, in Angstroms.

This function requires that you load one of the backends of the Makie/GLMakie/CairoMakie family.
"""
ballstick(
    args...;
    radii=DEFAULT_BALL_DIAMETER,
    bonddiameter=DEFAULT_WIRE_DIAMETER,
    kwargs...
) = moldisplay(
    args...;
    radii=radii,
    bonddiameter=bonddiameter,
    multiplebonds=true,
    kwargs...
)

ballstick!(
    args...;
    radii=DEFAULT_BALL_DIAMETER, 
    bonddiameter=DEFAULT_WIRE_DIAMETER,
    kwargs...
) = moldisplay!(
    args...;
    radii=radii,
    bonddiameter=bonddiameter,
    multiplebonds=true, 
    kwargs...
)


"""
    stick(mol::UndirectedGraph; size=0.3)

Represent `mol` as a stick model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. 3D SDF files can be downloaded from sites such as PubChem.

`size` optionally specifies the width of the sticks, in Angstroms.

This function requires that you load one of the backends of the Makie/GLMakie/CairoMakie family.
"""
stick(
    args...;
    size=DEFAULT_STICK_DIAMETER,
    kwargs...
) = moldisplay(
    args...;
    radii=size,
    bonddiameter=size,
    multiplebonds=false,
    kwargs...
)

stick!(
    args...;
    size=DEFAULT_STICK_DIAMETER,
    kwargs...
) = moldisplay!(
    args...;
    radii=size,
    bonddiameter=size,
    multiplebonds=false,
    kwargs...
)


"""
    wire(mol::UndirectedGraph; size=0.1)

Represent `mol` as a wire-frame model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. 3D SDF files can be downloaded from sites such as PubChem.

`size` optionally specifies the width of the bonds, in Angstroms.

This function requires that you load one of the backends of the Makie/GLMakie/CairoMakie family.
"""
wire(
    args...;
    size=DEFAULT_WIRE_DIAMETER,
    kwargs...
) = moldisplay(
    args...;
    bondwidth=size,
    showatoms=false,
    multiplebonds=true,
    kwargs...
)

wire!(
    args...;
    size=DEFAULT_WIRE_DIAMETER,
    kwargs...
) = moldisplay!(
    args...;
    bondwidth=size,
    showatoms=false,
    multiplebonds=true,
    kwargs...
)


@recipe(MolDisplay, mol) do scene
    Theme(
        radii=DEFAULT_BALL_DIAMETER,
        bonddiameter=DEFAULT_WIRE_DIAMETER,
        colortheme=RASMOL_ATOM_COLOR,
        multiplebonds=true,
        showbonds=true,
        showatoms=true,
        alpha=1.0,
    )
end


function MakieCore.plot!(md::MolDisplay{<:Tuple{<:SimpleMolGraph}})
    mols = [md[i][] for i=1:length(md)]
    radii = md.radii[]
    for mol in mols
        crds = coords3d(mol)
        col = atom_coloralpha(mol, alpha=md.alpha[], color_theme=md.colortheme[])
        if md.showatoms[]
            rd = atom_radius(mol; mapping=radii)
            drawatoms!(md, crds, col, rd)
        end
        if md.showbonds[]
            syms = atom_symbol(mol)
            nbrs = degree(mol)
            for e in edges(mol)
                drawbond!(
                    md, mol, e, crds, col, syms, nbrs;
                    bonddiameter=md.bonddiameter[], multiplebonds=md.multiplebonds[])
            end
        end
    end
    return md
end


drawatoms!(f, crds, col, rd; kwargs...) = meshscatter!(
    f, [c[1] for c in crds], [c[2] for c in crds], [c[3] for c in crds];
    color=col, markersize=rd)

function drawbond!(
        f, mol::SimpleMolGraph, e, crds, col, syms, nbrs;
        bonddiameter=DEFAULT_BOND_DIAMETER, multiplebonds=false, kwargs...)
    order = multiplebonds ? bond_order(props(mol, e)) : 1
    atomidx1, atomidx2 = e.src, e.dst
    pos1, pos2 = crds[atomidx1], crds[atomidx2]
    normaldir = Z_DIR
    if order > 1
        ng1, ng2 = nbrs[atomidx1], nbrs[atomidx2]
        # determine the plane for double bonds
        if ng1 == 3
            neighs = filter(x -> x != atomidx2, (neighbors(mol, atomidx1)))
            @assert length(neighs) == 2
            npos1, npos2 = crds[neighs[1]], crds[neighs[2]]
            normaldir = cross(npos1, npos2)
        elseif ng2 == 3
            neighs = filter(x -> x != atomidx1, (neighbors(mol, atomidx2)))
            @assert length(neighs) == 2
            npos1, npos2 = crds[neighs[1]], crds[neighs[2]]
            normaldir = cross(npos1, npos2)
        end
    end
    sepdir = normalize(cross(normaldir, pos2 .- pos1))
    dists = (bonddiameter * 2.5) .* collect(-0.5 * (order-1): 0.5 * (order-1))
    for dist in dists
        dvec = dist * sepdir
        p1, p2 = pos1 + dvec, pos2 + dvec
        if syms[atomidx1] == syms[atomidx2]
            cyl = Cylinder(p1, p2, bonddiameter)
            mesh!(f, cyl; color=col[atomidx1], kwargs...)
        else
            midpoint = 0.5 * (p1 + p2)
            pm = Point(midpoint...)
            cyl1 = Cylinder(p1, pm, bonddiameter)
            cyl2 = Cylinder(p2, pm, bonddiameter)
            mesh!(f, cyl1; color=col[atomidx1], kwargs...)
            mesh!(f, cyl2; color=col[atomidx2], kwargs...)
        end
    end
    return f
end
