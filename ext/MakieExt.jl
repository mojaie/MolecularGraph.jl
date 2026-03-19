#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MakieExt

using MolecularGraph:
    MolecularGraph, SimpleMolGraph,
    spacefilling, spacefilling!,
    ballstick, ballstick!,
    stick, stick!,
    wire, wire!,
    atom_radius, coords3d, atom_symbol, bond_order,
    ATOM_VANDERWAALS_RADII, ATOM_COVALENT_RADII,
    atom_coloralpha, RASMOL_ATOM_COLOR

using Makie:
    Makie, @recipe, Theme, meshscatter!, mesh!, plot!

using GeometryBasics: Cylinder, Point
using Graphs: edges, neighbors, degree
using LinearAlgebra: cross, normalize


const DEFAULT_BALL_DIAMETER = float(0.4)
const DEFAULT_STICK_DIAMETER = float(0.33)
const DEFAULT_WIRE_DIAMETER = float(0.1)

const Z_DIR = [0, 0, 1]

const RADII_TYPE = Dict(
    "van der Waals" => ATOM_VANDERWAALS_RADII,
    "covalent" => ATOM_COVALENT_RADII
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


function Makie.plot!(md::MolDisplay{<:Tuple{<:SimpleMolGraph}})
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
        bonddiameter=DEFAULT_WIRE_DIAMETER, multiplebonds=false, kwargs...)
    order = multiplebonds ? bond_order(mol[e]) : 1
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


function MolecularGraph.spacefilling(args...; radii="van der Waals", kwargs...)
    moldisplay(args...; radii=get(RADII_TYPE, radii, radii), showbonds=false, kwargs...)
end

function MolecularGraph.spacefilling!(args...; radii="van der Waals", kwargs...)
    moldisplay!(args...; radii=get(RADII_TYPE, radii, radii), showbonds=false, kwargs...)
end

function MolecularGraph.ballstick(
        args...; radii=DEFAULT_BALL_DIAMETER, bonddiameter=DEFAULT_WIRE_DIAMETER, kwargs...)
    moldisplay(args...; radii=radii, bonddiameter=bonddiameter, multiplebonds=true, kwargs...)
end

function MolecularGraph.ballstick!(
        args...; radii=DEFAULT_BALL_DIAMETER, bonddiameter=DEFAULT_WIRE_DIAMETER, kwargs...)
    moldisplay!(args...; radii=radii, bonddiameter=bonddiameter, multiplebonds=true, kwargs...)
end

function MolecularGraph.stick(args...; size=DEFAULT_STICK_DIAMETER, kwargs...)
    moldisplay(args...; radii=size, bonddiameter=size, multiplebonds=false, kwargs...)
end

function MolecularGraph.stick!(args...; size=DEFAULT_STICK_DIAMETER, kwargs...)
    moldisplay!(args...; radii=size, bonddiameter=size, multiplebonds=false, kwargs...)
end

function MolecularGraph.wire(args...; size=DEFAULT_WIRE_DIAMETER, kwargs...)
    moldisplay(args...; bondwidth=size, showatoms=false, multiplebonds=true, kwargs...)
end

function MolecularGraph.wire!(args...; size=DEFAULT_WIRE_DIAMETER, kwargs...)
    moldisplay!(args...; bondwidth=size, showatoms=false, multiplebonds=true, kwargs...)
end


end  # module MakieExt
