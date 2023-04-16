#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

import MakieCore: plot!
import GeometryBasics: Point

using MakieCore: @recipe, Theme, meshscatter!, lines!, mesh!
using Colors: RGB
using GeometryBasics: mesh, Cylinder

export atom_radius

colortype(c::Color) = RGB{Float32}(c.r/255, c.g/255, c.b/255)


function atom_radius(mol::SimpleMolGraph; mapping=ATOM_VANDERWAALS_RADII)
    desc = init_node_descriptor(Float64, mol)
    for i in vertices(mol)
        an = atomnumber(get_prop(mol, i, :symbol))
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


"""
    spacefilling(mol::SimpleMolGraph; radius="van der Waals", tform=identity, H::Bool=true)

Makie.jl recipe for drawing 3D space filling model.

Represent `mol` as a space-filling (Calotte) model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. (3D SDF files can be downloaded from sites such as PubChem.) The two supported options for `radius` are
`"van der Waals"` and `"covalent"`; the former are available only for main-group elements, and the latter are available for
all.


This function requires that you load one of the backends of the GLMakie/WGLMakie/CairoMakie family.

example:

scene = Scene(camera = cam3d!, show_axis = false)
"""
@recipe(SpaceFilling, mol) do scene
    Theme(
        colortheme=RASMOL_ATOM_COLOR,
        radii=ATOM_VANDERWAALS_RADII
    )
end

function plot!(sc::SpaceFilling{<:Tuple{<:SimpleMolGraph}})
    mol = sc[1][]
    crds = coords3d(mol)
    col = colortype.(atom_color(mol, color_theme=sc.colortheme[]))
    rd = atom_radius(mol, mapping=sc.radii[])
    meshscatter!(sc, crds[:, 1], crds[:, 2], crds[:, 3];
        color=col,
        markersize=rd
    )
end


"""
    ballstick(mol::SimpleMolGraph; tform=identity, H::Bool=true, markersize=1)

Makie.jl recipe for drawing 3D ball-and-stick model.

Represent `mol` as a ball-and-stick model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. 3D SDF files can be downloaded from sites such as PubChem.

`tform(xyz)` optionally transforms the positions before plotting them. `H=false` causes hydrogens to be omitted from the plot.
`markersize` specifies the diameter of the balls, in Angstroms.

This function requires that you load one of the backends of the GLMakie/WGLMakie/CairoMakie family.
"""
@recipe(BallStick, mol) do scene
    Theme(
        colortheme=RASMOL_ATOM_COLOR,
        bonddiameter=0.2
    )
end

function plot!(sc::BallStick{<:Tuple{<:SimpleMolGraph}})
    mol = sc[1][]
    crds = coords3d(mol)
    col = colortype.(atom_color(mol, color_theme=sc.colortheme[]))
    meshscatter!(sc, crds[:, 1], crds[:, 2], crds[:, 3];
        color=col,
        markersize=0.4
    )
    diam = sc.bonddiameter[]
    for e in edges(mol)
        b = Cylinder(
            Point(crds[src(e), :]...),
            Point(crds[dst(e), :]...), float(0.1*get_prop(mol, e, :order)))
        me = mesh(b)
        mesh!(sc, me, color=:gray)
        """
        lines!(sc,
            crds[[src(e), dst(e)], 1],
            crds[[src(e), dst(e)], 2],
            crds[[src(e), dst(e)], 3],
            linewidth=10*get_prop(mol, e, :order)
        )
        """
    end
end

# TODO: double/triple bond drawing
# TODO: bond color corresponds to connecting atom
"""
meshscatter!(sc, crds[:, 1], crds[:, 2], crds[:, 3];
        color=atomcolor,
        markersize=0.3
    )

sc = Scene(camera = cam3d!, show_axis = false)
center!(sc)


meshscatter(crds[:, 1], crds[:, 2], crds[:, 3];
    color=atomcolor,
    markersize=radius,
    axis=(;show_axis=false)
)

sc = Scene(camera=cam3d!, show_axis=false)
a = Cylinder(GeometryBasics.Point(0.0, 0.0, 0.0), GeometryBasics.Point(1.0, 1.0, 1.0), 0.2)
mesh = GeometryBasics.mesh(a)
mesh!(sc, mesh, color=:gray)
mol = sdftomol(joinpath(@__DIR__, "assets", "test", "aspirin_3d.sdf"))
"""
