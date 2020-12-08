#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# This is loaded via @require

using .AbstractPlotting: RGB, N0f8, Vec3f0, Scene, SceneSpace, Axis, meshscatter, scatter, lines!

export spacefilling, ballstick

colortype(c::Color) = RGB{N0f8}(c.r/255, c.g/255, c.b/255)

"""
    spacefilling(mol::UndirectGraph; radius="van der Waals", tform=identity, H::Bool=true)

Represent `mol` as a space-filling (Calotte) model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. (3D SDF files can be downloaded from sites such as PubChem.) The two supported options for `radius` are
`"van der Waals"` and `"covalent"`; the former are available only for main-group elements, and the latter are available for
all.

`tform(xyz)` optionally transforms the positions before plotting them. `H=false` causes hydrogens to be omitted from the plot.

This function requires that you load one of the backends of the Makie/WGLMakie/CairoMakie family.
"""
function spacefilling(mol::UndirectedGraph; tform=identity, H::Bool=true)
    syms = [a.symbol for a in nodeattrs(mol)]
    pos = reduce(hcat, [tform(atom.coords) for atom in nodeattrs(mol)])
    size(pos, 1) == 3 || error("this plots only in 3d")
    if !H
        keep = syms .!= :H
        syms, pos = syms[keep], pos[:,keep]
        atomidx = findall(keep)
    else
        atomidx = collect(1:length(syms))
    end
    isvdW = radius == "van der Waals"
    radii = map(atomidx, syms) do idx, sym
        an = ATOMSYMBOLMAP[string(sym)]
        if isvdW
            r = ATOM_VANDERWAALS_RADII[an]
        else
            r = ATOM_COVALENT_RADII[an]
            if !isa(r, Real)
                # Carbon and a few metals have multiple values
                if an == 6
                    nbrs = neighbors(mol, idx)
                    key = nbrs == 3 ? "Csp3" :
                        nbrs == 2 ? "Csp2" : "Csp"
                    r = r[key]
                else
                    # For metals, it's safest to choose the smallest radius
                    r = minimum(values(r))
                end
            end
        end
        r
    end
    return axisoff!(meshscatter(pos[1,:], pos[2,:], pos[3,:], color=colortype.(atomcolor(syms)), markersize=radii, markerspace=SceneSpace))
end

"""
    ballstick(mol::UndirectGraph; tform=identity, H::Bool=true, markersize=1)

Represent `mol` as a ball-and-stick model in three dimensions. `mol` should have 3d atom positions represented
in Angstroms. 3D SDF files can be downloaded from sites such as PubChem.

`tform(xyz)` optionally transforms the positions before plotting them. `H=false` causes hydrogens to be omitted from the plot.
`markersize` specifies the diameter of the balls, in Angstroms.

This function requires that you load one of the backends of the Makie/WGLMakie/CairoMakie family.
"""
function ballstick(mol::UndirectedGraph; tform=identity, H::Bool=true, markersize=1)
    syms = [a.symbol for a in nodeattrs(mol)]
    pos = reduce(hcat, [tform(atom.coords) for atom in nodeattrs(mol)])
    size(pos, 1) == 3 || error("this plots only in 3d")
    bondinfo = map(edgesiter(mol), edgeattrs(mol)) do e, attr
        ([Vec3f0(pos[:,e[1]]), Vec3f0(pos[:,e[2]])], attr.order)
    end
    if !H
        keep = syms .!= :H
        syms, pos = syms[keep], pos[:,keep]
        bondinfo = [bi for (bi, e) in zip(bondinfo, edgesiter(mol)) if (keep[e[1]] & keep[e[2]])]
    end
    scene = scatter(pos[1,:], pos[2,:], pos[3,:], color=colortype.(atomcolor(syms)), markersize=markersize, markerspace=SceneSpace)
    for (l, w) in bondinfo
        lines!(scene, l, linewidth=10*w)
    end
    return axisoff!(scene)
end

axisoff!(scene::Scene) = (axisoff!(scene[Axis]); return scene)
function axisoff!(ax)
    # ax.showgrid[] = (false, false, false)   # disabled due to https://github.com/JuliaPlots/Makie.jl/issues/774
    ax.showaxis[] = (false, false, false)
    ax.showticks[] = (false, false, false)
    return ax
end
