#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    printv2mol, printv2sdf, sdfilewriter


function printv2atoms(io::IO, mol::SimpleMolGraph, coords)
    for i in vertices(mol)
        x, y = coords[i, 1:2]
        z = 0.0  # TODO: keep 3D
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(get_prop(mol, i, :symbol))
        println(io, "$(xyzsym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
end


function printv2bonds(io::IO, mol::SimpleMolGraph, styles)
    for (i, e) in enumerate(edges(mol))
        u, v = styles[i] in (:revup, :revdown) ? (dst(e), src(e)) : (src(e), dst(e))
        # TODO: double bond, unknown
        direction = styles[i] in (:up, :revup) ? 1 : (styles in (:down, :revdown) ? 6 : 0)
        uv = @sprintf "%3d%3d%3d%3d  0  0  0" u v get_prop(mol, e, :order) direction
        println(io, uv)
    end
end


function printv2properties(io::IO, mol::SimpleMolGraph)
    charges = Tuple{Int,Int}[]
    radicals = Tuple{Int,Int}[]
    masses = Tuple{Int,Float64}[]
    for i in vertices(mol)
        get_prop(mol, i, :charge) == 0 || push!(charges, (i, get_prop(mol, i, :charge)))
        get_prop(mol, i, :multiplicity) == 1 || push!(radicals, (i, get_prop(mol, i, :multiplicity)))
        get_prop(mol, i, :mass) === nothing || push!(masses, (i, get_prop(mol, i, :mass)))
    end
    if !isempty(charges)
        head = @sprintf "M  CHG%3d" length(charges)
        print(io, head)
        for (i, chg) in charges
            rcd = @sprintf " %3d %3d" i chg
            print(io, rcd)
        end
        println(io)
    end
    if !isempty(radicals)
        head = @sprintf "M  RAD%3d" length(radicals)
        print(io, head)
        for (i, rad) in radicals
            rcd = @sprintf " %3d %3d" i rad
            print(io, rcd)
        end
        println(io)
    end
    if !isempty(masses)
        head = @sprintf "M  ISO%3d" length(masses)
        print(io, head)
        for (i, iso) in masses
            rcd = @sprintf " %3d %3d" i iso
            print(io, rcd)
        end
        println(io)
    end
end


function printv2data(io::IO, mol::SimpleMolGraph)
    for (key, val) in props(mol)
        println(io, "> <$(string(key))>")
        println(io, string(val))
        println(io, "")
    end
end


function printv2mol(io::IO, mol::SimpleMolGraph)
    # ver = VERSION
    program = "MGjlv$(string(MAJOR_VERSION)[end])$(string(MINOR_VERSION)[end-1:end])"
    datetime = Dates.format(Dates.now(), "mmddyyHHMM")
    println(io)
    println(io, "  $(program)$(datetime)2D            ")
    println(io)
    ncnt = nv(mol)
    ecnt = ne(mol)
    header = @sprintf "%3d%3d  0  0  0  0  0  0  0  0999 V2000" ncnt ecnt
    println(io, header)
    if !hasfield(vproptype(mol), :coords) && !has_state(mol, :v_coords2d)  # default SMILESAtom
        coords, styles = coordgen(mol)
        printv2atoms(io, mol, coords)
        printv2bonds(io, mol, styles)
    else
        printv2atoms(io, mol, coords2d(mol))
        printv2bonds(io, mol, single_bond_style(mol))
    end
    printv2properties(io, mol)
    println(io, "M  END")
end


function printv2mol(mol::SimpleMolGraph)
    buf = IOBuffer(write=true)
    printv2mol(buf, mol)
    res = String(take!(buf))
    close(buf)
    return res
end


function printv2sdf(io::IO, mol::SimpleMolGraph)
    printv2mol(io, mol)
    printv2data(io, mol)
    println(io, raw"$$$$")
end


"""
    sdfilewriter(io::IO, mols)
    sdfilewriter(filename::AbstractString, mols)

Write molecule data to the output stream as a SDFile format file.
"""
function sdfilewriter(io::IO, mols; writer=printv2sdf)
    cnt = length(writer.((io,), mols))
    @info "$(cnt) records exported."
end

function sdfilewriter(filename::AbstractString, mols; kwargs...)
    open(filename, "w") do io
        sdfilewriter(io, mols; kwargs...)
    end
end
