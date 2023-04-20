#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    printv2mol, printv2sdf, sdfilewriter


function sdf_bond_style(bondorder, bondstyle)
    arr = zeros(Int, length(bondorder))
    for i in 1:length(bondorder)
        if bondstyle[i] in (:up, :revup)
            arr[i] = 1
        elseif bondstyle[i] in (:down, :revdown)
            arr[i] = 6
        elseif bondstyle[i] === :unspecified
            arr[i] = bondorder[i] == 2 ? 3 : 4
        end
    end
    return arr
end


function printv2atoms(io::IO, g, atomsymbol, coords)
    for i in vertices(g)
        x, y = coords[i, 1:2]
        z = 0.0  # TODO: keep 3D
        xyzsym = @sprintf "%10.4f%10.4f%10.4f %-3s" x y z string(atomsymbol[i])
        println(io, "$(xyzsym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
end


function printv2bonds(io::IO, g, bondorder, styles)
    sdfstyles = sdf_bond_style(bondorder, styles)
    for (i, e) in enumerate(edges(g))
        u, v = styles[i] in (:revup, :revdown) ? (dst(e), src(e)) : (src(e), dst(e))
        uv = @sprintf "%3d%3d%3d%3d  0  0  0" u v bondorder[i] sdfstyles[i]
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
    bondorder = bond_order(mol)
    if !hasfield(vproptype(mol), :coords) && !has_cache(mol, :v_coords2d)  # default SMILESAtom
        coords, styles = coordgen(mol)
        printv2atoms(io, mol.graph, atom_symbol(mol), coords)
        # TODO: cis-trans unspecified double bond in SMILES
        printv2bonds(io, mol.graph, bondorder,
            bond_style(bondorder, styles, 
                double_bond_style(mol.graph, bond_order(mol), zeros(ne(mol)), coords, sssr(mol))
            )
        )
    else
        printv2atoms(io, mol.graph, atom_symbol(mol), coords2d(mol))
        printv2bonds(io, mol.graph, bondorder, bond_style(bondorder, single_bond_style(mol), double_bond_style(mol)))
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
