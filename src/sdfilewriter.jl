#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    sdfilewriter


function atomblock(io::IO, mol::GraphMol)
    for (i, atom) in enumerate(nodeattrs(mol))
        (x, y, z) = fmt.("10.4f", atom.coords)
        sym = fmt("-3s", string(atom.symbol))
        println(io, "$(x)$(y)$(z) $(sym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
    return
end


function bondblock(io::IO, mol::GraphMol)
    for (i, (u, v)) in enumerate(edgesiter(mol))
        ustr = fmt("3d", u)
        vstr = fmt("3d", v)
        bond = edgeattr(mol, i)
        order = fmt("3d", bond.order)
        # TODO: SmilesBond
        stereo = fmt("3d", bond.notation)
        println(io, "$(ustr)$(vstr)$(order)$(stereo)  0  0  0")
    end
    return
end


function propertyblock(io::IO, mol::GraphMol)
    charges = Tuple{Int,Int}[]
    radicals = Tuple{Int,Int}[]
    masses = Tuple{Int,Float64}[]
    for (i, atom) in enumerate(nodeattrs(mol))
        atom.charge == 0 || push!(charges, (i, atom.charge))
        atom.multiplicity == 1 || push!(radicals, (i, atom.multiplicity))
        atom.mass === nothing || push!(masses, (i, atom.mass))
    end
    if !isempty(charges)
        rcds = (format(" {:3d} {:3d}", i, chg) for (i, chg) in charges)
        println(io, format("M  CHG{:3d}{}", length(charges), join(rcds, "")))
    end
    if !isempty(radicals)
        rcds = (format(" {:3d} {:3d}", i, rad) for (i, rad) in radicals)
        println(io, format("M  RAD{:3d}{}", length(radicals), join(rcds, "")))
    end
    if !isempty(masses)
        rcds = (format(" {:3d} {:3d}", i, iso) for (i, iso) in masses)
        println(io, format("M  ISO{:3d}{}", length(masses), join(rcds, "")))
    end
    return
end


function datablock(io::IO, mol::GraphMol)
    for (key, val) in mol.attributes
        println(io, "> <$(string(key))>")
        println(io, string(val))
        println(io, "")
    end
    println(io, raw"$$$$")
    return
end


function molblock(io::IO, mol::GraphMol)
    # TODO: 2D coords for SmilesMol
    println(io)
    println(io, "MolecularGraph.jl version $(version())")
    println(io)
    chiral_flag = 0 # TODO: deplicated?
    println(io, format(
        "{:3d}{:3d}  0  0{:3d}  0  0  0  0  0999 V2000",
        nodecount(mol), edgecount(mol), chiral_flag
    ))
    atomblock(io, mol)
    bondblock(io, mol)
    propertyblock(io, mol)
    println(io, "M  END")
    datablock(io, mol)
    return
end


"""
    sdfilewriter(io::IO, mols)
    sdfilewriter(filename::AbstractString, mols)

Write molecule data to the output stream as a SDFile format file.
"""
sdfilewriter(io::IO, mols) = molblock.((io,), mols)
sdfilewriter(filename::AbstractString, mols
    ) = sdfilewriter(open(path, "w"), mols)
