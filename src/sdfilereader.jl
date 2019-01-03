#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFileReader,
    sdfilereader,
    nohaltsupplier,
    defaultpostprocess,
    sdftomol,
    sdfmol,
    sdfatom,
    sdfbond,
    sdfatomprops!,
    sdfoptions


import Base: parse, iterate, eltype, IteratorSize, IteratorEltype


struct SDFileReader
    lines::Base.EachLine
    parser::Function
end

sdfilereader(file::IO) = SDFileReader(eachline(file), nohaltsupplier)

function iterate(reader::SDFileReader, state=nothing)
    block = String[]
    next = iterate(reader.lines)
    while next !== nothing
        (line, state) = next
        if startswith(line, raw"$$$$")
            return (reader.parser(block), state)
        end
        push!(block, rstrip(line))
        next = iterate(reader.lines, state)
    end
    if !isempty(block)
        return (reader.parser(block), state)
    end
    return
end

eltype(::Type{SDFileReader}) = Vector{String}
IteratorSize(::Type{SDFileReader}) = Base.SizeUnknown()
IteratorEltype(::Type{SDFileReader}) = Base.HasEltype()


function nohaltsupplier(block)
    mol = try
        parse(SDFile, block)
    catch e
        if e isa MolParseError
            println("$(e.msg) (#$(i) in sdfilereader)")
            nullmol(SDFile)
        else
            throw(e)
        end
    end
    return defaultpostprocess(mol)
end


function defaultpostprocess(mol::SDFile)
    vmol = vectormol(mol)
    default_annotation!(vmol)
    return vmol
end


function sdftomol(lines)
    mol = parse(SDFile, lines)
    return defaultpostprocess(mol)
end

sdftomol(file::IO) = sdftomol(eachline(file))


# TODO: Type{SDFile}
function parse(::Type{T}, lines) where T <: SDFile
    lines = collect(lines)
    molend = findnext(x -> x == "M  END", lines, 1)
    mollines = view(lines, 1:(molend - 1))
    optlines = view(lines, (molend + 1):lastindex(lines))
    molobj = sdfmol(mollines)
    merge!(molobj.attribute, sdfoptions(optlines))
    return molobj
end


function sdfmol(lines)
    countline = lines[4]
    atomcount = parse(UInt16, countline[1:3])
    bondcount = parse(UInt16, countline[4:6])
    # chiralflag = countline[12:15] Not used
    # propcount = countline[30:33] No longer supported
    atomblock = @view lines[ 5 : atomcount + 4 ]
    atoms = sdfatom.(atomblock)
    bondblock = @view lines[ atomcount + 5 : atomcount + bondcount + 4 ]
    bonds = sdfbond.(bondblock)
    propblock = @view lines[ atomcount + bondcount + 5 : end ]
    sdfatomprops!(atoms, propblock)
    aobjs = [SDFileAtom(atom...) for atom in atoms]
    bobjs = [SDFileBond(bond...) for bond in bonds]
    return SDFile(aobjs, bobjs)
end


const SDF_CHARGE_TABLE = Dict(
    0 => 0, 1 => 3, 2 => 2, 3 => 1, 4 => 0, 5 => -1, 6 => -2, 7 => -3
)


function sdfatom(line)
    sym = Symbol(rstrip(line[32:34]))
    xpos = parse(Float64, line[1:10])
    ypos = parse(Float64, line[11:20])
    zpos = parse(Float64, line[21:30])
    coords = SVector(xpos, ypos, zpos)
    # atom.mass_diff = parse(Int, line[35:36]) use ISO property
    old_sdf_charge = parse(Int, line[37:39])
    charge = SDF_CHARGE_TABLE[old_sdf_charge]
    multi = old_sdf_charge == 4 ? 2 : 1
    # atom.stereo_flag = parse(Int, line[40:42])
    # valence = parse(Int, line[46:48])
    [sym, charge, multi, nothing, coords]
end


function sdfbond(line)
    u = parse(Int, line[1:3])
    v = parse(Int, line[4:6])
    order = parse(Int, line[7:9])
    notation = parse(Int, line[10:12])
    [u, v, order, notation]
end


function sdfatomprops!(atoms, lines)
    if isempty(lines)
        return
    end
    # props supersedes all charge and radical values in the atom block
    for atom in atoms
        atom[2] = 0
        atom[3] = 1
    end
    results = []
    for line in lines
        proptype = line[4:6]
        if !(proptype in ("CHG", "RAD", "ISO"))
            continue # Other properties are not supported yet
        end
        count = parse(Int, line[7:9])
        for i in 0 : count - 1
            idx = parse(Int, line[8i + 11 : 8i + 13])
            val = line[8i + 15 : 8i + 17]
            if proptype == "CHG"
                atoms[idx][2] = parse(Int, val)
            elseif proptype == "RAD"
                atoms[idx][3] = parse(Int, val)
            elseif proptype == "ISO"
                atoms[idx][4] = parse(Float64, val)
            end
        end
    end
    return
end


function sdfoptions(lines)
    data = Dict()
    for (i, line) in enumerate(lines)
        # Some inappropriate signs are accepted for practical use
        m = match(r">.*?<([\w -.%=]+)>", line)
        if m !== nothing
            data[m[1]] = lines[i + 1]
        end
    end
    return data
end
