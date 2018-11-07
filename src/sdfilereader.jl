#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    sdftomol,
    sdfilereader,
    sdfbatchsupplier,
    nohaltsupplier,
    parsesdf,
    parsesdfblocks,
    parsesdfmol,
    parsesdfatom,
    parsesdfbond,
    sdfatomprops!,
    parsesdfoptions


function sdftomol(lines)
    mol = parsesdf(lines, mutable=false)
    default_annotation!(mol)
    mol
end

sdftomol(file::IO) = sdftomol(eachline(file))


function sdfilereader(lines)
    sdfbatchsupplier(nohaltsupplier, sdftomol, lines)
end

sdfilereader(file::IO) = sdfilereader(eachline(file))


function sdfbatchsupplier(supplier, sdfparser, lines)
    # TODO: Return value of sdfparser should be Molecule (not MutableMolecule)
    Channel() do channel
        for (i, block) in enumerate(parsesdfblocks(lines))
            put!(channel, supplier(sdfparser, block, i))
        end
    end
end


function nohaltsupplier(sdfparser, block, i)
    try
        sdfparser(block)
    catch e
        if e isa GraphMolError
            println("$(e.msg) (#$(i) in sdfilereader)")
            sdfparser(nullmol())
        else
            throw(e)
        end
    end
end


function parsesdf(lines; mutable=false)
    lines = collect(lines)
    molend = findnext(x -> x == "M  END", lines, 1)
    mollines = view(lines, 1:(molend - 1))
    optlines = view(lines, (molend + 1):lastindex(lines))
    molobj = parsesdfmol(mollines, mutable)
    molobj.attribute[:sourcetype] = :sdfile
    merge!(molobj.attribute, parsesdfoptions(optlines))
    molobj
end


function parsesdfblocks(lines)
    Channel() do channel
        sdfblock = []
        for line in lines
            if startswith(line, raw"$$$$")
                put!(channel, sdfblock)
                empty!(sdfblock)
            else
                push!(sdfblock, rstrip(line))
            end
        end
        if !isempty(sdfblock)
            put!(channel, sdfblock)
        end
    end
end


function parsesdfmol(lines, mutable)
    countline = lines[4]
    atomcount = parse(UInt16, countline[1:3])
    bondcount = parse(UInt16, countline[4:6])
    # chiralflag = countline[12:15] Not used
    # propcount = countline[30:33] No longer supported
    atomblock = @view lines[ 5 : atomcount + 4 ]
    atoms = parsesdfatom.(atomblock)
    bondblock = @view lines[ atomcount + 5 : atomcount + bondcount + 4 ]
    bonds = parsesdfbond.(bondblock)
    propblock = @view lines[ atomcount + bondcount + 5 : end ]
    sdfatomprops!(atoms, propblock)
    aobjs = [sdfatom(atom...) for atom in atoms]
    bobjs = [sdfbond(bond...) for bond in bonds]
    mutable ? MutableMolecule(aobjs, bobjs) : Molecule(aobjs, bobjs)
end


const SDF_CHARGE_TABLE = Dict(
    0 => 0, 1 => 3, 2 => 2, 3 => 1, 4 => 0, 5 => -1, 6 => -2, 7 => -3
)


function parsesdfatom(line)
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


const SDF_STEREO_TABLE = Dict(0 => 0, 1 => 1, 3 => 3, 4 => 5, 6 => 3)


function parsesdfbond(line)
    """Bond

    * Notation
        * Single bond
            * 0: u - v
            * 1: u ◀ v (Up-arrow)
            * 2: u ▶ v
            * 3: u ◁ v (Down-arrow)
            * 4: u ▷ v
            * 5: u ~ v (Chiral)
        * Double bond
            * 0: v ニ u (clockwise, default)
            * 1: u ニ v (counter-clockwise)
            * 2: u ＝ v (equal length, for terminal bond by default)
            * 3: u × v (Cis-Trans Unknown)
    """
    u = parse(Int, line[1:3])
    v = parse(Int, line[4:6])
    """ TODO: no need to order?
    first = parse(UInt16, line[1:3])
    second = parse(UInt16, line[4:6])
    ordered = first < second
    u = ordered ? first : second
    v = ordered ? second : first
    """
    order = parse(Int, line[7:9])
    notation = parse(Int, line[10:12])
    "
    notation = SDF_STEREO_TABLE[parse(UInt8, line[10:12])]
    if !ordered && order == 1
        if notation == 1
            notation = 2
        elseif notation == 3
            notation = 4
        end
    end
    "
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


function parsesdfoptions(lines)
    data = Dict()
    for (i, line) in enumerate(lines)
        # Some inappropriate signs are accepted for practical use
        m = match(r">.*?<([\w -.%=]+)>", line)
        if m !== nothing
            data[m[1]] = lines[i + 1]
        end
    end
    data
end
