#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export loadsdfiter, loadsdfmol


function loadsdfiter(file::IO; nohalt=true, precalc=true)
    loadsdfiter(eachline(file), nohalt=false, precalc=precalc)
end


function loadsdfiter(data; nohalt=true, precalc=true)
    parseblock(data, false, precalc)
end


function loadsdfmol(file::IO; precalc=true)
    loadsdfmol(eachline(file), precalc=precalc)
end


function loadsdfmol(data; precalc=true)
    moliter = loadsdfiter(data, nohalt=false, precalc=precalc)
    iterate(moliter)[1]
end


function parseblock(lines, nohalt, precalc)
    sdfblock = Channel(ctype=Tuple, csize=0) do channel::Channel{Tuple}
        mol = String[]
        opt = String[]
        ismol = true
        for line in lines
            if startswith(line, raw"$$$$")
                put!(channel, (copy(mol), copy(opt)))
                ismol = true
                empty!(mol)
                empty!(opt)
            elseif startswith(line, "M  END")
                ismol = false
            elseif ismol
                push!(mol, rstrip(line))
            else
                push!(opt, rstrip(line))
            end
        end
        if !isempty(mol)
            put!(channel, (mol, opt))
        end
    end

    Channel(ctype=Molecule, csize=0) do channel::Channel{Molecule}
        for (i, (mol, opt)) in enumerate(sdfblock)
            molobj = try
                m = parsesdfmol(mol)
                if precalc
                    default_annotation!(m)
                end
                m
            catch e
                if !nohalt
                    # TODO: stacktrace
                    throw(e)
                elseif isa(e, AnnotationError)
                    print("Unsupported symbol: $(e) (#$(i+1) in sdfilereader)")
                    m = nullmol(precalc)
                elseif isa(e, OperationError)
                    print("Failed to minimize ring: $(e) (#$(i+1) in sdfilereader)")
                else
                    print("Unexpected error: (#$(i+1) in sdfilereader)")
                    m = nullmol(precalc)
                end
                m
            end
            molobj.attribute[:sourcetype] = :sdfile
            merge!(molobj.attribute, parsesdfoptions(opt))
            put!(channel, molobj)
        end
    end
end


function parsesdfmol(lines)
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
    Molecule(aobjs, bobjs)
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
            val = parse(Int, line[8i + 15 : 8i + 17])
            if proptype == "CHG"
                atoms[idx][2] = val
            elseif proptype == "RAD"
                atoms[idx][3] = val
            elseif proptype == "ISO"
                atoms[idx][4] = val
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
