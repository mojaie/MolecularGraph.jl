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

    Channel(ctype=MolecularGraph, csize=0) do channel::Channel{MolecularGraph}
        for (i, (mol, opt)) in enumerate(sdfblock)
            molobj = try
                m = parsemol(mol)
                if precalc
                    assign_descriptors!(m)
                end
                m
            catch e
                if !nohalt
                    # TODO: stacktrace
                    throw(e)
                elseif isa(e, DescriptorError)
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
            molobj.data = parseoptions(opt)
            put!(channel, molobj)
        end
    end
end


function parsemol(lines::AbstractArray{String})
    countline = lines[4]
    atomcount = parse(UInt16, countline[1:3])
    bondcount = parse(UInt16, countline[4:6])
    # chiralflag = countline[12:15] Not used
    # propcount = countline[30:33] No longer supported
    mol = MolecularGraph()
    atomblock = @view lines[ 5 : atomcount + 4 ]
    for atom in parseatoms(atomblock)
        newatom!(mol, atom)
    end
    bondblock = @view lines[ atomcount + 5 : atomcount + bondcount + 4 ]
    for bond in parsebonds(bondblock)
        updatebond!(mol, bond)
    end
    propblock = @view lines[ atomcount + bondcount + 5 : end ]
    props = parseprops(propblock)
    if length(props) > 0
        # props supersedes all charge and radical values in the atom block
        for atom in atomvector(mol)
            atom.charge = 0
            atom.multiplicity = 1
            atom.mass = nothing
        end
    end
    for (i, ptype, val) in props
        if ptype == "CHG"
            getatom(mol, i).charge = val
        elseif ptype == "RAD"
            getatom(mol, i).multiplicity = val
        elseif ptype == "ISO"
            getatom(mol, i).mass = val
        end
    end
    mol
end


function parseatoms(lines::AbstractArray{String})
    conv_charge_table = Dict([
        (0, 0), (1, 3), (2, 2), (3, 1), (4, 0), (5, -1), (6, -2), (7, -3)
    ])
    results = []
    for (i, line) in enumerate(lines)
        sym = line[32:34]
        symbol = rstrip(sym)
        atom = Atom(symbol)
        atom.index = i
        xpos = parse(Float32, line[1:10])
        ypos = parse(Float32, line[11:20])
        zpos = parse(Float32, line[21:30])
        atom.coords = [xpos, ypos, zpos]
        # atom.mass_diff = parse(Int, line[35:37]) use ISO property
        old_sdf_charge = parse(Int8, line[38:40])
        atom.charge = conv_charge_table[old_sdf_charge]
        if old_sdf_charge == 4
            atom.multiplicity = 2
        end
        # atom.stereo_flag = parse(Int, line[41:43])
        # valence = parse(Int, line[47:49])
        push!(results, atom)
    end
    results
end


function parsebonds(lines::AbstractArray{String})
    conv_stereo_table = Dict([
        (0, 0), (1, 1), (3, 3), (4, 5), (6, 3)
    ])
    results = []
    for line in lines
        bond = Bond()
        first = parse(UInt16, line[1:3])
        second = parse(UInt16, line[4:6])
        ordered = first < second
        bond.u = ordered ? first : second
        bond.v = ordered ? second : first
        bond.order = parse(UInt8, line[7:9])
        bond.notation = conv_stereo_table[parse(UInt8, line[10:12])]
        if !ordered && bond.order == 1
            if bond.notation == 1
                bond.notation = 2
            elseif bond.notation == 3
                bond.notation = 4
            end
        end
        push!(results, bond)
    end
    results
end


function parseprops(lines::AbstractArray{String})
    results = []
    for line in lines
        proptype = line[4:6]
        if proptype ∉ ("CHG", "RAD", "ISO")
            continue # Other properties are not supported yet
        end
        count = parse(UInt8, line[7:9])
        for i in 0 : count - 1
            idx = parse(UInt16, line[8i + 11 : 8i + 13])
            val = parse(Int16, line[8i + 15 : 8i + 17])
            push!(results, (idx, proptype, val))
        end
    end
    results
end


function parseoptions(lines::AbstractArray{String})
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
