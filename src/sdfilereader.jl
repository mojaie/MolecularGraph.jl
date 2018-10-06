#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export loadsdfiter, loadsdfmol


function loadsdfiter(file::IO, no_halt=true, precalc=true)
    loadsdfiter(eachline(file), false, precalc)
end


function loadsdfiter(data, no_halt=true, precalc=true)
    parseblock(data, false, precalc)
end


function loadsdfmol(file::IO, precalc=true)
    loadsdfmol(eachline(file), precalc)
end


function loadsdfmol(data, precalc=true)
    moliter = loadsdfiter(data, false, precalc)
    iterate(moliter)[0]
end


function parseblock(lines, nohalt, precalc)
    sdfblock = Channel(ctype=Tuple, csize=0) do channel::Channel{Tuple}
        mol = []
        opt = []
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
        if mol
            put!(channel, (mol, opt))
        end
    end

    Channel(ctype=MolecularGraph, csize=0) do channel::Channel{MolecularGraph}
        for (i, (mol, opt)) in enumerate(sdfblock)
            try
                c = parsemol(mol)
                if precalc
                    assign_descriptors!(c)
                end
            catch e
                if isa(e, ErrorException)
                    if nohalt
                        print("Unsupported symbol: $(e) (#$(i+1) in sdfilereader)")
                        c = nullmol(precalc)
                    else
                        throw(ErrorException(e, "Unsupported symbol: $(e)"))
                    end
                elseif isa(e, ErrorException)
                    if nohalt
                        print("Failed to minimize ring: $(e) (#$(i+1) in sdfilereader)")
                    else
                        throw(ErrorException(e, "Failed to minimize ring: $(e)"))
                    end
                else
                    if nohalt
                        print("Unexpected error: (#$(i+1) in sdfilereader)")
                        c = nullmol(precalc)
                        c.data = parseoption(opt)
                        put!(channel, c)
                        continue
                    else
                        # stacktrace
                        error("Unsupported Error")
                    end
                end
                c.data = parseoption(opt)
                put!(channel, c)
            end
        end
    end
end


function parsemol(lines::AbstractArray{String})
    countline = lines[4]
    atomcount = countline[1:3]
    bondcount = countline[4:6]
    # chiralflag = countline[12:15] Not used
    # propcount = countline[30:33] No longer supported
    mol = MolecularGraph()
    atomblock = @view lines[ 5 : atomcount + 4 ]
    mol.graph.nodes = parseatom(atomblock)
    bondblock = @view lines[ atomcount + 5 : atomcount + bondcount + 4 ]
    mol.graph.adjacency = parsebond(bondblock, mol.graph.nodes)
    propblock = @view lines[ atomcount + bondcount + 5 : end ]
    props = parseprop(propblock)
    addprops(props, mol)
    mol
end


function parseatom(lines::AbstractArray{String})
    conv_charge_table = Dict([
        (0, 0), (1, 3), (2, 2), (3, 1), (4, 0), (5, -1), (6, -2), (7, -3)
    ])
    results = Dict()
    for (i, line) in enumerate(lines)
        sym = line[32:34]
        symbol = rstrip(sym)
        atom = Atom(symbol)
        print(atom)
        atom = try
            Atom(symbol)
        catch e
            if isa(e, KeyError)
                throw(ErrorException(e, symbol))
            end
        end
        xpos = parse(Float64, line[1:10])
        ypos = parse(Float64, line[11:20])
        zpos = parse(Float64, line[21:30])
        atom.coords = (xpos, ypos, zpos)
        # atom.mass_diff = parse(Int, line[35:37]) use ISO property
        old_sdf_charge = parse(Int, line[38:40])
        atom.charge = conv_charge_table[old_sdf_charge]
        if old_sdf_charge == 4
            atom.multiplicity = 2
        end
        # atom.stereo_flag = parse(Int, line[41:43])
        # valence = parse(Int, line[47:49])
        results[i] = Dict([("atom", atom)])
    end
    results
end


function parsebond(lines::AbstractArray{String}, atoms::AbstractArray)
    conv_stereo_table = Dict([
        (0, 0), (1, 1), (3, 3), (4, 3), (6, 2)
    ])
    results = Dict([a => Dict() for a in atoms])

end
