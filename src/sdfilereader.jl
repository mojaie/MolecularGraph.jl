#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFileReader,
    sdfilereader,
    nohaltsupplier,
    sdftomol,
    RDFileReader,
    rdfilereader,
    nohalt_rxn_supplier,
    rxntoreaction

const SDF_CHARGE_TABLE = Dict(
    0 => 0, 1 => 3, 2 => 2, 3 => 1, 4 => 0, 5 => -1, 6 => -2, 7 => -3
)


function sdfatom(line)
    xpos = parse(Float64, line[1:10])
    ypos = parse(Float64, line[11:20])
    zpos = parse(Float64, line[21:30])
    coords = Float64[xpos, ypos, zpos]
    sym = Symbol(rstrip(line[32:34]))
    # atom.mass_diff = parse(Int, line[35:36]) use ISO property
    old_sdf_charge = parse(Int, line[37:39])
    charge = SDF_CHARGE_TABLE[old_sdf_charge]
    multi = old_sdf_charge == 4 ? 2 : 1
    # atom.stereo_flag = parse(Int, line[40:42])
    # valence = parse(Int, line[46:48])
    return (sym, charge, multi, coords)
end


function sdfbond(line)
    u = parse(Int, line[1:3])
    v = parse(Int, line[4:6])
    order = parse(Int, line[7:9])
    notation = parse(Int, line[10:12])
    return (u, v, order, notation)
end


function sdfprops(lines)
    props = Dict{Int,Dict{Symbol,Real}}() # atomindex, {type => value}
    for line in lines
        startswith(line, "M  ") || continue
        proptype = line[4:6]
        if !(proptype in ("CHG", "RAD", "ISO"))
            continue # Other properties are not supported yet
        end
        count = parse(Int, line[7:9])
        for c in 1:count
            i = c - 1
            idx = parse(Int, line[8i + 11 : 8i + 13])
            val = line[8i + 15 : 8i + 17]
            haskey(props, idx) || (props[idx] = Dict{Symbol,Real}())
            if proptype == "CHG"
                props[idx][:CHG] = parse(Int, val)
            elseif proptype == "RAD"
                props[idx][:RAD] = parse(Int, val)
            elseif proptype == "ISO"
                props[idx][:ISO] = parse(Float64, val)
            end
        end
    end
    return props
end


function sdfoptions(lines)
    function add_attribute!(d, key, list)
        if key != :NOTANOPTION
            while !isempty(list) && isempty(list[end])
                pop!(list)
            end
            d[key] = length(list) == 1 ? list[1] : copy(list)
            empty!(list)
        end
        return d
    end

    data = Dict{Symbol,Union{String,Vector{String}}}()
    lastoption = Ref(:NOTANOPTION)
    prevlines = String[]
    for (i, line) in enumerate(lines)
        # Some inappropriate signs are accepted for practical use
        m = match(r">.*?<([\w -.%=/]+)>", line)
        if m !== nothing
            add_attribute!(data, lastoption[], prevlines)
            lastoption[] = Symbol(m[1])
        elseif line == "\$\$\$\$"  # $$$$ indicates the next compound
            add_attribute!(data, lastoption[], prevlines)
            lastoption[] = :NOTANOPTION
            break
        else
            push!(prevlines, line)
        end
    end
    # Ensure that an unterminated file includes the final option
    add_attribute!(data, lastoption[], prevlines)

    return data
end


"""
    parse(::Type{SDFile}, lines)

Parse lines of a SDFile mol block data into a molecule object.
"""
function Base.parse(::Type{SDFile}, sdflines)
    sdflines = collect(sdflines)
    molend = findnext(x -> x == "M  END", sdflines, 1)    
    lines = @view sdflines[1:molend-1]
    optlines = @view sdflines[molend+1:end]

    # Get element blocks
    countline = lines[4]
    
    occursin("V3000", countline) && return sdftomol3000(join(sdflines, "\n"))

    atomcount = parse(UInt16, countline[1:3])
    bondcount = parse(UInt16, countline[4:6])
    # chiralflag = countline[12:15] Not used
    # propcount = countline[30:33] No longer supported
    atomoffset = 5
    bondoffset = atomcount + atomoffset
    propoffset = bondoffset + bondcount
    atomblock = @view lines[atomoffset:bondoffset-1]
    bondblock = @view lines[bondoffset:propoffset-1]
    propblock = @view lines[propoffset:end]

    # Parse atoms
    nodeattrs = SDFileAtom[]
    props = sdfprops(propblock)
    for (i, line) in enumerate(atomblock)
        (sym, charge, multi, coords) = sdfatom(line)
        mass = nothing
        if !isempty(props)
            # If prop block exists, any annotations in atom blocks are ignored
            charge = 0
            multi = 1
            if haskey(props, i)
                prop = props[i]
                haskey(prop, :CHG) && (charge = prop[:CHG])
                haskey(prop, :RAD) && (multi = prop[:RAD])
                haskey(prop, :ISO) && (mass = prop[:ISO])
            end
        end
        push!(nodeattrs, SDFileAtom(sym, charge, multi, mass, coords))
    end

    # Parse bonds
    edges = Tuple{Int,Int}[]
    edgeattrs = SDFileBond[]
    for line in bondblock
        (u, v, order, notation) = sdfbond(line)
        push!(edges, (u, v))
        push!(edgeattrs, SDFileBond(order, notation))
    end

    molobj = graphmol(edges, nodeattrs, edgeattrs)
    merge!(molobj.attributes, sdfoptions(optlines))
    return molobj
end

"""
    parse(::Type{GraphReaction}, lines)

Parse lines of a rxn file into a reaction object.
"""
function Base.parse(::Type{GraphReaction}, rxn::AbstractString)
    reactants, products = if startswith(rxn, raw"$RXN V3000")
        rr = String[m.captures[1] for m in eachmatch(blockregex("REACTANT"), rxn)]
        pp = String[m.captures[1] for m in eachmatch(blockregex("PRODUCT"), rxn)]
        sdftomol3000.(rr), sdftomol3000.(pp)
    else
        parts = split(rxn, r"\$MOL\r?\n")
        (ne, np) = parse.(Int, split(split(first(parts), r"\r?\n")[end-1]))
        rr = parts[2:1 + ne]
        pp = parts[(2 + ne):(1 + ne + np)]
        sdftomol.(IOBuffer.(rr)), sdftomol.(IOBuffer.(pp))
    end
    GraphReaction(reactants, products)
end

Base.parse(::Type{GraphReaction}, rxnlines) = parse(GraphReaction, join(rxnlines, "\n"))

function nohaltsupplier(block)
    return try
        mol = parse(SDFile, block)
        setdiastereo!(mol)
        setstereocenter!(mol)
        return mol
    catch e
        if e isa ErrorException
            println("$(e.msg) (#$(i) in sdfilereader)")
            return graphmol(SDFileAtom, SDFileBond)
        else
            throw(e)
        end
    end
end

function nohalt_rxn_supplier(block)
    return try
        reaction = parse(GraphReaction, block)
        setdiastereo!(reaction)
        setstereocenter!(reaction)
        return reaction
    catch e
        if e isa ErrorException
            println("$(e.msg) (#$(i) in rxnfilereader)")
            return Reaction()
        else
            throw(e)
        end
    end
end

struct SDFileReader
    lines::Base.EachLine
    parser::Function
end

struct RDFileReader
    io::IO
    parser::Function
end

function Base.iterate(reader::SDFileReader, state=nothing)
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

function Base.iterate(reader::RDFileReader, state=nothing)
    block = readuntil(reader.io, "\$RXN")
    while startswith(block, "\$RDFILE") && block != ""
        block = readuntil(reader.io, "\$RXN")
    end
    isempty(block) ? nothing : (reader.parser("\$RXN" * block), state)
end

Base.IteratorSize(::Type{SDFileReader}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{SDFileReader}) = Base.EltypeUnknown()

Base.IteratorSize(::Type{RDFileReader}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{RDFileReader}) = Base.EltypeUnknown()

"""
    sdfilereader(file::IO)
    sdfilereader(path::AbstractString)

Read SDFile data from input stream (or a file path as a string) and return a
lazy iterator that yields molecule objects.

`sdfilereader` does not stop and raise errors when an erroneous or incompatible
SDFile block is read but produces an error message and yields an empty molecule.
If this behavior is not desirable, you can use the customized supplier function
instead of default supplier `nohaltsupplier`

```
function customsupplier()
    mol = try
        parse(SDFile, block)
    catch e
        throw(ErrorException("incompatible molecule found, aborting..."))
    end
end

function sdfilereader(file::IO)
    return SDFileReader(eachline(file), customsupplier)
end
```
"""
sdfilereader(file::IO) = SDFileReader(eachline(file), nohaltsupplier)
sdfilereader(path::AbstractString) = sdfilereader(open(path))

rdfilereader(file::IO) = RDFileReader(file, nohalt_rxn_supplier)
rdfilereader(path::AbstractString) = rdfilereader(open(path))

"""
    sdftomol(lines) -> GraphMol{SDFileAtom,SDFileBond}
    sdftomol(file::IO) -> GraphMol{SDFileAtom,SDFileBond}
    sdftomol(path::AbstractString) -> GraphMol{SDFileAtom,SDFileBond}

Read a SDFile(.sdf or .mol) and parse it into a molecule object. The given
argument should be a file input stream, a file path as a string or an iterator
that yields each sdfile text lines.
"""
function sdftomol(lines)
    mol = parse(SDFile, lines)
    setdiastereo!(mol)
    setstereocenter!(mol)
    return mol
end
sdftomol(file::IO) = sdftomol(eachline(file))
sdftomol(path::AbstractString) = sdftomol(open(path))


"""
    rxntoreaction(lines) -> GraphReaction
    rxntoreaction(file::IO) -> GraphReaction
    sdftomol(path::AbstractString) -> GraphReaction

Read a RXN file and parse it into a reaction object. The given
argument should be a file input stream, a file path as a string or an iterator
that yields each sdfile text lines.
"""
function rxntoreaction(lines)
    reaction = parse(GraphReaction, lines)
    setdiastereo!(reaction)
    setstereocenter!(reaction)
    return reaction
end
rxntoreaction(file::IO) = rxntoreaction(eachline(file))
rxntoreaction(path::AbstractString) = rxntoreaction(open(path))

# support of extended mol file format (V3000)

sympair(s) = (p -> Symbol(p[1]) => p[2])(split(s,"="))
blockregex(s::AbstractString) = Regex("M  V30 BEGIN $(s)\r?\n(.*?)\r?\nM  V30 END $(s)", "s")

function parseatomblock3000(atomblock)
    nodeattrs = SDFileAtom[]
    for line in eachline(IOBuffer(atomblock))
        ss = split(line)
        coords = parse.(Float64, ss[5:7])
        sym = Symbol(ss[4])
        
        props = Dict(sympair.(ss[9:end])...)
        
        charge = parse(Int, get(props, :CHG, "0"))
        mass = tryparse(Int, get(props, :MASS, ""))

        push!(nodeattrs, SDFileAtom(sym, charge, 1, mass, coords))
    end
    nodeattrs
end

function parsebondblock3000(bondblock)
    edges = Tuple{Int,Int}[]
    edgeattrs = SDFileBond[]

    for line in eachline(IOBuffer(bondblock))
        ss = split(line)
        order, u, v = parse.(Int, ss[4:6])
        props = Dict(sympair.(ss[7:end])...)
        notation = get(props, :CFG, 0)
        push!(edges, (u, v))
        push!(edgeattrs, SDFileBond(order, notation))
    end
    edges, edgeattrs
end

function sdftomol3000(s::AbstractString)
    atomblock = match(blockregex("ATOM"), s).captures[1]
    bondblock = match(blockregex("BOND"), s).captures[1]

    nodeattrs = parseatomblock3000(atomblock)
    edges, edgeattrs = parsebondblock3000(bondblock)

    mol = graphmol(edges, nodeattrs, edgeattrs)
    # setdiastereo!(mol)
    # setstereocenter!(mol)
end
