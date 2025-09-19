#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const SDF_CHARGE_TABLE = Dict(
    0 => 0, 1 => 3, 2 => 2, 3 => 1, 4 => 0, 5 => -1, 6 => -2, 7 => -3
)


function sympair(s)::Pair{Symbol, Union{String, Int}}
    p = split(s,"=")
    length(p) != 2 && return :nothing => ""
    int = tryparse(Int, p[2])
    Symbol(p[1]) => isnothing(int) ? p[2] : int
end


function ctab_atom_v2(::Type{T}, line::AbstractString) where T <: AbstractElement
    d = Dict{String,Any}()
    xpos = parse(Float64, line[1:10])
    ypos = parse(Float64, line[11:20])
    zpos = parse(Float64, line[21:30])
    d["coords"] = Float64[xpos, ypos, zpos]
    d["symbol"] = rstrip(line[32:34])
    # atom.mass_diff = parse(Int, line[35:36])
    d["isotope"] = 0  # will be ignored, use ISO property instead
    old_sdf_charge = parse(Int, line[37:39])
    d["charge"] = SDF_CHARGE_TABLE[old_sdf_charge]
    d["multiplicity"] = old_sdf_charge == 4 ? 2 : 1
    # atom.stereo_flag = parse(Int, line[40:42])
    # valence = parse(Int, line[46:48])
    return T(d)
end

function ctab_atom_v3(::Type{T}, line::AbstractString) where T <: AbstractElement
    d = Dict{String,Any}()
    ss = split(line)
    d["coords"] = parse.(Float64, ss[5:7])
    d["symbol"] = ss[4]
    props = Dict(sympair.(ss[9:end])...)
    d["charge"] = get(props, :CHG, 0)
    d["isotope"] = get(props, :MASS, 0)
    d["multiplicity"] = get(props, :RAD, 1)
    return T(d)
end

function ctab_bond_v2(
        ::Type{T}, ::Type{E}, line::AbstractString) where {T<:Integer,E<:AbstractElement}
    d = Dict{String,Any}()
    u = parse(T, line[1:3])
    v = parse(T, line[4:6])
    d["order"] = parse(Int, line[7:9])
    d["notation"] = parse(Int, line[10:12])
    d["isordered"] = u < v
    u, v = d["isordered"] ? (u, v) : (v, u)
    return (u_edge(T, u, v), E(d))
end

function ctab_bond_v3(
        ::Type{T}, ::Type{E}, line::AbstractString) where {T<:Integer,E<:AbstractElement}
    d = Dict{String,Any}()
    ss = split(line)
    d["order"], u, v = parse.(T, ss[4:6])
    props = Dict(sympair.(ss[7:end])...)
    # d["notation"] = get(props, :CFG, 0)  # TODO: not compatible with v2
    d["isordered"] = u < v
    u, v = d["isordered"] ? (u, v) : (v, u)
    return (u_edge(T, u, v), E(d))
end

function ctab_props_v2(io::IO)
    props = Dict{Int,Dict{Symbol,Real}}() # atomindex, {type => value}
    while true
        line = readline(io)
        line == "M  END" && break
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
                props[idx][:ISO] = parse(Int, val)
            end
        end
    end
    return props
end


function sdf_on_init!(mol::SimpleMolGraph)
    coords_from_sdf!(mol)
    stereocenter_from_sdf2d!(mol)
    stereobond_from_sdf2d!(mol)
end

function sdf_on_update!(mol::SimpleMolGraph)
    # Preprocess
    default_atom_charge!(mol)
    default_bond_order!(mol)
    # Cache relatively expensive descriptors
    sssr!(mol)
    apparent_valence!(mol)
    valence!(mol)
    lone_pair!(mol)
    is_ring_aromatic!(mol)
end


function parse_ctab(::Type{T}, io::IO, config::Dict{Symbol,Any}) where T <: SimpleMolGraph
    line1 = readline(io)  # name line, not implemented
    if startswith(line1, "M  V30")  # v3 no header (rxnfile)
        startswith(line1, "M  V30 BEGIN CTAB") || (readuntil(io, "M  V30 BEGIN CTAB"); readline(io))  # skip END tags
        ctab_only = true
        ver = :v3
        ctab_atom = ctab_atom_v3
        ctab_bond = ctab_bond_v3
        count_line = split(readline(io))
        atomcount = parse(Int, count_line[4])
        bondcount = parse(Int, count_line[5])
    else
        ctab_only = false
        line2 = readline(io)  # format properties, not implemented
        # program_name = line2[3:10]
        # datetime = line2[11:20]
        # dimension = line2[21:22]
        line3 = readline(io)  # comment line, not implemented
        line4 = readline(io)
        if endswith(line4, "V2000")
            ver = :v2
            ctab_atom = ctab_atom_v2
            ctab_bond = ctab_bond_v2
            atomcount = parse(Int, line4[1:3])
            bondcount = parse(Int, line4[4:6])
            # chiral = countline[12:15]  # chiralflag, not implemented
        elseif endswith(line4, "V3000")
            ver = :v3
            ctab_atom = ctab_atom_v3
            ctab_bond = ctab_bond_v3
            readline(io)  # M  V30 BEGIN CTAB
            count_line = split(readline(io))
            atomcount = parse(Int, count_line[4])
            bondcount = parse(Int, count_line[5])
            # nsg = parse(Int, count_line[7])  # Sgroup, not implemented
            # n3d = parse(Int, count_line[8])  # 3D constraint, not implemented
            # chiral = parse(Int, count_line[9])  # chiralflag, not implemented
        else
            error("unsupported sdfile format - unknown ctab version")
        end
    end

    V = vproptype(T) === AbstractAtom ? SDFAtom : vproptype(T)
    E = eproptype(T) === AbstractBond ? SDFBond : eproptype(T)
    # Parse atoms
    ver === :v3 && (readuntil(io, "M  V30 BEGIN ATOM"); readline(io))
    vprops = V[]
    for _ in 1:atomcount
        push!(vprops, ctab_atom(V, readline(io)))
    end

    # Parse bonds
    ver === :v3 && (readuntil(io, "M  V30 BEGIN BOND"); readline(io))
    edges = edgetype(T)[]
    eprops = E[]
    for _ in 1:bondcount
        edge, eprop = ctab_bond(eltype(T), E, readline(io))
        push!(edges, edge)
        push!(eprops, eprop)
    end

    # Properties
    if ver === :v2
        for (i, ps) in ctab_props_v2(io)
            d = Dict{String,Any}()
            d["symbol"] = vprops[i][:symbol]
            d["coords"] = vprops[i][:coords]
            # If prop block exists, any annotations in atom blocks will be overwritten
            d["charge"] = get(ps, :CHG, 0)
            d["multiplicity"] = get(ps, :RAD, 1)
            d["isotope"] = get(ps, :ISO, 0)
            vprops[i] = V(d)
        end
    elseif ver === :v3
        readuntil(io, ctab_only ? "M  V30 END CTAB" : "M  END")
        readline(io)
    end
    mol = T(edges, vprops, eprops; NamedTuple((k, v) for (k, v) in config)...)
    # TODO: catch compatibility warnings
    if ver === :v3
        mol[:logs]["warning_sdfilereader"] = "SDFile V3000 stereochemistry is not supported"
    end
    return mol
end


function parse_rxn(::Type{T}, io::IO, config::Dict{Symbol,Any}) where T <: AbstractReaction
    line1 = readline(io)  # $RXN
    startswith(line1, "\$RXN") || error("rdfile parse error - \$RXN token not found")
    ver = line1 == "\$RXN V3000" ? :v3 : :v2
    line2 = readline(io)  # name line, not implemented
    line3 = readline(io)  # format properties, not implemented
    # program_name = line2[7:15]
    # datetime = line2[16:27]
    line4 = readline(io)  # comment line, not implemented
    count_line = readline(io)
    rxn = T()
    if ver === :v2
        rcount = parse(Int, count_line[1:3])
        pcount = parse(Int, count_line[4:6])
    elseif ver === :v3
        sp = split(count_line)
        rcount = parse(Int, sp[4])
        pcount = parse(Int, sp[5])
    end
    for _ in 1:rcount
        ver === :v2 && (readuntil(io, "\$MOL"); readline(io))
        push!(rxn.reactants, parse_ctab(eltype(T), io, config))
    end
    for _ in 1:pcount
        ver === :v2 && (readuntil(io, "\$MOL"); readline(io))
        push!(rxn.products, parse_ctab(eltype(T), io, config))
    end
    ver === :v3 && (readuntil(io, "M  END"); readline(io))
    return rxn
end


function parse_options(io::IO)
    options = OrderedDict{String,String}()
    while true
        eof(io) && break
        line = readline(io)
        line == "\$\$\$\$" && break
        # Some inappropriate signs are accepted for practical use
        m = match(r">.*?<([\w -.%=/]+)>", line)
        if m !== nothing
            lines = []
            while true
                line_ = readline(io)
                line_ == "" && break
                push!(lines, line_)
            end
            options[m[1]] = join(lines, "\n")
        end
    end
    return options
end

function parse_rdf_options(io::IO)
    options = OrderedDict{String,String}()
    while true
        eof(io) && break
        mark(io)
        line = readline(io)
        startswith(line, "\$") || continue  # multiline datum not supported
        if startswith(line, "\$MFMT") || startswith(line, "\$RFMT")
            reset(io)
            break
        end
        unmark(io)
        dtype = match(r"\$DTYPE (.*?)", line)[1]
        datum = match(r"\$DATUM (.*?)", readline(io))[1]
        options[dtype] = datum
    end
    return options
end


struct SDFileReader{T}
    io::IO
    unsupported::Symbol  # :error, :log, :ignore
    config::Dict{Symbol,Any}
end

Base.IteratorSize(::Type{SDFileReader{T}}) where T = Base.SizeUnknown()
Base.IteratorEltype(::Type{SDFileReader{T}}) where T = Base.EltypeUnknown()

function Base.iterate(reader::SDFileReader{T}, state=1) where T <: AbstractMolGraph
    eof(reader.io) && return nothing
    mol = try
        parse_ctab(T, reader.io, reader.config)
    catch e
        reader.unsupported === :error && throw(e)
        if e isa ErrorException  # Compatibility error
            reader.unsupported === :log && @info "$(e.msg) (#$(state) in sdfilereader)"
            nul = T()
            nul[:logs]["error_sdfilereader"] = e.msg
            nul
        else
            throw(e)
        end
    end
    merge!(mol.gprops.metadata, parse_options(reader.io))
    return (mol, state + 1)
end

function Base.iterate(reader::SDFileReader{T}, state=1) where T <: AbstractReaction
    eof(reader.io) && return nothing
    mark(reader.io)
    line1 = readline(reader.io)
    if line1 == "\$RDFILE 1"  # .rd
        line2 = readline(reader.io)  # $DATM
        fmt_line = readline(reader.io)  # $RFMT or $MFMT
    else
        fmt_line = line1
    end
    rxn = try
        if startswith(fmt_line, "\$MFMT")
            unmark(reader.io)
            parse_ctab(T, reader.io, reader.config)
        elseif startswith(fmt_line, "\$RFMT")
            unmark(reader.io)
            parse_rxn(T, reader.io, reader.config)
        elseif startswith(fmt_line, "\$RXN")  # .rxn
            reset(reader.io)
            parse_rxn(T, reader.io, reader.config)
        end
    catch e
        reader.unsupported === :error && throw(e)
        if e isa ErrorException  # Compatibility error
            reader.unsupported === :log && @info "$(e.msg) (#$(state) in rdfilereader)"
            T()
        else
            throw(e)
        end
    end
    rxn === nothing && error("rdfile parse error - invalid token $(fmt_line)")
    merge!(rxn.rprops.metadata, parse_rdf_options(reader.io))
    return (rxn, state + 1)
end


"""
    sdfilereader(file::IO)
    sdfilereader(path::AbstractString)

Read SDFile data from input stream (or a file path as a string) and return a
lazy iterator that yields molecule objects.

`sdfilereader` does not stop and raise errors when an erroneous or incompatible
SDFile block is read but produces an error message and yields an empty molecule.
If this behavior is not desirable, you can use the customized supplier function
instead of default supplier `nohaltsupplier`

"""
sdfilereader(
        ::Type{T}, file::IO; unsupported=:log,
        config=Dict{Symbol,Any}(:on_init => sdf_on_init!, :on_update => sdf_on_update!), kwargs...
    ) where T <: AbstractMolGraph = SDFileReader{T}(file, unsupported, config)
sdfilereader(file::IO; kwargs...) = sdfilereader(SDFMolGraph, file; kwargs...)
sdfilereader(::Type{T}, path::AbstractString; kwargs...
    ) where T <: AbstractMolGraph = sdfilereader(T, open(path); kwargs...)
sdfilereader(path::AbstractString; kwargs...) = sdfilereader(SDFMolGraph, open(path); kwargs...)

rdfilereader(
        ::Type{T}, file::IO; unsupported=:log,
        config=Dict{Symbol,Any}(:on_init => sdf_on_init!, :on_update => sdf_on_update!), kwargs...
    ) where T <: AbstractReaction = SDFileReader{T}(file, unsupported, config)
rdfilereader(file::IO; kwargs...) = rdfilereader(Reaction{SDFMolGraph}, file; kwargs...)
rdfilereader(::Type{T}, path::AbstractString; kwargs...
    ) where T <: AbstractReaction = rdfilereader(T, open(path); kwargs...)
rdfilereader(path::AbstractString; kwargs...) = rdfilereader(Reaction{SDFMolGraph}, open(path); kwargs...)


struct SDFileScanner{T<:AbstractMolGraph}
    io::IO
end

Base.IteratorSize(::Type{SDFileScanner{T}}) where T = Base.SizeUnknown()
Base.IteratorEltype(::Type{SDFileScanner{T}}) where T = Base.EltypeUnknown()

function Base.iterate(reader::SDFileScanner{T}, state=1) where T <: AbstractMolGraph
    eof(reader.io) && return nothing
    lines = []
    while true
        line = readline(reader.io)
        push!(lines, line)
        line == "M  END" && break
    end
    data = parse_options(reader.io)
    data["molblock"] = join(lines, "\n")
    return (data, state + 1)
end


"""
    sdfilescanner(file::IO)
    sdfilescanner(path::AbstractString)

Read SDFile data from input stream (or a file path as a string) and return a
lazy iterator that yields metadata.
"""
sdfilescanner(::Type{T}, file::IO; kwargs...
    ) where T <: AbstractMolGraph = SDFileScanner{T}(file)
sdfilescanner(file::IO; kwargs...) = sdfilescanner(SDFMolGraph, file; kwargs...)
sdfilescanner(::Type{T}, path::AbstractString; kwargs...
    ) where T <: AbstractMolGraph = sdfilescanner(T, open(path); kwargs...)
sdfilescanner(path::AbstractString; kwargs...) = sdfilescanner(SDFMolGraph, open(path); kwargs...)


"""
    sdftomol(::Type{T}, io::IO) -> T
    sdftomol(io::IO) -> SDFMolGraph
    sdftomol(::Type{T}, path::AbstractString) -> T
    sdftomol(path::AbstractString) -> SDFMolGraph

Read a SDFile(.sdf or .mol) and parse it into a molecule object with the given type. The given
argument should be a file input stream or a file path.
"""
sdftomol(::Type{T}, io::IO; kwargs...
    ) where T <: AbstractMolGraph = iterate(sdfilereader(T, io, unsupported=:error; kwargs...))[1]
sdftomol(file::IO; kwargs...) = sdftomol(SDFMolGraph, file; kwargs...)
sdftomol(::Type{T}, path::AbstractString; kwargs...) where T <: AbstractMolGraph = sdftomol(T, open(path); kwargs...)
sdftomol(path::AbstractString; kwargs...) = sdftomol(SDFMolGraph, open(path); kwargs...)


"""
    rxntoreaction(::Type{T}, io::IO) -> T
    rxntoreaction(io::IO) -> Reaction{SDFMolGraph}
    rxntoreaction(::Type{T}, path::AbstractString) -> T
    rxntoreaction(path::AbstractString) -> Reaction{SDFMolGraph}

Read a RXN file and parse it into a reaction object with the given type. The given
argument should be a file input stream or a file path.
"""
rxntoreaction(::Type{T}, io::IO; kwargs...
    ) where T <: AbstractReaction = iterate(rdfilereader(T, io, unsupported=:error; kwargs...))[1]
rxntoreaction(file::IO; kwargs...) = rxntoreaction(Reaction{SDFMolGraph}, file; kwargs...)
rxntoreaction(::Type{T}, path::AbstractString; kwargs...
    ) where T <: AbstractReaction = rxntoreaction(T, open(path); kwargs...)
rxntoreaction(path::AbstractString; kwargs...) = rxntoreaction(Reaction{SDFMolGraph}, open(path); kwargs...)
