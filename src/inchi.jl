#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#
const INCHI_StereoType_None       = 0
const INCHI_StereoType_DoubleBond = 1
const INCHI_StereoType_Tetrahedral= 2
const INCHI_StereoType_Allene     = 3

const INCHI_PARITY_NONE      = 0
const INCHI_PARITY_ODD       = 1   # 'o'
const INCHI_PARITY_EVEN      = 2   # 'e'
const INCHI_PARITY_UNKNOWN   = 3   # 'u'
const INCHI_PARITY_UNDEFINED = 4   # '?', rarely used

const NO_ATOM = -1             # from header: NO_ATOM == -1

const MAXVAL = 20
const ATOM_EL_LEN = 6
const NUM_H_ISOTOPES = 3
const AT_NUM = Int16
const S_CHAR = Int8
struct inchi_Atom
    x::Cdouble
    y::Cdouble
    z::Cdouble
    neighbor::NTuple{MAXVAL,AT_NUM}
    bond_type::NTuple{MAXVAL,S_CHAR}
    bond_stereo::NTuple{MAXVAL,S_CHAR}
    elname::NTuple{ATOM_EL_LEN,Cchar}
    num_bonds::AT_NUM
    num_iso_H::NTuple{NUM_H_ISOTOPES+1,S_CHAR}
    isotopic_mass::AT_NUM
    radical::S_CHAR
    charge::S_CHAR
end

mutable struct inchi_Output
    szInChI::Cstring
    szAuxInfo::Cstring
    szMessage::Cstring
    szLog::Cstring
    inchi_Output() = new(C_NULL, C_NULL, C_NULL, C_NULL)
end

mutable struct inchi_InputINCHI
    # the caller is responsible for the data allocation and deallocation
    szInChI::Cstring       # InChI ASCIIZ string to be converted to a strucure
    szOptions::Cstring     # InChI options: space-delimited; each is preceded by
                           # '/' or '-' depending on OS and compiler
end

struct inchi_Stereo0D
    neighbor::NTuple{4, AT_NUM} # always 4 atoms
    central_atom::AT_NUM        # central atom index (or NO_ATOM)
    type::S_CHAR                # inchi_StereoType0D
    parity::S_CHAR              # inchi_StereoParity0D (possibly combined)
end

mutable struct inchi_Input
    # the caller is responsible for the data allocation and deallocation
    atom::Ptr{inchi_Atom}            # array of num_atoms elements
    stereo0D::Ptr{inchi_Stereo0D}        # array of num_stereo0D 0D stereo elements or NULL
    szOptions::Cstring          # InChI options: space-delimited; each is preceded by
                                # '/' or '-' depending on OS and compiler
    num_atoms::Cshort           # number of atoms in the structure < MAX_ATOMS
    num_stereo0D::Cshort        # number of 0D stereo elements
end

mutable struct inchi_InputEx
    # the caller is responsible for the data allocation and deallocation
    atom::Ptr{inchi_Atom}          # array of num_atoms elements
    stereo0D::Ptr{inchi_Stereo0D}  # array of num_stereo0D 0D stereo elements or NULL
    szOptions::Cstring             # InChI options: space-delimited; each is preceded by
                                   # '/' or '-' depending on OS and compiler
    num_atoms::Cshort              # number of atoms in the structure < MAX_ATOMS
    num_stereo0D::Cshort           # number of 0D stereo elements
    polymer::Ptr{Cvoid}
    v3000::Ptr{Cvoid}
end

mutable struct inchi_OutputStruct
    # the caller is responsible for the data allocation and deallocation
    atom::Ptr{inchi_Atom}        # array of num_atoms elements
    stereo0D::Ptr{inchi_Stereo0D}    # array of num_stereo0D 0D stereo elements or NULL
    num_atoms::Cshort;      # number of atoms in the structure < MAX_ATOMS
    num_stereo0D::Cshort    # number of 0D stereo elements
    szMessage::Cstring      # Error/warning ASCIIZ message
    szLog::Cstring          # log-file ASCIIZ string, contains a human-readable list
                            # of recognized options and possibly an Error/warn message
    warningflags::NTuple{2, NTuple{2, Clong}} # warnings, see INCHIDIFF in inchicmp.h                        */
                            # [x][y]:
                            # x=0 => Reconnected if present in InChI
                            #         otherwise Disconnected/Normal
                            # x=1 => Disconnected layer if Reconnected layer is present
                            # y=1 => Main layer or Mobile-H
                            # y=0 => Fixed-H layer

    inchi_OutputStruct() = new(C_NULL, C_NULL, 0, 0, C_NULL, C_NULL, ((0, 0), (0, 0)))
end

mutable struct inchi_OutputStructEx
    # the caller is responsible for the data allocation and deallocation
    atom::Ptr{inchi_Atom}          # array of num_atoms elements
    stereo0D::Ptr{inchi_Stereo0D}  # array of num_stereo0D 0D stereo elements or NULL
    num_atoms::Cshort;             # number of atoms in the structure < MAX_ATOMS
    num_stereo0D::Cshort           # number of 0D stereo elements
    szMessage::Cstring             # Error/warning ASCIIZ message
    szLog::Cstring                 # log-file ASCIIZ string, contains a human-readable list
                                   # of recognized options and possibly an Error/warn message
    warningflags::NTuple{2, NTuple{2, Clong}} # warnings, see INCHIDIFF in inchicmp.h                        */
                            # [x][y]:
                            # x=0 => Reconnected if present in InChI
                            #         otherwise Disconnected/Normal
                            # x=1 => Disconnected layer if Reconnected layer is present
                            # y=1 => Main layer or Mobile-H
                            # y=0 => Fixed-H layer
    polymer::Ptr{Cvoid}
    v3000::Ptr{Cvoid}

    inchi_OutputStructEx() = new(C_NULL, C_NULL, 0, 0, C_NULL, C_NULL, ((0, 0), (0, 0)), C_NULL, C_NULL)
end

function opt_array(options::String)
    isempty(options) ? SubString{String}[] : lstrip.(split(options, ' ', keepempty = false), Ref(['-', '/']))
end

function opt_string(options::Vector{<:AbstractString})
    join((Sys.iswindows() ? '/' : '-') .* options, ' ')
end

function unsafe_info(cs::Cstring, title::String = "")
    if cs != C_NULL
        @info string(title, isempty(title) ? "" : ": ", unsafe_string(cs))
    end
end

function report_output(output, verbose::Bool)
    if output.szInChI == C_NULL || verbose
        output.szInChI == C_NULL && @info "InChI error with \$inchi or \$options"
        unsafe_info(output.szMessage, "message")
        unsafe_info(output.szMessage, "log")
    end
end


"""
    inchi(mol::SimpleMolGraph; options::String = "", verbose::Bool = false) -> Union{String,Nothing}
    inchi(molblock::String; options::String = "", verbose::Bool = false) -> Union{String,Nothing}

Generate InChI string from molblock string or molecule.

Options, e.g. "SNon" for 'no stereo information' are specified in https://github.com/mojaie/libinchi/blob/master/INCHIBASE/src/inchiapi.h
"""
inchi(mol::SimpleMolGraph; kwargs...) = inchi(mol, vproptype(mol), eproptype(mol); kwargs...)
inchi(
    mol::SimpleMolGraph, ::Type{<:StandardAtom}, ::Type{<:StandardBond}; kwargs...
) = inchi(printv2mol(mol); kwargs...)

function inchi(molblock::String; options::String = "", verbose::Bool = false)
    # support the correct options format depending on OS
    # add a timeout of 60s per molecule
    opts = opt_array(options)
    any(occursin.(r"^Wm?\d+$", opts)) || push!(opts, "W60")
    options = opt_string(opts)

    output = inchi_Output()
    @ccall libinchi.MakeINCHIFromMolfileText(
        molblock::Cstring, options::Cstring, output::Ref{inchi_Output})::Int32
    report_output(output, verbose)

    res = output.szInChI == C_NULL ? nothing : unsafe_string(output.szInChI)

    # Free string buffers allocated by MakeINCHIFromMolfileText
    @ccall libinchi.FreeINCHI(output::Ref{inchi_Output})::Cvoid
    return res
end



"""
    inchikey(inchi::String) -> Union{String,Nothing}
    inchikey(mol::MolGraph) -> Union{String,Nothing}

Generate InChI key from InChI string or molecule
"""
function inchikey(inchi::Union{String,Nothing})
    inchi === nothing && return nothing
    ikeybuf = pointer(Vector{UInt8}(undef, 256))
    # TODO: need extra buffer?
    xtra1buf = pointer(Vector{UInt8}(undef, 256))
    xtra2buf = pointer(Vector{UInt8}(undef, 256))
    @ccall libinchi.GetINCHIKeyFromINCHI(
        inchi::Cstring, 1::Int32, 1::Int32,
        ikeybuf::Cstring, xtra1buf::Cstring, xtra2buf::Cstring)::Int32
    return unsafe_string(ikeybuf)
end

inchikey(mol::SimpleMolGraph) = inchikey(inchi(mol))

"""
    inchitosdf(inchi::String; options::String = "") -> Union{String,Nothing}

Generate sdf string from inchi string, `options` are specified in https://github.com/mojaie/libinchi/blob/master/INCHI_BASE/src/inchi_api.h
"""
function inchitosdf(inchi::String; options::String = "", verbose::Bool = false)
    # support the correct options format depending on OS
    opts = opt_array(options)
    # add a timeout of 60s per molecule
    any(occursin.(r"^Wm?\d+$", opts)) || push!(opts, "W60")
    # switch output to sdf format
    "OutputSDF" ∈ opts || push!(opts, "OutputSDF")
    options = opt_string(opts)

    structure = inchi_OutputStructEx()
    inchi_input = inchi_InputINCHI(Base.unsafe_convert(Cstring, inchi), Base.unsafe_convert(Cstring, options))
    @ccall libinchi.GetStructFromINCHIEx(
        inchi_input::Ref{inchi_InputINCHI}, structure::Ref{inchi_OutputStructEx})::Int32

    input = inchi_InputEx(
        structure.atom,
        structure.stereo0D,
        Base.unsafe_convert(Cstring, options),
        structure.num_atoms,
        structure.num_stereo0D,
        structure.polymer,
        structure.v3000
    )
    output = inchi_Output()
    @ccall libinchi.GetINCHIEx(
        input::Ref{inchi_InputEx}, output::Ref{inchi_Output})::Cint
    report_output(output, verbose)

    res = output.szInChI == C_NULL ? nothing : unsafe_string(output.szInChI)

    # Free buffers allocated by GetStructFromINCHIEx and GetINCHI
    @ccall libinchi.FreeStructFromINCHIEx(structure::Ref{inchi_OutputStructEx})::Cvoid
    @ccall libinchi.FreeINCHI(output::Ref{inchi_Output})::Cvoid

    return res
end

function inchi_atom_symbol(a::inchi_Atom)
    # extract main atom symbol
    s = String(reinterpret(UInt8, collect(Iterators.takewhile(!=(0x00), a.elname))))
    return s
end

function bonds_from_atoms(atoms::Vector{inchi_Atom})
    bonds = Set{Tuple{Int,Int,Int}}()
    for (i, a) in enumerate(atoms)
        for k in 1:a.num_bonds
            j = a.neighbor[k] + 1        # C -> Julia 1-based
            btype = a.bond_type[k]
            if j > 0  # skip invalid neighbor
                push!(bonds, (min(i,j), max(i,j), btype))
            end
        end
    end
    return collect(bonds)
end

function expand_inchi_atoms(atoms::Vector{inchi_Atom})
    expanded_atoms = Vector{String}()
    bonds = Vector{Tuple{Int,Int,Int}}()

    for (i, a) in enumerate(atoms)
        # Add main atom
        push!(expanded_atoms, inchi_atom_symbol(a))

        # Keep track of next atom index
        parent_idx = length(expanded_atoms)

        # Add implicit H/D/T as explicit atoms
        for iso_index in 1:NUM_H_ISOTOPES
            n = a.num_iso_H[iso_index+1]
            iso_symbol = ["H", "D", "T"][iso_index]
            for _ in 1:n
                child_idx = length(expanded_atoms) + 1
                push!(expanded_atoms, iso_symbol)
                push!(bonds, (parent_idx, child_idx, 1))  # single bond to parent
            end
        end
    end

    # Add original bonds between heavy atoms
    heavy_bonds = bonds_from_atoms(atoms)
    append!(bonds, heavy_bonds)

    return expanded_atoms, bonds
end

function expand_inchi_atoms(atoms::Vector{inchi_Atom})
    expanded_atoms = Vector{String}()
    bonds = Vector{Tuple{Int,Int,Int}}()

    for (i, a) in enumerate(atoms)
        # Add main atom
        push!(expanded_atoms, inchi_atom_symbol(a))

        # Keep track of next atom index
        parent_idx = length(expanded_atoms)

        # Add implicit H/D/T as explicit atoms
        for iso_index in 1:NUM_H_ISOTOPES
            n = a.num_iso_H[iso_index+1]
            iso_symbol = ["H", "D", "T"][iso_index]
            for _ in 1:n
                child_idx = length(expanded_atoms) + 1
                push!(expanded_atoms, iso_symbol)
                push!(bonds, (parent_idx, child_idx, 1))  # single bond to parent
            end
        end
    end

    # Add original bonds between heavy atoms
    heavy_bonds = bonds_from_atoms(atoms)
    append!(bonds, heavy_bonds)

    return expanded_atoms, bonds
end

function process_inchi_stereo!(g::T, structure::inchi_OutputStructEx, orig_to_expanded::Vector{Int64}) where T <: MolGraph
    ET = edgetype(T)
    stereocenters = Dict{eltype(T), Tuple{Int64,Int64,Int64,Bool}}()
    stereobonds   = Dict{ET, Tuple{Int64,Int64,Bool}}()

    if structure.stereo0D != C_NULL && Int(structure.num_stereo0D) > 0
        stereos = unsafe_wrap(Vector{inchi_Stereo0D}, structure.stereo0D, Int(structure.num_stereo0D))

        # helper: map original 0-based atom idx -> expanded 1-based index (Int64)
        map_orig = orig_idx -> begin
            if orig_idx == NO_ATOM
                return nothing
            end
            oi = Int(orig_idx) + 1
            if oi < 1 || oi > length(orig_to_expanded)
                return nothing
            end
            return Int64(orig_to_expanded[oi])
        end

        for s in stereos
            stype = Int(s.type)
            # parity may be packed: lower 3 bits = connected parity, bits 3..5 = disconnected parity
            p_raw = Int(UInt8(s.parity))
            p_conn = p_raw & 0x7
            p_disc = (p_raw >> 3) & 0x7
            parity_used = (p_conn != 0 ? p_conn : p_disc)

            if stype == INCHI_StereoType_Tetrahedral || stype == INCHI_StereoType_Allene
                if s.central_atom == NO_ATOM
                    continue
                end
                center = map_orig(s.central_atom)
                if center === nothing
                    continue
                end
                # neighbors W,X,Y,Z -> we use W as looking_from, X/Y as v1/v2
                nbrs = (map_orig(s.neighbor[1]), map_orig(s.neighbor[2]),
                        map_orig(s.neighbor[3]), map_orig(s.neighbor[4]))
                # need at least W,X,Y present
                if any(x -> x === nothing, (nbrs[1], nbrs[2], nbrs[3]))
                    continue
                end
                looking_from = nbrs[1]::Int64
                v1 = nbrs[2]::Int64
                v2 = nbrs[3]::Int64
                # InChI docs: EVEN => clockwise when seen from W
                is_clockwise = (parity_used == INCHI_PARITY_EVEN)

                stereocenters[center] = (looking_from, v1, v2, is_clockwise)

            elseif stype == INCHI_StereoType_DoubleBond
                # neighbor order A1,A2,A3,A4 ; stereo about bond between A2-A3
                nbrs = (map_orig(s.neighbor[1]), map_orig(s.neighbor[2]),
                        map_orig(s.neighbor[3]), map_orig(s.neighbor[4]))
                if nbrs[2] === nothing || nbrs[3] === nothing
                    continue
                end
                nbrs = s.neighbor
                a1, a2, a3, a4 = Int.(nbrs)

                # canonical edge key:
                bkey = a2 < a3 ? ET(a2, a3) : ET(a3, a2)
                is_cis = (parity_used == INCHI_PARITY_ODD)  # odd => cis, even => trans
                stereobonds[bkey] = (a1 === nothing ? Int64(-1) : a1::Int64,
                                    a4 === nothing ? Int64(-1) : a4::Int64,
                                    is_cis)
            end
        end
    end

    merge!(g.gprops.stereocenter, stereocenters)
    merge!(g.gprops.stereobond, stereobonds)
    return g
end

"""
    function inchitomol(inchi::String; options = "", verbose = false)

Generate molecule from inchi string, `options` are specified in https://github.com/mojaie/libinchi/blob/master/INCHI_BASE/src/inchi_api.h
"""
function inchitomol(::Type{T}, inchi::AbstractString;
    options::String = "",
    config::Union{Nothing,Dict{Symbol,Any}} = nothing,
    stereo::Bool = true
) where T <: MolGraph
    inchi isa String || (inchi = "$inchi")

    structure = inchi_OutputStructEx()
    inchi_input = inchi_InputINCHI(Base.unsafe_convert(Cstring, inchi),
                                   Base.unsafe_convert(Cstring, options))

    ret = @ccall libinchi.GetStructFromINCHIEx(
        inchi_input::Ref{inchi_InputINCHI},
        structure::Ref{inchi_OutputStructEx}
    )::Int32

    if ret < 0
        error("libinchi failed with code $ret")
    end
    
    config = if config === nothing && T === SMILESMolGraph
        Dict{Symbol, Any}(:on_init => smiles_on_init!, :on_update => smiles_on_update!)
    elseif config === nothing && T === SDFMolGraph
        Dict{Symbol, Any}(:on_init => sdf_on_init!, :on_update => sdf_on_update!)
    else
        config
    end
    
    atoms = unsafe_wrap(Vector{inchi_Atom}, structure.atom, structure.num_atoms)

    # determine the MolecularGraph atom/bond types and constructors
    V = vproptype(T) === AbstractAtom ? SDFAtom : vproptype(T)
    E = eproptype(T) === AbstractBond ? SDFBond : eproptype(T)
    ET = edgetype(T)

    expanded_syms = String[]
    orig_to_expanded = Vector{Int}(undef, length(atoms))

    # heavy atoms
    for (i, a) in enumerate(atoms)
        push!(expanded_syms,
              String(reinterpret(UInt8, collect(Iterators.takewhile(!=(0x00), a.elname)))))
        orig_to_expanded[i] = length(expanded_syms)
    end

    # isotopic H/D/T
    for (i, a) in enumerate(atoms)
        parent_idx = orig_to_expanded[i]
        for iso_index in 1:NUM_H_ISOTOPES
            n = Int(a.num_iso_H[iso_index+1])
            iso_sym = ["H", "D", "T"][iso_index]
            for _ in 1:n
                push!(expanded_syms, iso_sym)
            end
        end
    end

    # build bonds
    bonds = Tuple{Int,Int,Int}[]
    first_isotope_index = length(atoms) + 1
    next_iso_idx = first_isotope_index
    for (i, a) in enumerate(atoms)
        parent_idx = orig_to_expanded[i]
        for iso_index in 1:NUM_H_ISOTOPES
            n = Int(a.num_iso_H[iso_index+1])
            for _ in 1:n
                push!(bonds, (parent_idx, next_iso_idx, 1))
                next_iso_idx += 1
            end
        end
    end

    for (i, a) in enumerate(atoms)
        for k in 1:Int(a.num_bonds)
            nb0 = Int(a.neighbor[k])
            if nb0 < 0
                continue
            end
            j = nb0 + 1
            if j > i
                push!(bonds, (orig_to_expanded[i], orig_to_expanded[j], Int(a.bond_type[k])))
            end
        end
    end

    N = length(expanded_syms)
    vprops = V[]
    edges = ET[]
    eprops = E[]

    for idx in 1:N
        if idx <= length(atoms)
            a = atoms[idx]
            symbol = expanded_syms[idx]
            d = Dict{String,Any}(
                "symbol" => symbol,
                "charge" => Int(a.charge),
                "multiplicity" => (Int(a.radical) == 0 ? 1 : Int(a.radical)),
                "isotope" => Int(a.isotopic_mass),
            )
            T === SDFMolGraph && push!(d, "coords" => [Float64(a.x), Float64(a.y), Float64(a.z)])
            push!(vprops, V(d))
        else
            symbol = expanded_syms[idx]
            d = Dict(
                "symbol" => symbol,
                "charge" => 0,
                "multiplicity" => 1,
                "isotope" => 0
            )
            T === SDFMolGraph && push!(d, "coords" => [0.0, 0.0, 0.0])
            push!(vprops, V(d))
        end
    end

    for (i, j, order) in bonds
        ekey = i < j ? ET(i, j) : ET(j, i)
        ed = Dict{String, Any}("order" => order)
        push!(edges, ekey)
        push!(eprops, E(ed))
    end

    g = T(edges, vprops, eprops; NamedTuple((k, v) for (k, v) in config)...)

    stereo && process_inchi_stereo!(g, structure, orig_to_expanded)
    
    coordgen!(g)
    return g
end

function inchitomol(inchi::String; options::String = "", stereo::Bool = true)
    inchitomol(SDFMolGraph, inchi; options, stereo)
end