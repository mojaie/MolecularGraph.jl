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
    inchitosdf(inchi::AbstractString; options::String = "", config::Union{Nothing,Dict{Symbol,Any}})

Generate sdf string from inchi string.
This new version goes via parsing to a MolGraph and exporting as sdf.
Parsing `options` are specified in https://github.com/mojaie/libinchi/blob/master/INCHI_BASE/src/inchi_api.h
Coordinates are generated via coordgen (Schrödinger coordgenlibs).
"""
function inchitosdf(inchi::AbstractString; options::String = "", config::Union{Nothing,Dict{Symbol,Any}} = nothing)
    printv2mol(inchitomol(inchi; options, config))
end

# decode parity byte (connected in low 3 bits, disconnected in bits 3..5)
# returns the parity value to use (0..4)
decode_parity(p::Integer) = begin
    p_raw  = UInt8(p)                    # make sure unsigned before bit ops
    p_conn = Int(p_raw & 0x07)           # low 3 bits
    p_disc = Int((p_raw >> 3) & 0x07)    # next 3 bits
    p_conn != 0 ? p_conn : p_disc        # prefer connected parity if present
end

function process_inchi_stereo!(g::T, structure::inchi_OutputStructEx) where T <: MolGraph
    ET = edgetype(T)
    
    stereocenters = Dict{Int64, Tuple{Int64,Int64,Int64,Bool}}()   # center => (looking_from,v1,v2,is_clockwise)
    stereobonds   = Dict{Edge{Int64}, Tuple{Int64,Int64,Bool}}()   # edge => (a1,a2,is_cis)

    atoms = unsafe_wrap(Vector{inchi_Atom}, structure.atom, structure.num_atoms)
    stereos = unsafe_wrap(Vector{inchi_Stereo0D}, structure.stereo0D, structure.num_stereo0D)
    n_heavy = length(atoms)  # number of heavy atoms (without isotopic H/D/T)
    is_heavy(idx) = 0 < idx <= n_heavy # for 1-based indices

    for s in stereos
        parity = decode_parity(s.parity)
        parity == INCHI_PARITY_NONE || parity == INCHI_PARITY_UNKNOWN && continue  # nothing to set
        parity == INCHI_PARITY_UNDEFINED && continue  # explicit undefined -> skip

        # convert to 1-based indexing
        nbrs = s.neighbor .+ 1 
        center = s.central_atom + 1

        if s.type == INCHI_StereoType_Tetrahedral
            center == 0 && continue  # sanity check

            # InChI neighbors order corresponds to W,X,Y,Z (see InChI comments).
            w, x, y, z = nbrs   # these are already expanded indices (or NO_ATOM)
            is_clockwise = (parity == INCHI_PARITY_EVEN)   # 'e' == clockwise per docs
            stereocenters[center] = (w, x, y, is_clockwise)
        elseif s.type == INCHI_StereoType_DoubleBond
            all(is_heavy.(nbrs)) || continue  # skip disconnected or missing substituents

            # the double bond is between the middle two atoms
            # by InChI convention, neighbor[1,2] on one end, neighbor[3,4] on the other
            # determine the actual bond by looking up the shared bond in the graph
            bond_candidates = [u_edge(g, nbrs[i], nbrs[j]) for i in 1:2, j in 3:4]
            # pick the one that actually exists in your bond list
            edges = collect(g.eprops)  # materialize iterator
            bond_key_index = findfirst(e -> e[1] in bond_candidates, edges)
            bond_key_index === nothing && continue
            bond_key = edges[bond_key_index][1]

            is_trans = if parity == INCHI_PARITY_ODD
                true      # InChI uses ODD -> trans
            elseif parity == INCHI_PARITY_EVEN
                false     # EVEN -> cis
            else
                # we currently don't get here, because we filtered out NONE/UNKNOWN/UNDEFINED above
                # once we can handle undetermined stereo, we need to implement that logic here
                continue   # unknown
            end
            
            # store stereobond
            stereobonds[bond_key] = (nbrs[1], nbrs[3], is_trans)
        elseif s.type == INCHI_StereoType_Allene
            @info "Allene stereo, not yet tested"
            # Tried testing with `mol = smilestomol("C[C@]=C=C(C)F")`,
            # but the InChI generated from it does not contain stereo information.
            # Allene (axial) handling
            # InChI gives neighbors as W,X,Y,Z where W/X are substituents on one end
            # and Y/Z on the other. central_atom is the central cumulated C.
            # We need to pick a viewing atom and two reference atoms to encode axial chirality.
            center == 0 && continue

            # For an allene to be chiral, both ends must have two substituents (i.e. W/X and Y/Z present)
            all(is_heavy.(nbrs)) || continue
            w, x, y, z = nbrs

            is_clockwise = (parity == INCHI_PARITY_EVEN)

            # store into same stereocenters dict so coordgen!/rendering code can use it
            stereocenters[center] = (w, x, y, is_clockwise)
        else
            @warn "Unknown stereo type $(s.type), skipping"
        end
    end

    # merge into mol.gprops the same way you already did:
    merge!(g.gprops.stereocenter, stereocenters)
    merge!(g.gprops.stereobond, stereobonds)
    return g
end

"""
    function inchitomol(inchi::String; options = "", config::Union{Nothing,Dict{Symbol,Any}} = nothing)

Generate molecule from inchi string, `options` are specified in https://github.com/mojaie/libinchi/blob/master/INCHI_BASE/src/inchi_api.h
`config` is for internal or advanced use only. Maybe removed in a future release.
"""
function inchitomol(::Type{T}, inchi::AbstractString;
    options::String = "",
    config::Union{Nothing,Dict{Symbol,Any}} = nothing
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
    
    # determine the MolecularGraph atom/bond types and constructors
    V = vproptype(T) === AbstractAtom ? SDFAtom : vproptype(T)
    E = eproptype(T) === AbstractBond ? SDFBond : eproptype(T)
    ET = edgetype(T)

    atoms = unsafe_wrap(Vector{inchi_Atom}, structure.atom, structure.num_atoms)

    # heavy atoms
    expanded_syms = String[
        String(reinterpret(UInt8, collect(Iterators.takewhile(!=(0x00), a.elname))))
    for a in atoms]

    # isotopic H/D/T
    for a in atoms
        # a.num_iso_H is a NTuple{4,S_CHAR}, where index 1 = non-specified, 2 = H, 3 = D, 4 = T
        # let's add symbols only for explicit isotopes [1H], [D] and [T]
        # as NUM_H_ISOTOPES = 3, we only need to check indices 2, 3, and 4
        for iso_index in 1:NUM_H_ISOTOPES
            n = a.num_iso_H[iso_index + 1]
            iso_sym = ["H", "D", "T"][iso_index]
            for _ in 1:n
                push!(expanded_syms, iso_sym)
            end
        end
    end

    # build bonds to isotopic H/D/T
    bonds = Tuple{Int,Int,Int}[]
    first_isotope_index = length(atoms) + 1
    next_iso_idx = first_isotope_index
    for (parent_idx, a) in enumerate(atoms)
        for iso_index in 1:NUM_H_ISOTOPES
            n = a.num_iso_H[iso_index + 1]
            for _ in 1:n
                push!(bonds, (parent_idx, next_iso_idx, 1))
                next_iso_idx += 1
            end
        end
    end

    # heavy atom bonds
    for (i, a) in enumerate(atoms)
        for k in 1:a.num_bonds
            j = a.neighbor[k] + 1
            j > i && push!(bonds, (i, j, Int(a.bond_type[k])))
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
            :coords in fieldnames(V) && push!(d, "coords" => [Float64(a.x), Float64(a.y), Float64(a.z)])
            push!(vprops, V(d))
        else
            symbol = expanded_syms[idx]
            d = Dict(
                "symbol" => symbol,
                "charge" => 0,
                "multiplicity" => 1,
                "isotope" => 0
            )
            :coords in fieldnames(V) && push!(d, "coords" => [0.0, 0.0, 0.0])
            push!(vprops, V(d))
        end
    end

    for (i, j, order) in bonds
        ekey = i < j ? ET(i, j) : ET(j, i)
        push!(edges, ekey)
        push!(eprops, E(; order))
    end

    g = T(edges, vprops, eprops; NamedTuple((k, v) for (k, v) in config)...)

    process_inchi_stereo!(g, structure)
    
    @ccall libinchi.FreeStructFromINCHIEx(structure::Ref{inchi_OutputStructEx})::Cvoid

    coordgen!(g)
    remove_hydrogens!(g)
    return g
end

function inchitomol(inchi::String; options::String = "", config::Union{Nothing,Dict{Symbol,Any}} = nothing)
    inchitomol(SDFMolGraph, inchi; options, config)
end