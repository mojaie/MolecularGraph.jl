# CDXML Parser for MolecularGraph.jl
# Based on ChemDraw CDXML format specification
#
# This parser reads CDXML (ChemDraw XML) files and converts them to MolecularGraph.jl
# MolGraph objects. It handles molecules with atoms, bonds, coordinates, charges,
# stereochemistry, and other chemical properties.

import EzXML

# =============================================================================
# Atom and Bond types for CDXML
# =============================================================================

"""CDXMLAtom represents an atom from CDXML with coordinates and properties"""
struct CDXMLAtom <: StandardAtom
    symbol::Symbol
    charge::Int
    multiplicity::Int
    isotope::Int
    isaromatic::Bool
    coords::Union{Vector{Float64}, Nothing}  # concrete element type (for `juliac --trim`)
    num_hs::Int

    function CDXMLAtom(
        symbol::Union{AbstractString, Symbol},
        charge::Int,
        multiplicity::Int,
        isotope::Int,
        isaromatic::Bool,
        coords::Union{Vector, Nothing},
        num_hs::Int
    )
        haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("sdfile parse error - unsupported atom symbol '$symbol'")
        new(Symbol(symbol), charge, multiplicity, isotope, isaromatic,
            isnothing(coords) ? nothing : Vector{Float64}(coords), num_hs)
    end
end

function CDXMLAtom(
    symbol::Union{AbstractString, Symbol};
    charge::Int=0, 
    multiplicity::Int=1, 
    isotope::Int=0,
    isaromatic::Bool=false,
    coords::Union{Vector, Nothing}=nothing,
    num_hs::Int=-1
)
    haskey(ATOMSYMBOLMAP, Symbol(symbol)) || error("sdfile parse error - unsupported atom symbol $(symbol)")
    CDXMLAtom(symbol, charge, multiplicity, isotope, isaromatic, coords, num_hs)
end

"""CDXMLBond represents a bond from CDXML with order and stereochemistry"""
struct CDXMLBond <: StandardBond
    order::Int
    isaromatic::Bool
    isordered::Bool
    stereo::Symbol
    notation::Int

    function CDXMLBond(
            order::Int;
            isaromatic::Bool=false,
            isordered::Bool=false,
            stereo::Union{AbstractString,Symbol} = :unspecified,
            notation::Int = 0,
        kwargs...
    )
        isempty(kwargs) || @warn "keyword(s) `$(join(keys(kwargs), "``, `"))` not supported!"
        order > 3 && error("sdfile parse error - unsupported bond order $(order)")
        new(order, isaromatic, isordered, Symbol(stereo), notation)
    end
end

# keyword-only constructor for compatibility with other readers
function CDXMLAtom(;
    symbol::Union{AbstractString, Symbol},
    charge::Int=0, 
    multiplicity::Int=1, 
    isotope::Int=0,
    isaromatic::Bool=false,
    coords::Union{Vector, Nothing}=nothing,
    num_hs::Int=-1,
    kwargs...
)
    isempty(kwargs) || @warn "keyword(s) `$(join(keys(kwargs), "``, `"))` not supported!"
    CDXMLAtom(symbol;
        charge,
        multiplicity,
        isotope,
        isaromatic,
        coords,
        num_hs
    )
end

"""
    CDXMLMolGraph

Alias of `MolGraph{Int,CDXMLAtom,CDXMLBond}`, the default molecule type of
[`cdxmltomol`](@ref) and [`cdxmltomols`](@ref).
"""
const CDXMLMolGraph = MolGraph{Int,CDXMLAtom,CDXMLBond}

# Implement required interface methods
atom_symbol(atom::CDXMLAtom) = atom.symbol
atom_charge(atom::CDXMLAtom) = atom.charge
multiplicity(atom::CDXMLAtom) = atom.multiplicity
atom_number(atom::CDXMLAtom) = atom_number(atom.symbol)
isotope(atom::CDXMLAtom) = atom.isotope
isaromatic(atom::CDXMLAtom) = atom.isaromatic
bond_order(bond::CDXMLBond) = bond.order
isaromatic(bond::CDXMLBond) = bond.isaromatic

has_isaromatic(::Type{CDXMLAtom}) = true

aromatic_atom(a::CDXMLAtom) = CDXMLAtom(
    a.symbol, a.charge, a.multiplicity, a.isotope, true, a.coords, a.num_hs)

# Post-processing function
function mark_aromatic_atoms!(mol::MolGraph{T,V,E}) where {T, V, E}
    aromatic_atoms = Set{Int}()
    
    for (edge_key, bond) in mol.eprops
        if bond.isaromatic
            push!(aromatic_atoms, edge_key.key.src)
            push!(aromatic_atoms, edge_key.key.dst)
        end
    end
    
    for atom_idx in aromatic_atoms
        atom = mol[atom_idx]
        new_atom = aromatic_atom(atom)
        mol[atom_idx] = new_atom
    end
    
    return mol
end


# =============================================================================
# Helper functions
# =============================================================================

"""Parse element from CDXML Element attribute"""
function parse_element(elem_str::String)::Symbol
    # Try parsing as number first
    elem_num = tryparse(Int, elem_str)
    elem_num === nothing || return atom_symbol(elem_num)

    # Otherwise treat as symbol string
    return Symbol(elem_str)
end

"""Parse coordinates from CDXML p attribute (format: \"x y\" or \"x y z\")"""
function parse_coords(coords_str::String)::Vector{Float64}
    parts = split(coords_str)
    x = parse(Float64, parts[1])
    y = parse(Float64, parts[2])
    z = length(parts) > 2 ? parse(Float64, parts[3]) : 0.0
    return [x, y, z]
end

"""Parse bond order from CDXML Order attribute"""
function parse_bond_order(order_str::String)::Tuple{Int, Bool}
    # CDXML: 1=single, 2=double, 3=triple, 1.5=aromatic, 0.5=half, etc.
    if order_str == "1.5"
        return (1, true)  # aromatic
    elseif order_str == "0.5"
        return 1, false  # half bond - not sure yet how to treat
    else
        return round(Int, parse(Float64, order_str)), false
    end
end

"""Parse bond stereochemistry from CDXML Display attribute"""
function parse_bond_stereo(display_str::String)::Int
    # Common stereochemistry display types
    if display_str ∈ ("Wedge", "WedgeBegin", "Bold")
        1
    elseif display_str ∈ ("Hash", "HashBegin", "WedgedHashBegin", "Dash")
        6
    elseif display_str == "WedgeHashEnd"
        @info "Found bond type 'WedgeHashEnd', please check carefully whether stereo is correct."
        # Tip at the end, the correct location of start and end atom is taken care of by `isordered`
        6 
    elseif display_str == "Wavy"
        4
    elseif display_str == "CisTransUnknown"
        3
    else
        0
    end
end

"""Get XML attribute value or `nothing` if missing"""
attribute(node::EzXML.Node, key::AbstractString) = haskey(node, key) ? node[key] : nothing
function attribute(node::EzXML.Node, key::AbstractString, default::T)::Union{String, T} where T
    haskey(node, key) ? node[key] : default
end

# =============================================================================
# Core parsing functions
# =============================================================================

"""Parse a single atom (node) from CDXML"""
function parse_cdxml_atom(node::EzXML.Node)::Tuple{String,CDXMLAtom}
    # Get atom ID
    id = attribute(node, "id")
    if id === nothing
        error("CDXML atom node missing 'id' attribute")
    end
    
    # Parse element (defaults to Carbon if not specified)
    elem = attribute(node, "Element")
    symbol = elem === nothing ? :C : parse_element(elem)
    
    # Parse coordinates
    coords_str = attribute(node, "p")
    coords = coords_str === nothing ? nothing : parse_coords(coords_str)
    
    # Parse charge
    charge_str = attribute(node, "Charge")
    charge = charge_str === nothing ? 0 : parse(Int, charge_str)
    
    # Parse isotope
    isotope_str = attribute(node, "Isotope")
    isotope = isotope_str === nothing ? 0 : parse(Int, isotope_str)
    
    # Parse radical (multiplicity)
    radical_str = attribute(node, "Radical")
    # CDXML Radical: 0=singlet, 1=doublet, 2=triplet
    multiplicity = if radical_str === nothing
        1
    else
        radical = parse(Int, radical_str)
        radical + 1  # Convert to multiplicity
    end

    num_hs_str = attribute(node, "NumHydrogens")
    num_hs = num_hs_str === nothing ? -1 : parse(Int, num_hs_str)
    
    isaromatic = false # will be set later in mark_aromatic_atoms!()
    atom = CDXMLAtom(symbol, charge, multiplicity, isotope, isaromatic, coords, num_hs)
    return (id, atom)
end

"""Parse a single bond from CDXML"""
function parse_cdxml_bond(node::EzXML.Node)::Tuple{String,String,String,CDXMLBond}
    # Get bond ID
    id = attribute(node, "id")
    
    # Get connected atoms
    begin_atom = attribute(node, "B")
    end_atom = attribute(node, "E")
    
    if begin_atom === nothing || end_atom === nothing
        error("CDXML bond missing 'B' or 'E' attribute")
    end
    
    # Parse bond order
    order_str = attribute(node, "Order")
    order, isaromatic = order_str === nothing ? (1, false) : parse_bond_order(order_str)
    
    # Parse stereochemistry
    display_str = attribute(node, "Display")
    notation = display_str === nothing ? 0 : parse_bond_stereo(display_str)
    isordered = begin_atom < end_atom
    bond = CDXMLBond(order; notation, isaromatic, isordered)
    
    return (id === nothing ? "" : id, begin_atom, end_atom, bond)
end

# Concrete molecule type returned by the CDXML reader
# (`MolGraph{Int,CDXMLAtom,CDXMLBond}` is not concrete. `juliac --trim` needs the concrete type)
cdxml_moltype() = MolGraph{Int,CDXMLAtom,CDXMLBond,
    MolState{Int,typeof(cdxml_on_init!),typeof(cdxml_on_update!)}}

# Message of an exception without `show` (not compilable with `juliac --trim`)
function exception_message(e)
    e isa ErrorException || return string(nameof(typeof(e)))
    m = e.msg
    return m isa String ? m : "error"
end

"""Parse a fragment (molecule) from CDXML"""
function parse_cdxml_fragment(frag_node::EzXML.Node)
    atom_map = Dict{String,Int}()  # CDXML id -> graph index
    atoms = CDXMLAtom[]
    bonds = Tuple{Int,Int,CDXMLBond}[]
    
    # Parse all child nodes
    for child in EzXML.eachelement(frag_node)
        tag = EzXML.nodename(child)
        
        if tag == "n"  # Node (Atom)
            id, atom = parse_cdxml_atom(child)
            push!(atoms, atom)
            atom_map[id] = length(atoms)
            
        elseif tag == "b"  # Bond
            bond_id, begin_id, end_id, bond = parse_cdxml_bond(child)
            
            # Look up atom indices
            if haskey(atom_map, begin_id) && haskey(atom_map, end_id)
                begin_idx = atom_map[begin_id]
                end_idx = atom_map[end_id]
                push!(bonds, (begin_idx, end_idx, bond))
            else
                @warn "Bond references unknown atoms: $begin_id -> $end_id"
            end
        end
    end

    # Build MolGraph using correct constructor pattern
    if isempty(atoms)
        # Return empty molecule with correct types
        return MolGraph{Int,CDXMLAtom,CDXMLBond}(
            Edge{Int}[],
            CDXMLAtom[],
            CDXMLBond[],
            on_init = cdxml_on_init!,
            on_update = cdxml_on_update!
        )::cdxml_moltype()
    end
    
    # Create edge list for SimpleGraph
    edge_list = Edge{Int}[]
    eprop_list = CDXMLBond[]
    
    for (i, j, bond) in bonds
        push!(edge_list, i < j ? Edge(i, j) : Edge(j, i))
        push!(eprop_list, bond)
    end

    # Construct MolGraph from edge_list, vprop_list (atoms), and eprop_list
    MolGraph{Int,CDXMLAtom,CDXMLBond}(edge_list, atoms, eprop_list;
        on_init = cdxml_on_init!,
        on_update = cdxml_on_update!
    )::cdxml_moltype()
end

# Find all fragments in the document
# Fragments can be at various levels: page > fragment, or directly in document
function find_cdxml_fragments!(fragments::Vector{EzXML.Node}, node::EzXML.Node)
    for child in EzXML.eachelement(node)
        tag = EzXML.nodename(child)
        if tag == "fragment"
            push!(fragments, child)
        elseif tag == "page" || tag == "group"
            # Recursively search in pages and groups
            find_cdxml_fragments!(fragments, child)
        end
    end
    return fragments
end

# =============================================================================
# Public API
# =============================================================================

"""
    cdxmltomols(io::IO) -> Vector{MolGraph}

Parse CDXML from an input stream and return vector of molecules.
Each fragment in the CDXML document becomes a separate molecule.
"""
function cdxmltomols(cdxml_str::IO)
    mols = cdxml_moltype()[]
    
    try
        doc = EzXML.parsexml(read(cdxml_str, String))
        root = EzXML.root(doc)
        fragments = find_cdxml_fragments!(EzXML.Node[], root)
        
        # Parse each fragment into a molecule
        for frag in fragments
            try
                mol = parse_cdxml_fragment(frag)
                if length(mol.vprops) > 0  # Check if molecule has atoms
                    push!(mols, mol)
                end
            catch e
                @warn "Failed to parse fragment: $(exception_message(e))"
            end
        end
        
    catch e
        error("Failed to parse CDXML: $(exception_message(e))")
    end
    
    return mols
end

"""
    cdxmltomols(filename::String) -> Vector{MolGraph}

Read and parse a CDXML file, returning vector of molecules.
"""
function cdxmltomols(filename::String)
    open(filename) do io
        cdxmltomols(io)
    end
end

"""
    cdxmltomol(input) -> MolGraph

Parse CDXML from a file path or an input stream and return the first molecule.
Raises error if no molecules found.
"""
function cdxmltomol(input::Union{AbstractString, IO})
    mols = cdxmltomols(input)
    
    if isempty(mols)
        error("No molecules found in CDXML")
    end
    
    length(mols) > 1 && @warn """
    More than one molecule found, returning only the first one.
    For obtaining all molecules use `cdxmltomols()`.
    """
    return mols[1]
end

function cdxml_on_init!(mol::SimpleMolGraph)
    check_valence!(mol)
    coords_from_sdf!(mol)
    stereocenter_from_sdf2d!(mol)
    stereobond_from_sdf2d!(mol)
end

function cdxml_on_update!(mol::SimpleMolGraph)
    # Preprocess
    mark_aromatic_atoms!(mol)
    default_atom_charge!(mol)
    default_bond_order!(mol)
    # seems not to work yet for cdxml aromaticity yet
    kekulize!(mol)

    # Cache relatively expensive descriptors
    sssr!(mol)
    apparent_valence!(mol)
    valence!(mol)
    lone_pair!(mol)
    is_ring_aromatic!(mol)

    desc = mol[:descriptors]
    
    for i in vertices(mol)
        atom = mol.vprops[i]
        
        if atom.num_hs >= 0
            desc.valence[i] = desc.apparent_valence[i] + atom.num_hs
        end
    end
end

function ctab_atom_v2(::Type{CDXMLAtom}, line::AbstractString)
    xpos = parse(Float64, line[1:10])
    ypos = parse(Float64, line[11:20])
    zpos = parse(Float64, line[21:30])
    crds = Float64[xpos, ypos, zpos]
    sym = Symbol(rstrip(line[32:34]))
    # atom.mass_diff = parse(Int, line[35:36])
    iso = 0  # will be ignored, use ISO property instead
    old_sdf_charge = parse(Int, line[37:39])
    chg = SDF_CHARGE_TABLE[old_sdf_charge]
    mul = old_sdf_charge == 4 ? 2 : 1
    # atom.stereo_flag = parse(Int, line[40:42])
    # valence = parse(Int, line[46:48])
    isaromatic = false
    # -1 means automatic, any value >= 0 is the explicit number of hydrogens
    num_hs = parse(Int, line[43:45]) - 1
    return CDXMLAtom(sym, chg, mul, iso, isaromatic, crds, num_hs)
end

function ctab_atom_v3(::Type{CDXMLAtom}, line::AbstractString)
    ss = split(line)
    crds = parse.(Float64, ss[5:7])
    sym = ss[4]
    props = Dict(sympair.(ss[9:end])...)
    chg = get(props, :CHG, 0)
    iso = get(props, :MASS, 0)
    mul = get(props, :RAD, 1)
    isaromatic = false
    # -1 means automatic, any value >= 0 is the explicit number of hydrogens
    num_hs = parse(Int, line[43:45]) - 1
    return CDXMLAtom(sym, chg, mul, iso, isaromatic, crds, num_hs)
end

function ctab_bond_v2(
        ::Type{T}, ::Type{CDXMLBond}, line::AbstractString) where {T<:Integer}
    u = parse(T, line[1:3])
    v = parse(T, line[4:6])
    ord = parse(Int, line[7:9])
    notation = parse(Int, line[10:12])
    isordered = u < v
    u, v = isordered ? (u, v) : (v, u)
    return (u_edge(T, u, v), CDXMLBond(ord; notation, isordered))
end