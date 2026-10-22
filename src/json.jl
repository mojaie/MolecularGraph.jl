#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# NamedTuple (instead of Dict{String,Any}) so that JSON writing is compilable with `juliac --trim`.
# Keys are in alphabetical order, the same order as the former sorted Dict output.
function StructUtils.lower(x::ReactiveMolGraph{T,V,E}) where {T,V,E}
    dispatch_update!(x)  # update descriptors before serialization
    return (
        eprops = x.eprops,
        eproptype = string(nameof(E)),
        gprops = x.gprops,
        graph = x.graph,
        state = Dict{String,Int}(),  # state is not serialized
        vprops = x.vprops,
        vproptype = string(nameof(V))
    )
end

StructUtils.structlike(::Type{<:MolGraph}) = false
StructUtils.structlike(::Type{<:QueryMolGraph}) = false


# Molecule JSON writer and reader
#
# JSON.jl serializes and deserializes nested types by generic recursion, which is not compilable
# with `juliac --trim` for deeply nested molecule data. The following writer and reader handle the
# molecule JSON format with concrete types only. The output is identical to JSON.jl
# (`Dict` keys sorted, struct fields and `OrderedDict` in order, shortest float representation).


# Writer

function write_jsonstring(io::IO, s::AbstractString)
    # same escaping as JSON.jl
    write(io, UInt8('"'))
    for b in codeunits(s)
        if b == UInt8('"') || b == UInt8('\\')
            write(io, UInt8('\\'), b)
        elseif b == UInt8('\b')
            write(io, UInt8('\\'), UInt8('b'))
        elseif b == UInt8('\f')
            write(io, UInt8('\\'), UInt8('f'))
        elseif b == UInt8('\n')
            write(io, UInt8('\\'), UInt8('n'))
        elseif b == UInt8('\r')
            write(io, UInt8('\\'), UInt8('r'))
        elseif b == UInt8('\t')
            write(io, UInt8('\\'), UInt8('t'))
        elseif b < 0x20 || b == 0x7f
            print(io, "\\u", string(b, base=16, pad=4))
        else
            write(io, b)
        end
    end
    write(io, UInt8('"'))
    return
end

write_jsonvalue(io::IO, x::AbstractString) = write_jsonstring(io, x)
write_jsonvalue(io::IO, x::Symbol) = write_jsonstring(io, string(x))
write_jsonvalue(io::IO, x::Bool) = print(io, x ? "true" : "false")
write_jsonvalue(io::IO, x::Integer) = print(io, x)

function write_jsonvalue(io::IO, x::Float64)
    isfinite(x) || throw(ArgumentError("$(x) not allowed to be written in JSON spec"))
    print(io, x)  # shortest representation (Ryu), same as JSON.jl
    return
end

function write_jsonvalue(io::IO, x::Union{AbstractVector,Tuple})
    print(io, "[")
    for (i, v) in enumerate(x)
        i == 1 || print(io, ",")
        write_jsonvalue(io, v)
    end
    print(io, "]")
    return
end

# object members: (key::String, value) pairs in output order
function write_jsonobject(io::IO, members::AbstractVector{<:Tuple{String,Any}})
    print(io, "{")
    for (i, (k, v)) in enumerate(members)
        i == 1 || print(io, ",")
        write_jsonstring(io, k)
        print(io, ":")
        write_jsonvalue(io, v)
    end
    print(io, "}")
    return
end

# Members of different value types (e.g. struct fields) are given as a Tuple.
# Unrolled at compile time so that each value is written by a statically dispatched method.
@generated function write_jsonobject(io::IO, members::Tuple)
    ex = Expr(:block, :(print(io, "{")))
    for i in 1:length(members.parameters)
        i == 1 || push!(ex.args, :(print(io, ",")))
        push!(ex.args, :(write_jsonstring(io, members[$i][1])))
        push!(ex.args, :(print(io, ":")))
        push!(ex.args, :(write_jsonvalue(io, members[$i][2])))
    end
    push!(ex.args, :(print(io, "}")))
    push!(ex.args, :(return))
    return ex
end

# Dict keys are sorted as JSON.jl does, OrderedDict keeps the order
jsonkey(k::VertexKey) = string(k.key)
jsonkey(k::EdgeKey) = string(k.key.src, "_", k.key.dst)
jsonkey(k::AbstractString) = String(k)
function json_members(d::AbstractDict{K,V}) where {K,V}
    members = Tuple{String,V}[]
    for (k, v) in d
        push!(members, (jsonkey(k), v))
    end
    return members
end
write_jsonvalue(io::IO, d::Dict) = write_jsonobject(io, sort!(json_members(d), by=first))
write_jsonvalue(io::IO, d::OrderedDict) = write_jsonobject(io, json_members(d))

# atom and bond properties (sparse, same as `StructUtils.lower`)
# Written by dedicated functions that do not call `write_jsonvalue` recursively
# (recursion with various argument types is not compilable with `juliac --trim`).
write_jsonprop(io::IO, x::AbstractString) = write_jsonstring(io, x)
write_jsonprop(io::IO, x::Bool) = print(io, x ? "true" : "false")
write_jsonprop(io::IO, x::Int) = print(io, x)

function write_jsonprop(io::IO, x::Vector{Float64})
    print(io, "[")
    for (i, v) in enumerate(x)
        i == 1 || print(io, ",")
        isfinite(v) || throw(ArgumentError("$(v) not allowed to be written in JSON spec"))
        print(io, v)
    end
    print(io, "]")
    return
end

function write_jsonprops(io::IO, d::Dict{String,V}) where V
    print(io, "{")
    for (i, k) in enumerate(sort!(collect(keys(d))))
        i == 1 || print(io, ",")
        write_jsonstring(io, k)
        print(io, ":")
        write_jsonprop(io, d[k])
    end
    print(io, "}")
    return
end

write_jsonvalue(io::IO, x::Union{SDFAtom,SMILESAtom,CommonChemAtom,SDFBond,SMILESBond,CommonChemBond,QueryNode}
    ) = write_jsonprops(io, StructUtils.lower(x))

write_jsonvalue(io::IO, x::Stereocenter) = write_jsonvalue(io, Tuple(x))
write_jsonvalue(io::IO, x::Stereobond) = write_jsonvalue(io, Tuple(x))
write_jsonvalue(io::IO, x::Union{StereocenterMap,StereobondMap}
    ) = write_jsonobject(io, (("mapping", x.mapping),))
write_jsonvalue(io::IO, x::Coords2d) = write_jsonvalue(io, StructUtils.lower(x))
write_jsonvalue(io::IO, x::Coords3d) = write_jsonvalue(io, StructUtils.lower(x))
write_jsonvalue(io::IO, x::Draw2dBondStyle) = write_jsonvalue(io, StructUtils.lower(x))

write_jsonvalue(io::IO, g::SimpleGraph) = write_jsonobject(io, (("ne", g.ne), ("fadjlist", g.fadjlist)))
write_jsonvalue(io::IO, g::SimpleDiGraph) = write_jsonobject(
    io, (("ne", g.ne), ("fadjlist", g.fadjlist), ("badjlist", g.badjlist)))
write_jsonvalue(io::IO, x::QueryTree) = write_jsonobject(io, (("graph", x.graph), ("vprops", x.vprops)))

write_jsonvalue(io::IO, x::MolDescriptor) = write_jsonobject(io, (
    ("atom_charge", x.atom_charge), ("bond_order", x.bond_order), ("sssr", x.sssr),
    ("lone_pair", x.lone_pair), ("apparent_valence", x.apparent_valence), ("valence", x.valence),
    ("is_ring_aromatic", x.is_ring_aromatic), ("coords2d", x.coords2d), ("coords3d", x.coords3d),
    ("draw2d_bond_style", x.draw2d_bond_style)))
write_jsonvalue(io::IO, x::QueryMolDescriptor) = write_jsonobject(io, (
    ("coords2d", x.coords2d), ("coords3d", x.coords3d)))

write_jsonvalue(io::IO, x::MolProperty) = write_jsonobject(io, (
    ("stereocenter", x.stereocenter), ("stereobond", x.stereobond), ("pyrrole_like", x.pyrrole_like),
    ("smarts_input", x.smarts_input), ("smarts_lexical_succ", x.smarts_lexical_succ),
    ("descriptors", x.descriptors), ("metadata", x.metadata), ("logs", x.logs)))
write_jsonvalue(io::IO, x::QueryMolProperty) = write_jsonobject(io, (
    ("stereocenter", x.stereocenter), ("stereobond", x.stereobond),
    ("smarts_input", x.smarts_input), ("smarts_lexical_succ", x.smarts_lexical_succ),
    ("smarts_connectivity", x.smarts_connectivity), ("descriptors", x.descriptors),
    ("metadata", x.metadata), ("logs", x.logs)))

# molecule types supported by the writer and the reader
const JSONMolGraph = Union{
    MolGraph{<:Integer,SDFAtom,SDFBond}, MolGraph{<:Integer,SMILESAtom,SMILESBond},
    MolGraph{<:Integer,CommonChemAtom,CommonChemBond}, QueryMolGraph{<:Integer,<:QueryAtom,<:QueryBond}}

function write_jsonvalue(io::IO, mol::JSONMolGraph)
    V, E = vproptype(typeof(mol)), eproptype(typeof(mol))
    dispatch_update!(mol)  # update descriptors before serialization
    print(io, "{")
    write_jsonstring(io, "eprops"); print(io, ":"); write_jsonvalue(io, mol.eprops)
    print(io, ",\"eproptype\":"); write_jsonstring(io, string(nameof(E)))
    print(io, ",\"gprops\":"); write_jsonvalue(io, mol.gprops)
    print(io, ",\"graph\":"); write_jsonvalue(io, mol.graph)
    print(io, ",\"state\":{}")  # state is not serialized
    print(io, ",\"vprops\":"); write_jsonvalue(io, mol.vprops)
    print(io, ",\"vproptype\":"); write_jsonstring(io, string(nameof(V)))
    print(io, "}")
    return
end

"""
    JSON.json(mol::ReactiveMolGraph) -> String

Serialize the molecule into JSON text. Molecules of SDFile, SMILES, SMARTS and CommonChem types
are written by the concretely typed writer (compilable with `juliac --trim`).
"""
function JSON.json(mol::JSONMolGraph; kwargs...)
    isempty(kwargs) || return invoke(JSON.json, Tuple{Any}, mol; kwargs...)  # e.g. pretty=true
    io = IOBuffer()
    write_jsonvalue(io, mol)
    return String(take!(io))
end


# Reader

"""
    JSONNode

Concretely typed JSON value tree used by the molecule JSON reader.

`kind`: :null, :bool, :number, :string, :array or :object. `str` holds the string value,
the number text or "true"/"false". `items` holds array items or object values,
and `keys` holds object keys.
"""
struct JSONNode
    kind::Symbol
    str::String
    keys::Vector{String}
    items::Vector{JSONNode}
end

mutable struct JSONCursor
    data::Vector{UInt8}
    pos::Int
end

json_error(c::JSONCursor, msg) = throw(ArgumentError("invalid JSON at byte $(c.pos): $(msg)"))

function skip_ws!(c::JSONCursor)
    while c.pos <= length(c.data) && c.data[c.pos] in (UInt8(' '), UInt8('\t'), UInt8('\n'), UInt8('\r'))
        c.pos += 1
    end
    return
end

function peek_byte(c::JSONCursor)
    skip_ws!(c)
    c.pos <= length(c.data) || json_error(c, "unexpected end")
    return c.data[c.pos]
end

function expect_byte!(c::JSONCursor, b::UInt8)
    peek_byte(c) == b || json_error(c, "expected '$(Char(b))'")
    c.pos += 1
    return
end

function parse_hex4!(c::JSONCursor)
    c.pos + 3 <= length(c.data) || json_error(c, "invalid unicode escape")
    v = parse(UInt16, String(c.data[c.pos:c.pos + 3]); base=16)
    c.pos += 4
    return v
end

function parse_jsonstring!(c::JSONCursor)
    expect_byte!(c, UInt8('"'))
    io = IOBuffer()
    while true
        c.pos <= length(c.data) || json_error(c, "unterminated string")
        b = c.data[c.pos]
        c.pos += 1
        b == UInt8('"') && break
        if b != UInt8('\\')
            write(io, b)
            continue
        end
        c.pos <= length(c.data) || json_error(c, "unterminated string")
        e = c.data[c.pos]
        c.pos += 1
        if e == UInt8('"') || e == UInt8('\\') || e == UInt8('/')
            write(io, e)
        elseif e == UInt8('b')
            write(io, UInt8('\b'))
        elseif e == UInt8('f')
            write(io, UInt8('\f'))
        elseif e == UInt8('n')
            write(io, UInt8('\n'))
        elseif e == UInt8('r')
            write(io, UInt8('\r'))
        elseif e == UInt8('t')
            write(io, UInt8('\t'))
        elseif e == UInt8('u')
            u = parse_hex4!(c)
            if 0xd800 <= u <= 0xdbff  # surrogate pair
                (c.pos + 1 <= length(c.data) && c.data[c.pos] == UInt8('\\') && c.data[c.pos + 1] == UInt8('u')
                    ) || json_error(c, "invalid surrogate pair")
                c.pos += 2
                l = parse_hex4!(c)
                print(io, Char(0x10000 + ((UInt32(u) - 0xd800) << 10) + (UInt32(l) - 0xdc00)))
            else
                print(io, Char(u))
            end
        else
            json_error(c, "invalid escape")
        end
    end
    return String(take!(io))
end

function parse_jsonnode!(c::JSONCursor)
    b = peek_byte(c)
    if b == UInt8('{')
        c.pos += 1
        ks = String[]
        vs = JSONNode[]
        if peek_byte(c) == UInt8('}')
            c.pos += 1
            return JSONNode(:object, "", ks, vs)
        end
        while true
            push!(ks, parse_jsonstring!(c))
            expect_byte!(c, UInt8(':'))
            push!(vs, parse_jsonnode!(c))
            nb = peek_byte(c)
            c.pos += 1
            nb == UInt8('}') && break
            nb == UInt8(',') || json_error(c, "expected ',' or '}'")
        end
        return JSONNode(:object, "", ks, vs)
    elseif b == UInt8('[')
        c.pos += 1
        vs = JSONNode[]
        if peek_byte(c) == UInt8(']')
            c.pos += 1
            return JSONNode(:array, "", String[], vs)
        end
        while true
            push!(vs, parse_jsonnode!(c))
            nb = peek_byte(c)
            c.pos += 1
            nb == UInt8(']') && break
            nb == UInt8(',') || json_error(c, "expected ',' or ']'")
        end
        return JSONNode(:array, "", String[], vs)
    elseif b == UInt8('"')
        return JSONNode(:string, parse_jsonstring!(c), String[], JSONNode[])
    end
    # literal or number
    start = c.pos
    while c.pos <= length(c.data) && !(c.data[c.pos] in (
            UInt8(','), UInt8('}'), UInt8(']'), UInt8(' '), UInt8('\t'), UInt8('\n'), UInt8('\r')))
        c.pos += 1
    end
    token = String(c.data[start:c.pos - 1])
    token == "null" && return JSONNode(:null, "", String[], JSONNode[])
    token in ("true", "false") && return JSONNode(:bool, token, String[], JSONNode[])
    isnothing(tryparse(Float64, token)) && json_error(c, "invalid value")
    return JSONNode(:number, token, String[], JSONNode[])
end

function parse_jsonnode(json::AbstractString)
    c = JSONCursor(Vector{UInt8}(codeunits(String(json))), 1)
    node = parse_jsonnode!(c)
    skip_ws!(c)
    c.pos > length(c.data) || json_error(c, "unexpected trailing characters")
    return node
end

# accessors

function jget(node::JSONNode, key::String)
    node.kind === :object || return nothing
    i = findfirst(==(key), node.keys)
    return isnothing(i) ? nothing : node.items[i]
end

jmembers(node::JSONNode) = zip(node.keys, node.items)  # object members

function jint(node::JSONNode)
    node.kind === :number || throw(ArgumentError("JSON number expected"))
    return parse(Int, node.str)
end

function jfloat(node::JSONNode)
    node.kind === :number || throw(ArgumentError("JSON number expected"))
    return parse(Float64, node.str)
end

function jbool(node::JSONNode)
    node.kind === :bool || throw(ArgumentError("JSON boolean expected"))
    return node.str == "true"
end

function jstring(node::JSONNode)
    node.kind === :string || throw(ArgumentError("JSON string expected"))
    return node.str
end

jarray(node::JSONNode) = node.kind === :array ? node.items : throw(ArgumentError("JSON array expected"))

# with default values (atom and bond properties are written sparsely)
function jint(node::JSONNode, key::String, default::Int)
    v = jget(node, key)
    return isnothing(v) ? default : jint(v)
end
function jbool(node::JSONNode, key::String, default::Bool)
    v = jget(node, key)
    return isnothing(v) ? default : jbool(v)
end
function jstring(node::JSONNode, key::String, default::String)
    v = jget(node, key)
    return isnothing(v) ? default : jstring(v)
end
function jitems(node::JSONNode, key::String)
    v = jget(node, key)
    return isnothing(v) ? JSONNode[] : jarray(v)
end

jintvec(node::JSONNode) = Int[jint(v) for v in jarray(node)]
jintvecs(nodes::Vector{JSONNode}) = Vector{Int}[jintvec(v) for v in nodes]
jfloatvec(node::JSONNode) = Float64[jfloat(v) for v in jarray(node)]

# conversion

fromjson(::Type{SDFAtom}, x::JSONNode) = SDFAtom(
    jstring(x, "symbol", "C"), jint(x, "charge", 0), jint(x, "multiplicity", 1),
    jint(x, "isotope", 0), isnothing(jget(x, "coords")) ? nothing : jfloatvec(something(jget(x, "coords"))),
    jint(x, "num_hs", -1))
fromjson(::Type{SMILESAtom}, x::JSONNode) = SMILESAtom(
    jstring(x, "symbol", "C"), jint(x, "charge", 0), jint(x, "multiplicity", 1),
    jint(x, "isotope", 0), jbool(x, "isaromatic", false), jstring(x, "stereo", "unspecified"))
fromjson(::Type{CommonChemAtom}, x::JSONNode) = CommonChemAtom(
    jint(x, "z", 6), jint(x, "chg", 0), jint(x, "impHs", 0), jint(x, "isotope", 0),
    jint(x, "nRad", 0), jstring(x, "stereo", "unspecified"))
fromjson(::Type{SDFBond}, x::JSONNode) = SDFBond(
    jint(x, "order", 1), jint(x, "notation", 0), jbool(x, "isordered", true))
fromjson(::Type{SMILESBond}, x::JSONNode) = SMILESBond(
    jint(x, "order", 1), jbool(x, "isaromatic", false), jstring(x, "direction", "unspecified"))
# `bo` in RDKit notation, `type` in CommonChem specification
fromjson(::Type{CommonChemBond}, x::JSONNode) = CommonChemBond(jint(x, "bo", jint(x, "type", 1)))
fromjson(::Type{QueryNode}, x::JSONNode) = QueryNode(
    Symbol(jstring(x, "operator", "eq")), Symbol(jstring(x, "key", "none")), jstring(x, "value", ""))

function fromjson(::Type{Q}, x::JSONNode) where {T,U,Q<:QueryTree{T,U}}
    gnode = something(jget(x, "graph"))
    g = SimpleDiGraph{T}(jint(gnode, "ne", 0),
        Vector{T}.(jintvecs(jitems(gnode, "fadjlist"))), Vector{T}.(jintvecs(jitems(gnode, "badjlist"))))
    vprops = Dict{VertexKey{T},U}()
    vnode = jget(x, "vprops")
    if !isnothing(vnode)
        for (k, v) in jmembers(vnode)
            vprops[VertexKey{T}(parse(T, k))] = fromjson(U, v)
        end
    end
    return Q(g, vprops)
end

function edgekey(::Type{T}, s::AbstractString) where T
    i = something(findfirst('_', s))
    return EdgeKey{T}(u_edge(T, parse(T, s[1:i - 1]), parse(T, s[i + 1:end])))
end

function stereocenters_fromjson(::Type{T}, x::Union{JSONNode,Nothing}) where T
    m = StereocenterMap{T}()
    mapping = isnothing(x) ? nothing : jget(x, "mapping")
    isnothing(mapping) && return m
    for (k, v) in jmembers(mapping)
        a = jarray(v)
        m.mapping[VertexKey{T}(parse(T, k))] = Stereocenter{T}(
            jint(a[1]), jint(a[2]), jint(a[3]), jbool(a[4]))
    end
    return m
end

function stereobonds_fromjson(::Type{T}, x::Union{JSONNode,Nothing}) where T
    m = StereobondMap{T}()
    mapping = isnothing(x) ? nothing : jget(x, "mapping")
    isnothing(mapping) && return m
    for (k, v) in jmembers(mapping)
        a = jarray(v)
        m.mapping[edgekey(T, k)] = Stereobond{T}(jint(a[1]), jint(a[2]), jbool(a[3]))
    end
    return m
end

coords2d_fromjson(nodes::Vector{JSONNode}) = Coords2d[
    Coords2d([Point2d(jfloat(jarray(p)[1]), jfloat(jarray(p)[2])) for p in jarray(c)]) for c in nodes]
coords3d_fromjson(nodes::Vector{JSONNode}) = Coords3d[
    Coords3d([Point3d(jfloat(jarray(p)[1]), jfloat(jarray(p)[2]), jfloat(jarray(p)[3])) for p in jarray(c)])
    for c in nodes]

function strdict_fromjson(::Type{D}, x::Union{JSONNode,Nothing}) where D <: AbstractDict
    d = D()
    isnothing(x) && return d
    for (k, v) in jmembers(x)
        d[k] = jstring(v)
    end
    return d
end

function fromjson(::Type{MolDescriptor{T}}, x::Union{JSONNode,Nothing}) where T
    isnothing(x) && return MolDescriptor{T}()
    return MolDescriptor{T}(
        atom_charge=Int[jint(v) for v in jitems(x, "atom_charge")],
        bond_order=Int[jint(v) for v in jitems(x, "bond_order")],
        sssr=Vector{T}.(jintvecs(jitems(x, "sssr"))),
        lone_pair=Int[jint(v) for v in jitems(x, "lone_pair")],
        apparent_valence=Int[jint(v) for v in jitems(x, "apparent_valence")],
        valence=Int[jint(v) for v in jitems(x, "valence")],
        is_ring_aromatic=Bool[jbool(v) for v in jitems(x, "is_ring_aromatic")],
        coords2d=coords2d_fromjson(jitems(x, "coords2d")),
        coords3d=coords3d_fromjson(jitems(x, "coords3d")),
        draw2d_bond_style=Draw2dBondStyle[
            Draw2dBondStyle(Symbol[Symbol(jstring(s)) for s in jarray(v)])
            for v in jitems(x, "draw2d_bond_style")])
end

function fromjson(::Type{QueryMolDescriptor{T}}, x::Union{JSONNode,Nothing}) where T
    isnothing(x) && return QueryMolDescriptor{T}()
    return QueryMolDescriptor{T}(
        coords2d=coords2d_fromjson(jitems(x, "coords2d")),
        coords3d=coords3d_fromjson(jitems(x, "coords3d")))
end

fromjson(::Type{MolProperty{T}}, x::JSONNode) where T = MolProperty{T}(
    stereocenter=stereocenters_fromjson(T, jget(x, "stereocenter")),
    stereobond=stereobonds_fromjson(T, jget(x, "stereobond")),
    pyrrole_like=T[jint(v) for v in jitems(x, "pyrrole_like")],
    smarts_input=jstring(x, "smarts_input", ""),
    smarts_lexical_succ=Vector{T}.(jintvecs(jitems(x, "smarts_lexical_succ"))),
    descriptors=fromjson(MolDescriptor{T}, jget(x, "descriptors")),
    metadata=strdict_fromjson(OrderedDict{String,String}, jget(x, "metadata")),
    logs=strdict_fromjson(Dict{String,String}, jget(x, "logs")))

fromjson(::Type{QueryMolProperty{T}}, x::JSONNode) where T = QueryMolProperty{T}(
    stereocenter=stereocenters_fromjson(T, jget(x, "stereocenter")),
    stereobond=stereobonds_fromjson(T, jget(x, "stereobond")),
    smarts_input=jstring(x, "smarts_input", ""),
    smarts_lexical_succ=Vector{T}.(jintvecs(jitems(x, "smarts_lexical_succ"))),
    smarts_connectivity=Vector{T}.(jintvecs(jitems(x, "smarts_connectivity"))),
    descriptors=fromjson(QueryMolDescriptor{T}, jget(x, "descriptors")),
    metadata=strdict_fromjson(OrderedDict{String,String}, jget(x, "metadata")),
    logs=strdict_fromjson(Dict{String,String}, jget(x, "logs")))


gproptype(::Type{<:MolGraph{T}}) where T = MolProperty{T}
gproptype(::Type{<:QueryMolGraph{T}}) where T = QueryMolProperty{T}

json_callbacks(::Type{<:MolGraph{T,SDFAtom,SDFBond}}) where T = (on_init=sdf_on_init!, on_update=sdf_on_update!)
json_callbacks(::Type{<:MolGraph{T,SMILESAtom,SMILESBond}}) where T = (on_init=smiles_on_init!, on_update=smiles_on_update!)
# Same callbacks as `smartstomol` so that the molecule has the same type
json_callbacks(::Type{<:QueryMolGraph}) = (on_init=smarts_on_init!, on_update=default_on_update!)
json_callbacks(::Type{<:ReactiveMolGraph}) = (on_init=default_on_init!, on_update=default_on_update!)

function read_moljson(::Type{G}, data::JSONNode; kwargs...
        ) where {T,V,E,G<:ReactiveMolGraph{T,V,E}}
    gnode = something(jget(data, "graph"))
    g = SimpleGraph{T}(jint(gnode, "ne", 0), Vector{T}.(jintvecs(jitems(gnode, "fadjlist"))))
    vprops = Dict{VertexKey{T},V}()
    for (k, v) in jmembers(something(jget(data, "vprops")))
        vprops[VertexKey{T}(parse(T, k))] = fromjson(V, v)
    end
    eprops = Dict{EdgeKey{T},E}()
    for (k, v) in jmembers(something(jget(data, "eprops")))
        eprops[edgekey(T, k)] = fromjson(E, v)
    end
    gprops = fromjson(gproptype(G), something(jget(data, "gprops")))
    return G(reactive_molgraph(g, vprops, eprops, gprops;
        initialized=true,  # Skip initialization
        has_updates=false,  # Do not update ready-to-use descriptors!
        json_callbacks(G)..., kwargs...)...)
end

"""
    read_moljson(::Type{G}, json::AbstractString; kwargs...) where G <: ReactiveMolGraph -> G

Read a molecule of the given type from the JSON text written by `JSON.json(mol)`.
"""
read_moljson(::Type{G}, json::AbstractString; kwargs...
    ) where G <: ReactiveMolGraph = read_moljson(G, parse_jsonnode(json); kwargs...)


# `JSON.parse(json, T)` and `StructUtils.make(T, x)` (not compilable with `juliac --trim`)
StructUtils.lift(::Type{G}, x) where {G<:Union{MolGraph{<:Integer,SDFAtom,SDFBond},
    MolGraph{<:Integer,SMILESAtom,SMILESBond}, QueryMolGraph}} = read_moljson(G, JSON.json(x))


MolGraph{T,V,E}(json::AbstractString; on_init=nothing, on_update=nothing) where {T,V,E} = with_callbacks(
    read_moljson(MolGraph{T,V,E}, json); on_init=on_init, on_update=on_update)

QueryMolGraph{T,V,E}(json::AbstractString; on_init=nothing, on_update=nothing) where {T,V,E} = with_callbacks(
    read_moljson(QueryMolGraph{T,V,E}, json); on_init=on_init, on_update=on_update)


# JSON auto detect

"""
    mol_from_json(json::AbstractString; on_init=nothing, on_update=nothing) -> ReactiveMolGraph
    mol_from_json(f, json::AbstractString; on_init=nothing, on_update=nothing)

Read a molecule from the JSON text (`JSON.json(mol)` output or RDKit CommonChem JSON).
The molecule type is detected from the JSON.

The second form calls `f(mol)` and returns its result. `mol` has a concrete type in each
branch of the type detection, so this form is compilable with `juliac --trim`
(the return type of the first form cannot be inferred).

```julia
n = mol_from_json(json) do mol
    nv(mol)
end
```
"""
function mol_from_json(f::F, json::AbstractString; on_init=nothing, on_update=nothing) where F
    data = parse_jsonnode(json)
    if !isnothing(jget(data, "commonchem"))
        return f(with_callbacks(read_commonchem(MolGraph{Int,CommonChemAtom,CommonChemBond}, data);
            on_init=on_init, on_update=on_update))
    end
    vtype = jstring(data, "vproptype", "")
    etype = jstring(data, "eproptype", "")
    if vtype == "SDFAtom" && etype == "SDFBond"
        return f(with_callbacks(read_moljson(MolGraph{Int,SDFAtom,SDFBond}, data);
            on_init=on_init, on_update=on_update))
    elseif vtype == "SMILESAtom" && etype == "SMILESBond"
        return f(with_callbacks(read_moljson(MolGraph{Int,SMILESAtom,SMILESBond}, data);
            on_init=on_init, on_update=on_update))
    elseif vtype == "QueryAtom" && etype == "QueryBond"
        return f(with_callbacks(read_moljson(
            QueryMolGraph{Int,QueryAtom{Int,QueryNode},QueryBond{Int,QueryNode}}, data);
            on_init=on_init, on_update=on_update))
    end
    error("Invalid JSON format")
end

mol_from_json(json::AbstractString; kwargs...) = mol_from_json(identity, json; kwargs...)
