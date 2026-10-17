#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function rdk_on_init!(mol::SimpleMolGraph)
    # Do nothing
end


function rdk_on_update!(mol::SimpleMolGraph)
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


# RDKit CommonChem JSON reading (concretely typed, see json.jl)

"""
    read_commonchem(::Type{G}, json::AbstractString) where G <: MolGraph -> G

Read a molecule from RDKit CommonChem JSON text.
"""
read_commonchem(::Type{G}, json::AbstractString
    ) where G <: MolGraph = read_commonchem(G, parse_jsonnode(json))

function read_commonchem(::Type{G}, data::JSONNode
        ) where {T,G<:MolGraph{T,CommonChemAtom,CommonChemBond}}
    jint(something(jget(data, "commonchem")), "version", 0) == 10 || error(
        "CommonChem version other than 10 is not supported")
    molecules = jitems(data, "molecules")
    length(molecules) == 1 || error("Only single molecule data is supported")
    moldata = molecules[1]
    extensions = jitems(moldata, "extensions")
    (!isempty(extensions) && jstring(extensions[1], "name", "") == "rdkitRepresentation"
        ) || error("Invalid RDKit CommonChem file format")
    jint(extensions[1], "formatVersion", 0) == 2 || error("Unsupported RDKit CommonChem version")
    gps = MolProperty{T}()
    # edges
    es = Edge{T}[]
    eps = Dict{EdgeKey{T},CommonChemBond}()
    for b in jitems(moldata, "bonds")
        atoms = jintvec(something(jget(b, "atoms")))
        edge = u_edge(T, atoms[1] + 1, atoms[2] + 1)
        push!(es, edge)
        stereo = jstring(b, "stereo", "unspecified")
        if stereo in ("cis", "trans")
            sa = jintvec(something(jget(b, "stereoAtoms")))
            gps.stereobond[edge] = Stereobond{T}(sa[1] + 1, sa[2] + 1, stereo == "cis")
        end
        eps[edge] = fromjson(CommonChemBond, b)
    end
    g = SimpleGraph(es)
    # vertices
    vps = Dict{VertexKey{T},CommonChemAtom}()
    for (i, a) in enumerate(jitems(moldata, "atoms"))
        vps[i] = fromjson(CommonChemAtom, a)
        stereo = jstring(a, "stereo", "unspecified")
        if stereo in ("cw", "ccw")
            nbrs = ordered_neighbors(g, i)
            # TODO: ambiguity in implicit H
            gps.stereocenter[i] = Stereocenter{T}(nbrs[1], nbrs[2], nbrs[3], stereo == "cw")
        end
    end
    # expand fadjlist for vprops of isolated nodes
    for _ in nv(g):(length(vps) - 1)
        push!(g.fadjlist, T[])
    end
    for cds in jitems(moldata, "conformers")
        # Wedges cannot be preserved. Use coordgen to generate 2D coords manually
        jint(cds, "dim", 2) == 3 || continue
        push!(gps.descriptors.coords3d, Coords3d(
            [Point3d(jfloat(jarray(p)[1]), jfloat(jarray(p)[2]), jfloat(jarray(p)[3]))
                for p in jitems(cds, "coords")]))
    end
    return G(reactive_molgraph(g, vps, eps, gps; on_init=rdk_on_init!, on_update=rdk_on_update!)...)
end

# `JSON.parse(json, T)` and `StructUtils.make(T, x)` (not compilable with `juliac --trim`)
StructUtils.lift(::Type{G}, x) where {G<:MolGraph{<:Integer,CommonChemAtom,CommonChemBond}
    } = read_commonchem(G, JSON.json(x))


"""
    to_rdkdict(mol::MolGraph) -> Dict{String,Any}

Convert the molecule object into `Dict` in RDKit compatible CommonChem JSON format.
"""
to_rdkdict(mol::SimpleMolGraph) = JSON.parse(to_rdkjson(mol); dicttype=Dict{String,Any})


"""
    to_rdkjson(mol::MolGraph) -> String

Convert the molecule object into RDKit compatible CommonChem JSON text.
"""
function to_rdkjson(mol::SimpleMolGraph)
    io = IOBuffer()
    write_rdkjson(io, mol)
    return String(take!(io))
end

# The JSON text is written directly with concrete types
# (serializing `Dict{String,Any}` is not compilable with `juliac --trim`).
function write_rdkjson(io::IO, mol::SimpleMolGraph)
    print(io, """{"commonchem":{"version":10},""")
    print(io, """"defaults":{"atom":{"z":6,"impHs":0,"chg":0,"nRad":0,"isotope":0,"stereo":"unspecified"},""")
    print(io, """"bond":{"bo":1,"stereo":"unspecified"}},""")  # `bo` is `type` in CommonChem specification
    print(io, """"molecules":[{"atoms":[""")
    atomnum = atom_number(mol)
    implh = implicit_hydrogens(mol)
    chg = atom_charge(mol)
    mul = multiplicity(mol)
    iso = isotope(mol)
    for i in vertices(mol)
        i == 1 || print(io, ",")
        print(io, "{")
        sep = ""  # separator between members
        if atomnum[i] != 6
            print(io, sep, "\"z\":", atomnum[i]); sep = ","
        end
        if implh[i] != 0
            print(io, sep, "\"impHs\":", implh[i]); sep = ","
        end
        if chg[i] != 0
            print(io, sep, "\"chg\":", chg[i]); sep = ","
        end
        if mul[i] != 1
            print(io, sep, "\"nRad\":", mul[i] - 1); sep = ","
        end
        if iso[i] != 0
            print(io, sep, "\"isotope\":", iso[i]); sep = ","
        end
        if haskey(mol[:stereocenter], i)
            center = mol[:stereocenter][i]
            nbrs = ordered_neighbors(mol, i)
            cw = isclockwise(center, nbrs[1], nbrs[2], nbrs[3])
            print(io, sep, "\"stereo\":\"", cw ? "cw" : "ccw", "\"")
        end
        print(io, "}")
    end
    print(io, "],\"bonds\":[")
    bondorder = bond_order(mol)
    for (i, e) in enumerate(edges(mol))
        i == 1 || print(io, ",")
        print(io, "{\"atoms\":[", src(e) - 1, ",", dst(e) - 1, "]")
        bondorder[i] == 1 || print(io, ",\"bo\":", bondorder[i])
        if haskey(mol[:stereobond], e)
            v1, v2, is_cis = mol[:stereobond][e]
            print(io, ",\"stereoAtoms\":[", v1 - 1, ",", v2 - 1, "]")
            print(io, ",\"stereo\":\"", is_cis ? "cis" : "trans", "\"")
        end
        print(io, "}")
    end
    print(io, "],")
    print(io, """"extensions":[{"name":"rdkitRepresentation","formatVersion":2,"toolkitVersion":"2022.09.5"}],""")
    # Wedges cannot be preserved. Use coordgen to generate 2D coords manually
    print(io, "\"conformers\":[")
    for (k, crds) in enumerate(mol[:descriptors].coords3d)
        k == 1 || print(io, ",")
        print(io, "{\"dim\":3,\"coords\":[")
        for (n, p) in enumerate(crds)
            n == 1 || print(io, ",")
            print(io, "[", p[1], ",", p[2], ",", p[3], "]")
        end
        print(io, "]}")
    end
    print(io, "]}]}")
    return
end


"""
    rdkitmol(mol::MolGraph) -> RDKitMinimalLib.Mol

Convert the molecule object into a RDKit molecule object that can be used in RDKitMinimalLib.jl
"""
function rdkitmol end


"""
    smiles(mol::MolGraph) -> String

Return a SMILES string of the molecule object.
"""
function smiles end


"""
    morgan_fp_vector(mol::MolGraph, details=nothing) -> BitArray

Return a Morgan fingerprint bit array
"""
function morgan_fp_vector end

"""
    morgan_fp_string(mol::MolGraph, details=nothing) -> String

Return the fingerprint of [`morgan_fp_vector`](@ref) as a string.
"""
function morgan_fp_string end


"""
    rdkit_fp_vector(mol::MolGraph, details=nothing) -> BitArray

Return a RDKit fingerprint bit array
"""
function rdkit_fp_vector end

"""
    rdkit_fp_string(mol::MolGraph, details=nothing) -> String

Return the fingerprint of [`rdkit_fp_vector`](@ref) as a string.
"""
function rdkit_fp_string end


"""
    pattern_fp_vector(mol::MolGraph, details=nothing) -> BitArray

Return a pattern fingerprint bit array, a topological fingerprint
optimized for substructure screening
"""
function pattern_fp_vector end

"""
    pattern_fp_string(mol::MolGraph, details=nothing) -> String

Return the fingerprint of [`pattern_fp_vector`](@ref) as a string.
"""
function pattern_fp_string end


"""
    atom_pair_fp_vector(mol::MolGraph, details=nothing) -> BitArray

Return a atom pairs fingerprint bit array
"""
function atom_pair_fp_vector end

"""
    atom_pair_fp_string(mol::MolGraph, details=nothing) -> String

Return the fingerprint of [`atom_pair_fp_vector`](@ref) as a string.
"""
function atom_pair_fp_string end


"""
    topological_torsion_fp_vector(mol::MolGraph, details=nothing) -> BitArray

Return a topological torsions fingerprint bit array
"""
function topological_torsion_fp_vector end

"""
    topological_torsion_fp_string(mol::MolGraph, details=nothing) -> String

Return the fingerprint of [`topological_torsion_fp_vector`](@ref) as a string.
"""
function topological_torsion_fp_string end



jaccard_index(a::BitVector, b::BitVector
    ) = sum(a .& b) / sum(a .| b)  # Not exported
