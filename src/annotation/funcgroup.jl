#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    funcgrouptable,
    FUNC_GROUP_TABLE,
    FunctionalGroup,
    functionalgroup!,
    fgrouprecord,
    fgroupcond,
    fgroupquery

a = 1

function funcgrouptable()
    table = []
    files = [
        "funcgroup.yaml",
        "ring.yaml",
        "biomolecule.yaml"
    ]
    dir = joinpath(dirname(@__FILE__), "..", "..", "assets", "funcgroup")
    for f in files
        src = joinpath(dir, f)
        data = YAML.load(open(src))
        println("loading: $(f)")
        append!(table, data)
    end
    return table
end

FUNC_GROUP_TABLE = funcgrouptable()


struct FunctionalGroup <: Annotation
    nodeset::Dict{Symbol,Set{Set{Int}}}

    function FunctionalGroup()
        new(Dict())
    end
end


function functionalgroup!(mol::VectorMol)
    required_annotation(mol, :Topology)
    required_annotation(mol, :Elemental)
    required_annotation(mol, :Aromatic)
    mol.annotation[:FunctionalGroup] = FunctionalGroup()
    for rcd in FUNC_GROUP_TABLE
        fgset = fgrouprecord(mol, rcd)
        fgkey = Symbol(rcd["key"])
        mol.annotation[:FunctionalGroup].nodeset[fgkey] = fgset
    end
end


function fgrouprecord(mol::VectorMol, rcd)
    fgsetmap = mol.annotation[:FunctionalGroup].nodeset
    newset = Set{Set{Int}}()
    # Membership filter
    if "have" in keys(rcd)
        for k in rcd["have"]
            if isempty(fgsetmap[Symbol(k)])
                return newset
            end
        end
        # TODO: update ontology graph
    end
    # Substructure match
    preprocess!(mol)
    if "isa" in keys(rcd)
        q = parse(SMARTS, rcd["query"])
        refset = fgsetmap[Symbol(rcd["isa"])]
        for s in refset
            # TODO: SubstructureView for isomorph mapping
            subg = nodesubgraph(mol.graph, s)
            state = molidentstate(
                subg, q.graph, atommatch(mol, q), bondmatch(mol, q))
            mappings = substructmap!(mol, q, state)
            for (emap, nmap) in mappings
                esub = edgesubgraph(mol.graph, keys(emap))
                # nodes = union(nodekeys(esub), keys(nmap))
                # push!(newset, nodes)
                push!(newset, nodekeys(esub))
            end
        end
        # TODO: update ontology graph
    else
        union!(newset, fgroupcond(mol, rcd))
    end
    return newset
end


function fgroupcond(mol::VectorMol, rcd)
    newset = Set{Set{Int}}()
    if "any" in keys(rcd)
        for query in rcd["any"]
            union!(newset, fgroupquery(mol, query))
        end
    else
        union!(newset, fgroupquery(mol, rcd["query"]))
    end
    return newset
end


function fgroupquery(mol::VectorMol, query)
    q = parse(SMARTS, query)
    newset = Set{Set{Int}}()
    state = molsubstrstate(
        mol.graph, q.graph, atommatch(mol, q), bondmatch(mol, q))
    mappings = substructmap!(mol, q, state)
    for (emap, nmap) in mappings
        esub = edgesubgraph(mol.graph, keys(emap))
        nodes = union(nodekeys(esub), keys(nmap))
        push!(newset, nodes)
    end
    return newset
end
