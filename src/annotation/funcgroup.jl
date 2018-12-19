#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    funcgrouptable,
    FUNC_GROUP_TABLE,
    FGTermNode,
    FGRelationEdge,
    FunctionalGroup,
    functionalgroup!,
    fgrouprecord,
    fgroupcond,
    fgroupquery,
    largestcomponents


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


struct FGTermNode <: AbstractNode
    term::Symbol
end


struct FGRelationEdge <: AbstractDirectedEdge
    source::Int
    target::Int
    relation::Symbol
end


struct FunctionalGroup <: Annotation
    nodeset::Dict{Symbol,Set{Set{Int}}}
    graph::DAG{FGTermNode,FGRelationEdge}

    function FunctionalGroup()
        new(Dict(), DAG{FGTermNode,FGRelationEdge}())
    end
end


function functionalgroup!(mol::VectorMol)
    required_annotation(mol, :Topology)
    required_annotation(mol, :Elemental)
    required_annotation(mol, :Aromatic)
    mol.annotation[:FunctionalGroup] = fg = FunctionalGroup()
    fggraphidx = Dict{Symbol,Int}()
    ncnt = 0
    ecnt = 0
    for rcd in FUNC_GROUP_TABLE
        fgset = fgrouprecord(mol, rcd)
        fgkey = Symbol(rcd["key"])
        fg.nodeset[fgkey] = fgset
        # Update ontology graph
        if !isempty(fgset)
            ncnt += 1
            fggraphidx[fgkey] = ncnt
            updatenode!(fg.graph, FGTermNode(fgkey), ncnt)
            if "have" in keys(rcd)
                for k in rcd["have"]
                    ecnt += 1
                    e = FGRelationEdge(fggraphidx[Symbol(k)], ncnt, :partof)
                    updateedge!(fg.graph, e, ecnt)
                end
            end
            if "isa" in keys(rcd)
                ecnt += 1
                e = FGRelationEdge(ncnt, fggraphidx[Symbol(rcd["isa"])], :isa)
                updateedge!(fg.graph, e, ecnt)
            end
        end
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



function largestcomponents(fg::FunctionalGroup)
    components = Dict{Symbol,Set{Set{Int}}}()
    ontnodes = topologicalsort(fg.graph)
    for n in ontnodes
        rmset = Set{Set{Int}}()
        nterm = getnode(fg.graph, n).term
        for (s, e) in successors(fg.graph, n)
            rel = getedge(fg.graph, e)
            if rel.relation != :partof
                continue
            end
            ansterm = getnode(fg.graph, s).term
            ansset = union(fg.nodeset[ansterm]...)
            for nset in fg.nodeset[nterm]
                if isempty(setdiff(nset, ansset))
                    push!(rmset, nset)
                end
            end
        end
        components[nterm] = setdiff(fg.nodeset[nterm], rmset)
    end
    for n in ontnodes
        rmset = Set{Set{Int}}()
        nterm = getnode(fg.graph, n).term
        for (s, e) in successors(fg.graph, n)
            rel = getedge(fg.graph, e)
            if rel.relation != :isa
                continue
            end
            ansterm = getnode(fg.graph, s).term
            intersect!(components[nterm], components[ansterm])
        end
        for (p, e) in predecessors(fg.graph, n)
            rel = getedge(fg.graph, e)
            if rel.relation != :isa
                continue
            end
            descterm = getnode(fg.graph, p).term
            setdiff!(components[nterm], components[descterm])
        end
    end
    return components
end
