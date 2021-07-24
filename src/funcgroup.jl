#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    FunctionalGroupClassifier,
    functionalgroupgraph,
    functionalgroupmap


const FUNCGROUPTABLE = let
    fgroupfile = joinpath(dirname(@__FILE__), "../assets/const/functionalgroup.yaml")
    include_dependency(fgroupfile)
    YAML.load(open(fgroupfile))
end


struct FGTermNode <: AbstractNode
    term::Symbol
end


struct FGRelationEdge <: DirectedEdge
    relation::Symbol
end


struct FunctionalGroupClassifier <: OrderedDiGraph
    outneighbormap::Vector{Dict{Int,Int}}
    inneighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{FGTermNode}
    edgeattrs::Vector{FGRelationEdge}
    componentmap::Dict{Symbol,Set{Set{Int}}}

    FunctionalGroupClassifier() = new([], [], [], [], [], Dict())
end

FGC = FunctionalGroupClassifier


getterm(graph::FGC, term::Symbol) = graph.componentmap[term]
hasterm(graph::FGC, term::Symbol) = haskey(graph.componentmap, term)

function setterm!(graph::FGC, term::Symbol, component::Set{Set{Int}})
    graph.componentmap[term] = component
end


"""
    functionalgroupgraph(mol::GraphMol) -> FunctionalGroupClassifier

Generate functional group graph that is a directed acyclic graph similar to
an ontology graph.
"""
function functionalgroupgraph(mol::GraphMol; funcgrouptable=FUNCGROUPTABLE)
    fgc = FunctionalGroupClassifier()
    termidmap = Dict{Symbol,Int}()
    for rcd in funcgrouptable
        components = fgrouprecord(mol, fgc, rcd)
        isempty(components) && continue
        term = Symbol(rcd["key"])
        setterm!(fgc, term, components)
        # Update class graph
        termid = addnode!(fgc, FGTermNode(term))
        termidmap[term] = termid
        if haskey(rcd, "has")
            for k in rcd["has"]
                child = termidmap[Symbol(k)]
                addedge!(fgc, termid, child, FGRelationEdge(:has))
            end
        end
        if haskey(rcd, "isa")
            for k in rcd["isa"]
                child = termidmap[Symbol(k)]
                addedge!(fgc, termid, child, FGRelationEdge(:isa))
            end
        end
    end
    return fgc
end


function fgrouprecord(mol::GraphMol, fgc::FGC, rcd)
    newset = Set{Set{Int}}()
    # Membership filter
    if haskey(rcd, "has")
        for k in rcd["has"]
            if !hasterm(fgc, Symbol(k))
                return Set{Set{Int}}()
            end
        end
    end
    # Substructure match
    if haskey(rcd, "isa")
        q = parse(SMARTS, rcd["query"])
        allset = Set{Set{Int}}[]
        for k in rcd["isa"]
            hasterm(fgc, Symbol(k)) || continue
            refset = getterm(fgc, Symbol(k))
            eachset = Set{Set{Int}}()
            for s in refset
                subst = nodesubgraph(mol, s)
                if hasexactmatch(subst, q)
                    push!(eachset, s)
                end
            end
            push!(allset, eachset)
        end
        isempty(allset) && return Set{Set{Int}}()
        return intersect(allset...)
    end
    return fgroupcond(mol, rcd)
end


function fgroupcond(mol::GraphMol, rcd)
    newset = Set{Set{Int}}()
    if haskey(rcd, "any")
        for query in rcd["any"]
            union!(newset, fgroupquery(mol, query))
        end
    else
        union!(newset, fgroupquery(mol, rcd["query"]))
    end
    return newset
end


function fgroupquery(mol::GraphMol, query)
    q = parse(SMARTS, query)
    newset = Set{Set{Int}}()
    for nmap in substructmatches(mol, q, fastsingleton=true)
        isempty(nmap) || push!(newset, keys(nmap))
    end
    return newset
end


"""
    functionalgroupmap(graph::FGC) -> Dict{Symbol,Set{Set{Int}}}

Return mapping between terms of FGC leafs (roots of the DAG) and their corresponding atom indices.
"""
function functionalgroupmap(graph::FGC)
    mapping = Dict()
    for r in roots(graph)
        term = nodeattr(graph, r).term
        mapping[term] = graph.componentmap[term]
    end
    return mapping
end