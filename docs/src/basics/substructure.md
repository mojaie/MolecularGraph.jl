# Substructure and query

This page covers the following database search and substructure mining tasks:

- Substructure match
- InChI and InChIKey
- SMARTS query
- Structural alerts (e.g. PAINS)
- Functional group analysis
- Query containment

```@example substructure
using MolecularGraph, Graphs
```

The examples use molfiles bundled in the documentation sources. The same files can be
obtained from [PubChem](https://pubchem.ncbi.nlm.nih.gov/) (e.g. the SDF download of
Thiostrepton, PubChem SID 175266353, which has 2D coordinates).

```@example substructure
datadir = joinpath(pkgdir(MolecularGraph), "docs", "src", "assets", "data")
thiostrepton = sdftomol(joinpath(datadir, "Thiostrepton.mol"))
```


## Substructure match

- [`has_exact_match(mol1, mol2)`](@ref has_exact_match) returns whether the structures of
  `mol1` and `mol2` match exactly.
- [`has_substruct_match(mol1, mol2)`](@ref has_substruct_match) returns whether `mol2` is a
  substructure of `mol1`.

```@example substructure
has_exact_match(thiostrepton, thiostrepton)
```

```@example substructure
thiazole = smilestomol("s1cncc1")
```

```@example substructure
has_exact_match(thiostrepton, thiazole)
```

```@example substructure
has_substruct_match(thiostrepton, thiazole)
```

[`exact_matches(mol1, mol2)`](@ref exact_matches) and
[`substruct_matches(mol1, mol2)`](@ref substruct_matches) return lazy iterators that yield
vertex mappings (`Dict` of `mol1` vertex => `mol2` vertex) between the molecules.

```@example substructure
matches = collect(substruct_matches(thiostrepton, thiazole))
length(matches)
```

```@example substructure
first(matches)
```

The matched atoms can be highlighted in the structure drawing. `induced_subgraph_edges`
returns the bonds between the given atoms.

```@example substructure
function highlight_matches(mol, query; width=300, height=300)
    nodes = Set{Int}()
    for mapping in substruct_matches(mol, query)
        union!(nodes, keys(mapping))
    end
    nodes = collect(nodes)
    edges = induced_subgraph_edges(mol.graph, nodes)
    return html_fixed_size(mol, width, height, atomhighlight=nodes, bondhighlight=edges)
end

highlight_matches(thiostrepton, thiazole)
```

The substructure match functions accept the following keyword options (see
[`substruct_matches`](@ref) for details):

- `vmatchgen`, `ematchgen`: custom generators of atom and bond attribute matching functions
- `stereo`: if `true`, stereochemistry is also considered (default: `false`)
- `mandatory`, `forbidden`: node mappings that must (not) be included in the result
- `timeout`: abort the calculation after the given seconds (default: `nothing`, no timeout)

[`node_substruct_matches`](@ref) (node-induced subgraph match) and
[`edge_substruct_matches`](@ref) (edge-induced subgraph match, yielding bond mappings) are also
available. The bond mapping can be converted into an atom mapping by [`emaptonmap`](@ref).


## InChI and InChIKey

- [`inchi(mol)`](@ref inchi) and [`inchikey(mol)`](@ref inchikey) generate InChI and
  InChIKey by using the InChI library (libinchi) published by the InChI Trust and IUPAC.
- InChI and InChIKey provide structure-based identifiers for chemical databases and enable
  efficient exact structure match. See the
  [InChI Trust technical FAQ](https://www.inchi-trust.org/technical-faq/) for usage and
  limitations of InChI and InChIKey.

```@example substructure
inchi(thiazole)
```

```@example substructure
inchikey(thiazole)
```


## SMARTS query

SMARTS queries parsed by [`smartstomol`](@ref) can also be used as an argument of the
substructure match functions.

Nitrogens with one or two hydrogens:

```@example substructure
highlight_matches(thiostrepton, smartstomol(raw"[#7;H1,H2]"))
```

Amides next to thiazoles (recursive SMARTS):

```@example substructure
amide_query = smartstomol(raw"[$(Cc1ncsc1)](=O)N")
highlight_matches(thiostrepton, amide_query)
```

```@example substructure
has_substruct_match(thiostrepton, amide_query)
```


## Structural alerts

- MolecularGraph.jl has default query sets for functional group detection and structural
  alerts (substructures often avoided in medicinal chemistry; these are often called "PAINS",
  but PAINS is just one of them).
- Note that structural alerts are not a universal structural filter. Some of them are based
  solely on chemists' intuition, others flag compounds that interfere with only some assay
  systems.
- [`query_containment_diagram(; sources=[])`](@ref query_containment_diagram) generates the
  query containment diagram (described later) from the default query sets. It returns a tuple
  of a directed graph and the vertex and edge property dictionaries.
- The structural alert datasets were collected from ChEMBL. Which datasets to use can be
  specified by the `sources` keyword argument with the following keywords that indicate their
  origin: `BMS`, `Dundee`, `Glaxo`, `MLSMR` and `PAINS`. If `sources` is not specified or is an
  empty vector, all default datasets (including the functional group sets described below)
  are used.

```@example substructure
sa_diagram = query_containment_diagram()
nv(sa_diagram[1])
```

```@example substructure
curcumin = smilestomol(raw"O=C(\C=C\c1ccc(O)c(OC)c1)CC(=O)\C=C\c2cc(OC)c(O)cc2")
html_fixed_size(curcumin, 500, 300, atomindex=true)
```

[`find_queries(mol, query_diagram; sources=[])`](@ref find_queries) returns the queries
matched by the molecule, as a subgraph of the query containment diagram (a tuple of the graph,
vertex properties and edge properties). Each vertex property has the query key, the query
information and the matched atom indices.

```@example substructure
curcumin_sa = find_queries(
    curcumin, sa_diagram; sources=["BMS", "Dundee", "Glaxo", "MLSMR", "PAINS"])
for rcd in values(curcumin_sa[2])
    info = rcd["info"][1]
    println(rcd["key"], " | ", info["name"], " | ", info["source"], " | ",
        info["query"], " | ", rcd["matched"])
end
```


## Functional group analysis

Manually curated default functional group datasets are available for efficient calculation
and systematic annotation of functional groups. These datasets can be specified in the
`sources` keyword of [`query_containment_diagram`](@ref) with the following keywords:
`default`, `scaffold` and `biomolecule`.

```@example substructure
fg_diagram = query_containment_diagram(sources=["default", "scaffold", "biomolecule"])
nv(fg_diagram[1])
```

```@example substructure
cefditoren = sdftomol(joinpath(datadir, "Cefditoren Pivoxil.mol"))
remove_hydrogens!(cefditoren)
html_fixed_size(cefditoren, 500, 300, atomindex=true)
```

Thanks to the containment diagram, only the most specific functional groups (vertices without
incoming edges, i.e. queries not contained in other matched queries) can be listed.

```@example substructure
cefditoren_fg = find_queries(cefditoren, fg_diagram)
for n in vertices(cefditoren_fg[1])
    indegree(cefditoren_fg[1], n) == 0 || continue
    rcd = cefditoren_fg[2][n]
    println(rcd["key"], " | ", rcd["info"][1]["query"], " | ", rcd["matched"])
end
```


## Query containment

A carbonyl group query (`[#6]=O`) always matches molecules that match a carboxyl query
(`O=CO`). Similarly, a carboxyl query always matches molecules that match an ester query
(`[$(C[#6])](=O)[$(O(C)[#6])]`). Thus, if a particular query clearly returns a subset of the
search results of another query, the calculation of the former can be skipped when the
latter does not match. MolecularGraph.jl implements a mechanism to resolve such query
containment, which is used in [`find_queries`](@ref) (`filtering=true` by default).

Substructure match between two query molecules tells whether the second query is contained in
the first one.

```@example substructure
has_substruct_match(
    smartstomol(raw"[$(C[#6])](=O)[$(O(C)[#6])]"),
    smartstomol("O=CO")
)
```

At the atom (or bond) level, `issubset(q, r)` returns whether all the results of the atom
query `q` are included in the results of the atom query `r`.

```@example substructure
nh2 = smartstomol("[NH2]")[1]
nitrogen = smartstomol("[#7]")[1]
issubset(nh2, nitrogen), issubset(nitrogen, nh2)
```

```@example substructure
issubset(smartstomol(raw"[$(C=O)]")[1], smartstomol("[#6]")[1])
```

Atom query containment is evaluated by converting the logical query trees into a satisfiability
problem, so its cost may grow exponentially with the number of query literals of an atom.
