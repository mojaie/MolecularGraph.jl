# Maximum common substructure

MolecularGraph.jl implements the following essential maximum common substructure (MCS)
methods for cheminformatics applications.

- Maximum common induced substructure (MCIS)
- Maximum common edge-induced substructure (MCES)
- Connected or disconnected MCS
- Working with larger molecules
- Topological constraint (known as TD-MCS)

```@example mcs
using MolecularGraph, Graphs
```

The examples use molfiles bundled in the documentation sources, which were downloaded from
[PubChem](https://pubchem.ncbi.nlm.nih.gov/) (Levodopa: CID 6047, 3-Aminocoumarin: CID 74217,
Cefditoren Pivoxil: CID 6437877, Ceftazidime: CID 5481173).

```@example mcs
datadir = joinpath(pkgdir(MolecularGraph), "docs", "src", "assets", "data")

function loadmol(name)
    mol = sdftomol(joinpath(datadir, "$(name).mol"))
    remove_all_hydrogens!(mol)
    return mol
end

ldopa = loadmol("Levodopa")
aminocoumarin = loadmol("3-Aminocoumarin")
nothing # hide
```

To compare the results, the following helper draws two molecules side by side with the
matched atoms and bonds highlighted.

```@example mcs
function draw_pair(mol1, mol2, nodes1, edges1, nodes2, edges2; size=250, atomindex=false)
    svgs = [
        drawsvg(mol1; atomhighlight=nodes1, bondhighlight=edges1, atomindex=atomindex),
        drawsvg(mol2; atomhighlight=nodes2, bondhighlight=edges2, atomindex=atomindex)
    ]
    cells = ["<div style=\"width:$(size)px;height:$(size)px\">$(s)</div>" for s in svgs]
    return HTML("<div style=\"display:flex\">$(join(cells))</div>")
end

# highlight from a vertex mapping (MCIS)
function draw_nmap(mol1, mol2, nmap; kwargs...)
    nodes1, nodes2 = collect(keys(nmap)), collect(values(nmap))
    edges1 = induced_subgraph_edges(mol1.graph, nodes1)
    edges2 = induced_subgraph_edges(mol2.graph, nodes2)
    return draw_pair(mol1, mol2, nodes1, edges1, nodes2, edges2; kwargs...)
end

# highlight from an edge mapping (MCES)
function draw_emap(mol1, mol2, emap; kwargs...)
    edges1, edges2 = collect(keys(emap)), collect(values(emap))
    nodes1 = unique(vcat(src.(edges1), dst.(edges1)))
    nodes2 = unique(vcat(src.(edges2), dst.(edges2)))
    return draw_pair(mol1, mol2, nodes1, edges1, nodes2, edges2; kwargs...)
end
nothing # hide
```


## Maximum common induced substructure (MCIS)

- [`disconnected_mcis(mol1, mol2)`](@ref disconnected_mcis) calculates the maximum common
  induced (node-induced) substructure (MCIS).
- "Induced" (or specifically "node-induced") means that the subgraph is defined by a subset
  of the vertices of the original graph, and the presence or absence of edges between all
  vertices of the subgraph must match that of the original graph. That also means there is at
  least one vertex mapping (injection) from the subgraph to the original graph.
- MCS functions in this library return a result as a `Tuple{Dict,Symbol}`.
  - The first element is a mapping of vertex indices (or edges for MCES) between the two
    molecules.
  - The second element is the calculation status. `:done` means that the calculation
    successfully finished.
- Unlike substructure match, an MCS result has only one of the largest mappings. This is
  because of the high computational cost of MCS calculation, as described below.
- Like substructure match, MCS calculation considers the match of the atom symbol and the
  number of pi electrons on the atom (not the bond order). This is why the nitrogen atom in
  the following example does not match (N adjacent to a multiple or aromatic bond is treated
  as sp2).

```@example mcs
mcis_result, mcis_status = disconnected_mcis(ldopa, aminocoumarin)
```

```@example mcs
draw_nmap(ldopa, aminocoumarin, mcis_result, atomindex=true)
```


## Maximum common edge-induced substructure (MCES)

- [`disconnected_mces(mol1, mol2)`](@ref disconnected_mces) calculates the maximum common
  edge-induced substructure (MCES).
- Similar to the node-induced subgraph, "edge-induced" means that the subgraph is defined by a
  subset of the edges of the original graph. The vertex set of the edge-induced subgraph is a
  subset of the vertices of the original graph. Also, there is at least one edge mapping
  (injection) from the subgraph to the original graph.
- In many cases the MCIS size is sufficient as an indicator of the degree of structural match.
  On the other hand, many chemists may find MCES more intuitive than MCIS, so MCES may be
  worth using for visualization at a slightly higher computational cost.

```@example mcs
mces_result, mces_status = disconnected_mces(ldopa, aminocoumarin)
mces_status, length(mces_result)
```

```@example mcs
draw_emap(ldopa, aminocoumarin, mces_result, atomindex=true)
```


## Connected or disconnected MCS

- [`connected_mcis(mol1, mol2)`](@ref connected_mcis) and
  [`connected_mces(mol1, mol2)`](@ref connected_mces) also calculate MCIS and MCES, but with a
  connection constraint (the obtained MCS must be a connected component). RDKit FMCS also
  adopts this approach.
- Connected MCS is much less computationally expensive than disconnected MCS.
- As disconnected MCS does not care about the spatial relationship (distance) among matched
  fragments, connected MCS can be a more appropriate option for some kinds of applications.

```@example mcs
conn_mcis_result, conn_mcis_status = connected_mcis(ldopa, aminocoumarin)
draw_nmap(ldopa, aminocoumarin, conn_mcis_result, atomindex=true)
```

```@example mcs
conn_mces_result, conn_mces_status = connected_mces(ldopa, aminocoumarin)
draw_emap(ldopa, aminocoumarin, conn_mces_result, atomindex=true)
```


## Working with larger molecules

Then, let's try it with slightly larger molecules.

```@example mcs
cefditoren = loadmol("Cefditoren Pivoxil")
ceftazidime = loadmol("Ceftazidime")
nv(cefditoren), nv(ceftazidime)
```

As disconnected MCS may take a very long time, the calculation is interrupted when the time
specified by the keyword argument `timeout` (in seconds, default: 10) has been reached. Here
`timeout=5` is used.

```@example mcs
larger_mces, larger_mces_status = disconnected_mces(cefditoren, ceftazidime, timeout=5)
larger_mces_status, length(larger_mces)
```

If it timed out, the status will be `:timedout` and the result will be the mapping of the
suboptimal common substructure found so far.

```@example mcs
draw_emap(cefditoren, ceftazidime, larger_mces)
```

- Profiling the calculation shows that the most costly task is the maximum clique detection
  in the modular product graph, which is used to solve the maximum common subgraph problem.
  This problem is known to be NP-hard, which means that the computation time may increase
  exponentially with the number of nodes in the graph.
- As described above, connected MCS is far faster than disconnected MCS (but there is still a
  possibility of hang-up).

```@example mcs
larger_conn_mces, larger_conn_mces_status = connected_mces(cefditoren, ceftazidime, timeout=5)
larger_conn_mces_status, length(larger_conn_mces)
```

```@example mcs
draw_emap(cefditoren, ceftazidime, larger_conn_mces)
```

Alternatively, the `targetsize` keyword argument terminates the calculation with the
`:targetreached` status when the MCS reaches the target size. This feature is quite useful in
library screening.

```@example mcs
targeted_mces, targeted_mces_status = disconnected_mces(cefditoren, ceftazidime, targetsize=15)
targeted_mces_status, length(targeted_mces)
```

```@example mcs
draw_emap(cefditoren, ceftazidime, targeted_mces)
```


## Topological constraint (known as TD-MCS)

- Disconnected MCS methods can detect as many matched fragments as possible, but they do not
  reflect the spatial relationship of each fragment.
- On the other hand, connected MCS is too strict to allow distant matches.
- The graph distance can be used as a constraint on the spatial arrangement of matched
  substructure fragments (topological constraint). This feature seems to be preferable in
  pharmacophore matching and structural similarity screening.
- The topological constraint also effectively decreases the calculation cost (but note that,
  mathematically, this approach still has a possibility of hang-up, so the default timeout is
  set to 10 seconds).
- [`tdmcis`](@ref) and [`tdmces`](@ref) calculate MCIS and MCES with the topological
  constraint. The distance mismatch tolerance parameter ``θ`` is available as the keyword
  argument `tolerance`, and the distance cutoff as `diameter`.

```@example mcs
tdmces_result, tdmces_status = tdmces(cefditoren, ceftazidime, tolerance=1)
tdmces_status, length(tdmces_result)
```

```@example mcs
draw_emap(cefditoren, ceftazidime, tdmces_result)
```
