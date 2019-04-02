var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#MolecularGraph.jl-1",
    "page": "Home",
    "title": "MolecularGraph.jl",
    "category": "section",
    "text": "README.md on GitHubMolecular graph modeling and chemoinformatics toolkit\nFully implemented in Julia\nMIT License"
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": " (v1.0) pkg> add MolecularGraph"
},

{
    "location": "#Usage-1",
    "page": "Home",
    "title": "Usage",
    "category": "section",
    "text": "Jupyter notebook tutorials"
},

{
    "location": "moleculargraph/io/#",
    "page": "I/O",
    "title": "I/O",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/io/#Molecule-I/O-1",
    "page": "I/O",
    "title": "Molecule I/O",
    "category": "section",
    "text": ""
},

{
    "location": "moleculargraph/io/#Base.parse-Tuple{Type{MapMol{SDFileAtom,SDFileBond}},Any}",
    "page": "I/O",
    "title": "Base.parse",
    "category": "method",
    "text": "parse(::Type{SDFile}, lines)\n\nParse lines of a SDFile mol block data into a molecule object.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdfilereader-Tuple{AbstractString}",
    "page": "I/O",
    "title": "MolecularGraph.sdfilereader",
    "category": "method",
    "text": "sdfilereader(path::AbstractString)\n\nRead a SDFile and return a lazy iterator which yields molecule objects.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdfilereader-Tuple{IO}",
    "page": "I/O",
    "title": "MolecularGraph.sdfilereader",
    "category": "method",
    "text": "sdfilereader(file::IO)\n\nRead SDFile data from input stream and return a lazy iterator which yields molecule objects.\n\nsdfilereader does not stop and raise errors when an erroneous or incompatible SDFile block is read but produces an error message and yields an empty molecule. If this behavior is not desirable, you can use the customized supplier function instead of default supplier nohaltsupplier\n\nfunction customsupplier()\n    mol = try\n        parse(SDFile, block)\n    catch e\n        throw(ErrorException(\"incompatible molecule found, aborting...\"))\n    end\n    return defaultpostprocess(mol)\nend\n\nfunction sdfilereader(file::IO)\n    return SDFileReader(eachline(file), customsupplier)\nend\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdftomol-Tuple{AbstractString}",
    "page": "I/O",
    "title": "MolecularGraph.sdftomol",
    "category": "method",
    "text": "sdftomol(path::AbstractString)\n\nRead a SDFile and parse it into a molecule object. Single mol block files without optional information are often provided as a .mol file.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdftomol-Tuple{IO}",
    "page": "I/O",
    "title": "MolecularGraph.sdftomol",
    "category": "method",
    "text": "sdftomol(file::IO)\n\nRead a SDFile mol block from the input stream and parse it into a molecule object.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#SDFile-reader-1",
    "page": "I/O",
    "title": "SDFile reader",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"sdfilereader.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdfilewriter-Tuple{AbstractString,Any}",
    "page": "I/O",
    "title": "MolecularGraph.sdfilewriter",
    "category": "method",
    "text": "sdfilewriter(path::AbstractString)\n\nRead a SDFile and return a lazy iterator which yields molecule objects.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdfilewriter-Tuple{IO,Any}",
    "page": "I/O",
    "title": "MolecularGraph.sdfilewriter",
    "category": "method",
    "text": "sdfilewriter(mols::MolGraph, file::IO)\n\nWrite molecule data to the output stream as a SDFile format file.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#SDFile-writer-1",
    "page": "I/O",
    "title": "SDFile writer",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"sdfilewriter.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/draw/#",
    "page": "Structure drawing",
    "title": "Structure drawing",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/draw/#MolecularGraph.DRAW_SETTING",
    "page": "Structure drawing",
    "title": "MolecularGraph.DRAW_SETTING",
    "category": "constant",
    "text": "DRAW_SETTING\n\nDefault setting parameters of the molecule drawing canvas.\n\nRequired fields\n\n:display_terminal_carbon(Bool) whether to display terminal C atom or not\n:double_bond_notation(Symbol)\n:alongside: all double bonds are represented as a carbon skeleton and               a segment alongside it.\n:dual: all double bonds are represented as two equal length parallel          segments.\n:chain: :dual for chain bonds and :alongside for ring bonds\n:terminal: :dual for terminal bonds (adjacent to degree=1 node) and              :alongside for others (default)\n:atomcolor(Dict{Symbol,Color}) atom symbol and bond colors for organic atoms\n:defaul_atom_color(Dict{Symbol,Color}) colors for other atoms\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.atomcolor-Tuple{VectorMol}",
    "page": "Structure drawing",
    "title": "MolecularGraph.atomcolor",
    "category": "method",
    "text": "atomcolor(mol::VectorMol; setting=DRAW_SETTING)\n\nReturn atom colors for molecule 2D drawing\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.chargesign-Tuple{Int64}",
    "page": "Structure drawing",
    "title": "MolecularGraph.chargesign",
    "category": "method",
    "text": "chargesign(charge::Int) -> String\n\nGet a charge sign.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.draw2d!-Tuple{Canvas,VectorMol}",
    "page": "Structure drawing",
    "title": "MolecularGraph.draw2d!",
    "category": "method",
    "text": "draw2d!(canvas::Canvas, mol::VectorMol;\n        setting=copy(DRAW_SETTING), recalculate=false)\n\nDraw molecular image to the canvas.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.drawsvg-Tuple{VectorMol,Int64,Int64}",
    "page": "Structure drawing",
    "title": "MolecularGraph.drawsvg",
    "category": "method",
    "text": "drawsvg(mol::VectorMol, width::Int, height::Int)\n\nGenerate molecular structure image as a SVG format string.\n\nwidth and height specifies the size of the image (width and height attribute of svg tag).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.initcanvas!-Tuple{Canvas,VectorMol}",
    "page": "Structure drawing",
    "title": "MolecularGraph.initcanvas!",
    "category": "method",
    "text": "initcanvas!(canvas::Canvas, mol::VectorMol)\n\nMove and adjust the size of the molecule for drawing.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#Molecular-structure-drawing-1",
    "page": "Structure drawing",
    "title": "Molecular structure drawing",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"./draw/base.jl\", \"./draw/draw2d.jl\", \"./draw/svg.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/properties/#",
    "page": "Basic chemical properties",
    "title": "Basic chemical properties",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/properties/#MolecularGraph.hydrogen_acceptor_count-Tuple{VectorMol}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.hydrogen_acceptor_count",
    "category": "method",
    "text": "hydrogen_acceptor_count(mol::VectorMol) -> Int\n\nReturn the number of hydrogen bond acceptors (N, O and F).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.hydrogen_donor_count-Tuple{VectorMol}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.hydrogen_donor_count",
    "category": "method",
    "text": "hydrogen_donor_count(mol::VectorMol) -> Int\n\nReturn the number of hydrogen bond donors (O and N attached to hydrogens).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.molweight-Tuple{VectorMol}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.molweight",
    "category": "method",
    "text": "molweight(mol::VectorMol; digits=2) -> Float64\n\nReturn standard molecular weight.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.rotatable_count-Tuple{VectorMol}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.rotatable_count",
    "category": "method",
    "text": "rotatable_count(mol::VectorMol) -> Int\n\nReturn the number of rotatable bonds.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.wildman_crippen_logp-Tuple{VectorMol}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.wildman_crippen_logp",
    "category": "method",
    "text": "wildman_crippen_logp(mol::VectorMol) -> Float64\n\nReturn predicted logP value calculated by using Wildman and Crippen method.\n\nReference\n\nWildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical Parameters by Atomic Contributions. Journal of Chemical Information and Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#Molecular-properties-1",
    "page": "Basic chemical properties",
    "title": "Molecular properties",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"properties.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/descriptor/#",
    "page": "Molecular descriptor",
    "title": "Molecular descriptor",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/descriptor/#Molecular-descriptor-1",
    "page": "Molecular descriptor",
    "title": "Molecular descriptor",
    "category": "section",
    "text": ""
},

{
    "location": "moleculargraph/descriptor/#Property-vector-1",
    "page": "Molecular descriptor",
    "title": "Property vector",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\n  \"./annotation/base.jl\", \"./annotation/aromatic.jl\",\n  \"./annotation/elemental.jl\", \"./annotation/rotatable.jl\",\n  \"./annotation/topology.jl\", \"./annotation/wclogp.jl\"\n]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/preprocess/#",
    "page": "Preprocessing",
    "title": "Preprocessing",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.all_hydrogens-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.all_hydrogens",
    "category": "method",
    "text": "all_hydrogens(mol::VectorMol) -> Set{Int}\n\nReturn a set of hydrogen nodes.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.canonicalize!-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.canonicalize!",
    "category": "method",
    "text": "canonicalize!(mol::VectorMol)\n\nCanonicalize molecule notation and apply the changes to the molecular property vector.\n\nNeutralize oxo acid, 1-3° ammonium and polarized carbonyls except in the case that polarization is required for aromaticity.\nCanonicalize anions next to triple bonds (ex. [C-][N+]#N -> C=[N+]=[N-])\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.depolarize!-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.depolarize!",
    "category": "method",
    "text": "depolarize!(mol::VectorMol)\n\nDepolarize oxo groups except in the case that polarization is required for aromaticity.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.largest_component_nodes-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.largest_component_nodes",
    "category": "method",
    "text": "largest_component_nodes(mol::VectorMol) -> Set{Int}\n\nReturn a set of nodes in the largest connected component.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.largestcomponent-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.largestcomponent",
    "category": "method",
    "text": "largestcomponent(mol::VectorMol) -> VectorMol\n\nReturn largest connected component of the molecular graph.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.make_hydrogens_explicit-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.make_hydrogens_explicit",
    "category": "method",
    "text": "make_hydrogens_explicit(mol::VectorMol) -> VectorMol\n\nReturn molecule whose hydrogens are fully attached. If option all is set to false, only trivial hydrogens are removed (see trivialhydrogens).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.make_hydrogens_implicit-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.make_hydrogens_implicit",
    "category": "method",
    "text": "make_hydrogens_implicit(mol::VectorMol) -> VectorMol\n\nReturn molecule whose hydrogen nodes are removed. If option all is set to false, only trivial hydrogens are removed (see trivialhydrogens).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.neutralize_acids!-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.neutralize_acids!",
    "category": "method",
    "text": "neutralize_acids!(mol::VectorMol)\n\nNeutralize oxo(thio) acids.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.neutralize_oniums!-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.neutralize_oniums!",
    "category": "method",
    "text": "neutralize_oniums!(mol::VectorMol)\n\nNeutralize 1-3° oniums. Permanently charged quart-oniums are not neutralized.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.triplebond_anion!-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.triplebond_anion!",
    "category": "method",
    "text": "triplebond_anion!(mol::VectorMol)\n\nCanonicalize anions next to triple bonds (ex. [C-][N+]#N -> C=[N+]=[N-]).\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.trivialhydrogens-Tuple{VectorMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.trivialhydrogens",
    "category": "method",
    "text": "trivialhydrogens(mol::VectorMol) -> Set{Int}\n\nReturn a set of trivial hydrogen nodes (light hydrogens which are uncharged, non-radical, non-stereospecific and attached to organic heavy atoms)\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#Preprocessing-1",
    "page": "Preprocessing",
    "title": "Preprocessing",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"preprocess.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/structure/#",
    "page": "Structure match",
    "title": "Structure match",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/structure/#MolecularGraph.is_identical-Tuple{GeneralMol,GeneralMol}",
    "page": "Structure match",
    "title": "MolecularGraph.is_identical",
    "category": "method",
    "text": "is_identical(mol1::GeneralMol, mol2::GeneralMol)\n\nReturn whether mol1 and mol2 are identical in chemical structure.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/structure/#MolecularGraph.is_querymatch-Tuple{Any,Any}",
    "page": "Structure match",
    "title": "MolecularGraph.is_querymatch",
    "category": "method",
    "text": "is_querymatch(mol, query; kwargs...)\n\nReturn whether mol matches with the query.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/structure/#MolecularGraph.is_substruct-Tuple{GeneralMol,GeneralMol}",
    "page": "Structure match",
    "title": "MolecularGraph.is_substruct",
    "category": "method",
    "text": "is_substruct(mol1::GeneralMol, mol2::GeneralMol)\n\nReturn whether mol1 is a substructure of mol2.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/structure/#MolecularGraph.is_superstruct-Tuple{Any,Any}",
    "page": "Structure match",
    "title": "MolecularGraph.is_superstruct",
    "category": "method",
    "text": "is_superstruct(mol1, mol2)\n\nReturn whether mol1 is a superstructure of mol2.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/structure/#Structure-match-1",
    "page": "Structure match",
    "title": "Structure match",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"substructure.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/model/#",
    "page": "Models",
    "title": "Models",
    "category": "page",
    "text": ""
},

{
    "location": "graph/model/#Graph-models-1",
    "page": "Models",
    "title": "Graph models",
    "category": "section",
    "text": ""
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.mapgraph-Tuple{Any,Any}",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.mapgraph",
    "category": "method",
    "text": "mapgraph(nodes, edges) -> MapGraph{Node,Edge}\n\nGenerate map graph that has given nodes and edges represented by the list of node indices in integer and the list of pairs of node indices, respectively.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.mapgraph-Tuple{Union{Graph, GraphView}}",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.mapgraph",
    "category": "method",
    "text": "mapgraph(graph::UndirectedGraph; clone=false) -> MapGraph\n\nConvert the given graph into a new MapGraph. The node type and edge type are inherited from the original graph. If the given graph is MapGraph, return a copy of the graph.\n\nIf you really need mutable nodes and edges, and want the graph with deepcopied elements, implement clone and setnodes as deepcopy methods for nodes and edges, respectively.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.mapgraph-Union{Tuple{E}, Tuple{N}, Tuple{Type{N},Type{E}}} where E<:MolecularGraph.MolecularGraphModel.UndirectedEdge where N<:MolecularGraph.MolecularGraphModel.AbstractNode",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.mapgraph",
    "category": "method",
    "text": "mapgraph(::Type{N}, ::Type{E}\n    ) where {N<:AbstractNode,E<:UndirectedEdge} -> MapGraph{N,E}()\n\nGenerate empty map graph that has nodes and edges with the given types.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.vectorgraph-Tuple{Int64,Any}",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.vectorgraph",
    "category": "method",
    "text": "vectorgraph(nodes, edges) -> VectorGraph{Node,Edge}\n\nGenerate vector graph that have given nodes and edges represented by the list of node indices in integer and the list of pairs of node indices, respectively.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.vectorgraph-Tuple{Union{Graph, GraphView}}",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.vectorgraph",
    "category": "method",
    "text": "vectorgraph(graph::UndirectedGraph; clone=false) -> VectorGraph\n\nConvert the given graph into a new VectorGraph. The node type and edge type are inherited from the original graph.\n\nNode indices are sorted in ascending order and are re-indexed. This behavior is intended for some cannonicalization operations (ex. chirality flag).\n\nIf you really need mutable nodes and edges, and want the graph with deepcopied elements, implement clone and setnodes as deepcopy methods for nodes and edges, respectively.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.vectorgraph-Union{Tuple{E}, Tuple{N}, Tuple{Type{N},Type{E}}} where E<:MolecularGraph.MolecularGraphModel.UndirectedEdge where N<:MolecularGraph.MolecularGraphModel.AbstractNode",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.vectorgraph",
    "category": "method",
    "text": "vectorgraph(::Type{N}, ::Type{E}\n    ) where {N<:AbstractNode,E<:UndirectedEdge} -> VectorGraph{N,E}()\n\nGenerate empty vector graph that have nodes and edges with the given types.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#Undirected-graph-1",
    "page": "Models",
    "title": "Undirected graph",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/ugraph.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.mapdigraph-Tuple{Any,Any}",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.mapdigraph",
    "category": "method",
    "text": "mapdigraph(nodes, edges) -> MapDiGraph{Node,Arrow}\n\nGenerate MapDiGraph that have given nodes and edges represented by the list of node indices in integer and the list of pairs of node indices, respectively.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.mapdigraph-Tuple{Union{DiGraph, DiGraphView}}",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.mapdigraph",
    "category": "method",
    "text": "mapdigraph(graph::DirectedGraph; clone=false) -> MapDiGraph\n\nConvert the given graph into a new MapDiGraph. The node type and edge type are inherited from the original graph. If the given graph is MapDiGraph, return a copy of the graph.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.mapdigraph-Union{Tuple{E}, Tuple{N}, Tuple{Type{N},Type{E}}} where E<:MolecularGraph.MolecularGraphModel.DirectedEdge where N<:MolecularGraph.MolecularGraphModel.AbstractNode",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.mapdigraph",
    "category": "method",
    "text": "mapdigraph(::Type{N}, ::Type{E}\n    ) where {N<:AbstractNode,E<:DirectedEdge} -> MapDiGraph{N,E}()\n\nGenerate empty MapDiGraph that have nodes and edges with the given types.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#Directed-graph-1",
    "page": "Models",
    "title": "Directed graph",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/dgraph.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.mapmultigraph-Tuple{Any,Any}",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.mapmultigraph",
    "category": "method",
    "text": "mapmultigraph(nodes, edges) -> MapMultiGraph{Node,Edge}\n\nGenerate MapMultiGraph that have given nodes and edges represented by the list of node indices in integer and the list of pairs of node indices, respectively.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#MolecularGraph.MolecularGraphModel.mapmultigraph-Union{Tuple{E}, Tuple{N}, Tuple{Type{N},Type{E}}} where E<:MolecularGraph.MolecularGraphModel.UndirectedEdge where N<:MolecularGraph.MolecularGraphModel.AbstractNode",
    "page": "Models",
    "title": "MolecularGraph.MolecularGraphModel.mapmultigraph",
    "category": "method",
    "text": "mapmultigraph(::Type{N}, ::Type{E}\n    ) where {N<:AbstractNode,E<:UndirectedEdge} -> MapMultiGraph{N,E}()\n\nGenerate empty MapMultiGraph that have nodes and edges with the given types.\n\n\n\n\n\n"
},

{
    "location": "graph/model/#Multi-graph-1",
    "page": "Models",
    "title": "Multi graph",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/multigraph.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/interface/#",
    "page": "Interface",
    "title": "Interface",
    "category": "page",
    "text": ""
},

{
    "location": "graph/interface/#Interface-1",
    "page": "Interface",
    "title": "Interface",
    "category": "section",
    "text": ""
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgecount",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgecount",
    "category": "function",
    "text": "edgecount(graph) -> Int\n\nReturn the number of graph edges.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgekeys",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgekeys",
    "category": "function",
    "text": "edgekeys(graph) -> Vector{Int}\n\nReturn graph edge keys. If the given graph is a vector graph, the keys are in ascending order, whereas the order of indices in map graph is not guaranteed.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgeset",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgeset",
    "category": "function",
    "text": "edgeset(graph) -> Set{Int}\n\nReturn the set of edge keys.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgesiter",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgesiter",
    "category": "function",
    "text": "edgesiter(graph)\n\nAn iterator that yields (i, e) where i is the edge index, and e is the edge object at the index i within the graph.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgetype",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgetype",
    "category": "function",
    "text": "edgetype(graph)\n\nReturn the edge type of the graph\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.getedge-Tuple{Any,Any,Any}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.getedge",
    "category": "method",
    "text": "getedge(graph, u, v) -> AbstractEdge\n\nRetrieve an edge object which connects the given nodes.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.getedge-Tuple{Any,Any}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.getedge",
    "category": "method",
    "text": "getedge(graph, index) -> AbstractEdge\n\nRetrieve the edge object at the given index within the graph.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.getnode",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.getnode",
    "category": "function",
    "text": "getnode(graph, index) -> AbstractNode\n\nRetrieve the node object at the given index within the graph.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.hasedge",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.hasedge",
    "category": "function",
    "text": "hasedge(graph, u, v) -> AbstractEdge\n\nReturn whether the given two nodes are connected by at least one edge.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.indegree",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.indegree",
    "category": "function",
    "text": "indegree(graph, n) -> Int\n\nReturn the number of inneighbors of the node \'n\'.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.inneighbors",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.inneighbors",
    "category": "function",
    "text": "inneighbors(graph, n) -> Dict{Int,Int}\n\nReturn the mapping of predecessor node keys and in edge keys connected to the given node.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.neighborcount",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.neighborcount",
    "category": "function",
    "text": "neighborcount(graph, n) -> Int\ndegree(graph, n) -> Int\n\nReturn the number of adjacent nodes of the node \'n\'.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.neighbors",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.neighbors",
    "category": "function",
    "text": "neighbors(graph, n) -> Dict{Int,Int}\n\nReturn the mapping of adjacent node keys and incident edge keys connected to the given node. If the graph is directed graph, both outneighbors and inneighbors are mapped.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodecount",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodecount",
    "category": "function",
    "text": "nodecount(graph) -> Int\n\nReturn the number of graph nodes.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodekeys",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodekeys",
    "category": "function",
    "text": "nodekeys(graph) -> Vector{Int}\n\nReturn graph node keys. If the given graph is a vector graph, the keys are in ascending order, whereas the order of indices in map graph is not guaranteed.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodeset",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodeset",
    "category": "function",
    "text": "nodeset(graph) -> Set{Int}\n\nReturn the set of node keys.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodesiter",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodesiter",
    "category": "function",
    "text": "nodesiter(graph)\n\nAn iterator that yields (i, n) where i is the node index, and n is the node object at the index i within the graph.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodetype",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodetype",
    "category": "function",
    "text": "nodetype(graph)\n\nReturn the node type of the graph\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.outdegree",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.outdegree",
    "category": "function",
    "text": "outdegree(graph, n) -> Int\n\nReturn the number of outneighbors of the node \'n\'.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.outneighbors",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.outneighbors",
    "category": "function",
    "text": "outneighbors(graph, n) -> Dict{Int,Int}\n\nReturn the mapping of successor node keys and out edge keys connected to the given node.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.unlinkedge!-Tuple{Any,Any,Any}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.unlinkedge!",
    "category": "method",
    "text": "unlinkedge!(graph, u, v)\n\nDelete the edge that connect nodes u and v.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.unlinkedge!-Tuple{Any,Any}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.unlinkedge!",
    "category": "method",
    "text": "unlinkedge!(graph, e)\n\nDelete the edge at the index of e.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.unlinknode!",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.unlinknode!",
    "category": "function",
    "text": "unlinknode!(graph, n)\n\nDelete the node at the index of n and its incident edges.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.updateedge!-NTuple{4,Any}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.updateedge!",
    "category": "method",
    "text": "updateedge!(graph, edge, u, v)\n\nRebind the edge that connects nodes u and v by the given edge object. If the nodes do not exist, throws KeyError.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.updateedge!-Tuple{Any,Any,Any}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.updateedge!",
    "category": "method",
    "text": "updateedge!(graph, edge, e)\n\nRebind the edge object stored at e of the graph by the given edge object. If the index does not exist, add new edge to the position e.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.updatenode!",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.updatenode!",
    "category": "function",
    "text": "updatenode!(graph, node, n)\n\nRebind the node object stored at n of the graph by the given node object. If the index does not exist, add new node to the position n.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#Interface-2",
    "page": "Interface",
    "title": "Interface",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/interface.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/generator/#",
    "page": "Generator",
    "title": "Generator",
    "category": "page",
    "text": ""
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.bipartitegraph-Tuple{Int64,Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.bipartitegraph",
    "category": "method",
    "text": "bipartitegraph(m::Int,n::Int) -> MapGraph{Node,Edge}\n\nGenerate bipartite graph K_mn.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.circularladder-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.circularladder",
    "category": "method",
    "text": "circularladder(n::Int) -> MapGraph{Node,Edge}\n\nGenerate circular ladder graph CL_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.completegraph-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.completegraph",
    "category": "method",
    "text": "completegraph(length::Int) -> MapGraph{Node,Edge}\n\nGenerate complete graph K_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.cyclegraph-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.cyclegraph",
    "category": "method",
    "text": "cyclegraph(length::Int) -> MapGraph{Node,Edge}\n\nGenerate cycle graph C_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.laddergraph-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.laddergraph",
    "category": "method",
    "text": "laddergraph(n::Int) -> MapGraph{Node,Edge}\n\nGenerate ladder graph L_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.moebiusladder-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.moebiusladder",
    "category": "method",
    "text": "moebiusladder(n::Int) -> MapGraph{Node,Edge}\n\nGenerate Möbius ladder graph ML_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.pathgraph-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.pathgraph",
    "category": "method",
    "text": "pathgraph(n::Int) -> MapGraph{Node,Edge}\n\nGenerate path graph P_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#Graph-generator-1",
    "page": "Generator",
    "title": "Graph generator",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/generator.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/shortestpath/#",
    "page": "Shortest path",
    "title": "Shortest path",
    "category": "page",
    "text": ""
},

{
    "location": "graph/shortestpath/#Shortest-path-1",
    "page": "Shortest path",
    "title": "Shortest path",
    "category": "section",
    "text": ""
},

{
    "location": "graph/shortestpath/#MolecularGraph.MolecularGraphModel.diameter-Tuple{Union{Graph, GraphView}}",
    "page": "Shortest path",
    "title": "MolecularGraph.MolecularGraphModel.diameter",
    "category": "method",
    "text": "diameter(graph::UndirectedGraph) => Int\n\nCompute the diameter of the graph (the largest eccentricity of any nodes).\n\n\n\n\n\n"
},

{
    "location": "graph/shortestpath/#MolecularGraph.MolecularGraphModel.distance-Tuple{Union{Graph, GraphView},Any}",
    "page": "Shortest path",
    "title": "MolecularGraph.MolecularGraphModel.distance",
    "category": "method",
    "text": "distance(graph::UndirectedGraph, root) -> Dict{Int,Union{Int,Nothing}}\n\nCompute the distance from root to any other nodes. If the nodes are not reachable each other, the value will be nothing.\n\n\n\n\n\n"
},

{
    "location": "graph/shortestpath/#MolecularGraph.MolecularGraphModel.distancematrix-Tuple{Union{Graph, GraphView}}",
    "page": "Shortest path",
    "title": "MolecularGraph.MolecularGraphModel.distancematrix",
    "category": "method",
    "text": "distancematrix(graph::UndirectedGraph) -> Matrix{Float64}\n\nCompute the distance among each other nodes.\n\nNote that the type of the generated matrix will be Float64. If the nodes are not reachable each other, the distance value will be Inf.\n\n\n\n\n\n"
},

{
    "location": "graph/shortestpath/#MolecularGraph.MolecularGraphModel.eccentricity-Tuple{Union{Graph, GraphView},Any}",
    "page": "Shortest path",
    "title": "MolecularGraph.MolecularGraphModel.eccentricity",
    "category": "method",
    "text": "eccentricity(graph::UndirectedGraph, v) => Int\n\nCompute the eccentricity of the graph (the largest distance between v and any other nodes).\n\n\n\n\n\n"
},

{
    "location": "graph/shortestpath/#MolecularGraph.MolecularGraphModel.longestshortestpath-Tuple{Union{Graph, GraphView}}",
    "page": "Shortest path",
    "title": "MolecularGraph.MolecularGraphModel.longestshortestpath",
    "category": "method",
    "text": "longestshortestpath(graph::UndirectedGraph) -> Vector{Int}\n\nCompute the longest shortest path in the graph (a path between two arbitrary peripheral nodes) as a vector of nodes that starts with one of the peripheral node and ends with the other side.\n\n\n\n\n\n"
},

{
    "location": "graph/shortestpath/#MolecularGraph.MolecularGraphModel.shortestpath-Tuple{Union{Graph, GraphView},Any,Any}",
    "page": "Shortest path",
    "title": "MolecularGraph.MolecularGraphModel.shortestpath",
    "category": "method",
    "text": "shortestpath(graph::UndirectedGraph, u, v) -> Vector{Int}\n\nCompute the shortest path between u and v as a vector of the nodes that starts with u and ends withv. Return nothing if u == v or not reachable.\n\n\n\n\n\n"
},

{
    "location": "graph/shortestpath/#Shortest-path-2",
    "page": "Shortest path",
    "title": "Shortest path",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/shortestpath.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/linegraph/#",
    "page": "Line graph",
    "title": "Line graph",
    "category": "page",
    "text": ""
},

{
    "location": "graph/linegraph/#Line-graph-1",
    "page": "Line graph",
    "title": "Line graph",
    "category": "section",
    "text": ""
},

{
    "location": "graph/linegraph/#MolecularGraph.MolecularGraphModel.linegraph-Tuple{Union{Graph, GraphView}}",
    "page": "Line graph",
    "title": "MolecularGraph.MolecularGraphModel.linegraph",
    "category": "method",
    "text": "linegraph(G::UndirectedGraph) -> MapGraph{LineGraphNode,LineGraphEdge}\n\nGenerate line graph.\n\n\n\n\n\n"
},

{
    "location": "graph/linegraph/#Line-graph-2",
    "page": "Line graph",
    "title": "Line graph",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/linegraph.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/connectivity/#",
    "page": "Connectivity",
    "title": "Connectivity",
    "category": "page",
    "text": ""
},

{
    "location": "graph/connectivity/#Connectivity-1",
    "page": "Connectivity",
    "title": "Connectivity",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/connectivity.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/planarity/#",
    "page": "Planarity",
    "title": "Planarity",
    "category": "page",
    "text": ""
},

{
    "location": "graph/planarity/#MolecularGraph.MolecularGraphModel.is_outerplanar-Tuple{Union{Graph, GraphView}}",
    "page": "Planarity",
    "title": "MolecularGraph.MolecularGraphModel.is_outerplanar",
    "category": "method",
    "text": "is_outerplanar(graph::UndirectedGraph) -> Bool\n\nReturn whether the graph is outerplanar. The outerplanarity test is based on a planarity test (see is_planar).\n\n\n\n\n\n"
},

{
    "location": "graph/planarity/#MolecularGraph.MolecularGraphModel.is_planar-Tuple{Union{Graph, GraphView}}",
    "page": "Planarity",
    "title": "MolecularGraph.MolecularGraphModel.is_planar",
    "category": "method",
    "text": "is_planar(graph::UndirectedGraph) -> Bool\n\nReturn whether the graph is planar.\n\nReference\n\nde Fraysseix, H., & Ossona de Mendez, P. (2012). Trémaux trees and planarity. European Journal of Combinatorics, 33(3), 279–293. https://doi.org/10.1016/j.ejc.2011.09.012\n\n\n\n\n\n"
},

{
    "location": "graph/planarity/#Planarity-1",
    "page": "Planarity",
    "title": "Planarity",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/planarity.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/clique/#",
    "page": "Clique",
    "title": "Clique",
    "category": "page",
    "text": ""
},

{
    "location": "graph/clique/#Clique-1",
    "page": "Clique",
    "title": "Clique",
    "category": "section",
    "text": ""
},

{
    "location": "graph/clique/#MolecularGraph.MolecularGraphModel.maxclique-Tuple{Union{Graph, GraphView}}",
    "page": "Clique",
    "title": "MolecularGraph.MolecularGraphModel.maxclique",
    "category": "method",
    "text": "maxclique(graph::UndirectedGraph; kwargs...) -> Set{Int}\n\nCompute maximum clique of the graph. For details, see maximalcliques.\n\n\n\n\n\n"
},

{
    "location": "graph/clique/#MolecularGraph.MolecularGraphModel.maximalcliques-Tuple{Union{Graph, GraphView}}",
    "page": "Clique",
    "title": "MolecularGraph.MolecularGraphModel.maximalcliques",
    "category": "method",
    "text": "maximalcliques(graph::UndirectedGraph; kwargs...)\n\nReturn Channel which generates maximal cliques of the graph. Each cliques are represented as a Set of member nodes.\n\nReference\n\nTomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time complexity for generating all maximal cliques and computational experiments. Theoretical Computer Science, 363(1), 28–42. https://doi.org/10.1016/J.TCS.2006.06.015\nCazals, F., & Karande, C. (2005). An algorithm for reporting maximal c-cliques. Theoretical Computer Science, 349(3), 484–490. https://doi.org/10.1016/j.tcs.2005.09.038\n\n\n\n\n\n"
},

{
    "location": "graph/clique/#Clique-2",
    "page": "Clique",
    "title": "Clique",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/clique.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "python/#",
    "page": "Python interface",
    "title": "Python interface",
    "category": "page",
    "text": ""
},

{
    "location": "python/#Python-interface-1",
    "page": "Python interface",
    "title": "Python interface",
    "category": "section",
    "text": ""
},

{
    "location": "design/#",
    "page": "Design of molecular graph models",
    "title": "Design of molecular graph models",
    "category": "page",
    "text": ""
},

{
    "location": "design/#Design-of-molecular-graph-models-1",
    "page": "Design of molecular graph models",
    "title": "Design of molecular graph models",
    "category": "section",
    "text": ""
},

{
    "location": "design/#Graph-types-1",
    "page": "Design of molecular graph models",
    "title": "Graph types",
    "category": "section",
    "text": "(AbstractGraph)\n(Graph)\nMapGraph\nVectorGraph\n(DiGraph)\nMapDiGraph\nVectorDiGraph\n(MultiGraph)\nMapMultiGraph\nVectorMultiGraph\n(MultiDiGraph)\nMapMultiDiGraph\nVectorMultiDiGraph\n(GraphView)\nSubgraphView\nComplementGraphView\nMolGraph\n(DiGraphView)\nDiSubgraphView\nComplementDiGraphView\nReverseGraphView\n(UndirectedGraph)\n(Graph)\n(MultiGraph)\n(GraphView)\n(DirectedGraph)\n(DiGraph)\n(MultiDiGraph)\n(DiGraphView)"
},

{
    "location": "design/#Molecule-1",
    "page": "Design of molecular graph models",
    "title": "Molecule",
    "category": "section",
    "text": "MolGraph\nMapMolGraph\nGeneralMapMol{A<:Atom, B<:Bond}\nSDFile (Alias of GeneralMapMol{SDFileAtom, SDFileBond})\nSMILES (Alias of GeneralMapMol{SmilesAtom, SmilesBond})\nQueryMolGraph\nConnectedQueryMol{A<:QueryAtom, B<:QueryBond}\nConnectedSMARTS (Alias of ConnectedQueryMol{SmartsAtom, SmartsBond})\nDisconnectedQueryMol{A<:QueryAtom, B<:QueryBond}\nSMARTS (Alias of DisconnectedQueryMol{SmartsAtom, SmartsBond})\nVectorMolGraph\nGeneralVectorMol{A<:Atom, B<:Bond}\nMapMolView\nQueryMolView\nVectorMolView"
},

{
    "location": "design/#Methods-1",
    "page": "Design of molecular graph models",
    "title": "Methods",
    "category": "section",
    "text": "getatom\ngetbond\nneighbors\nneighborcount (or degree)\natomcount\nbondcount\nupdateatom!\nupdatebond!\nunlinkatom!\nunlinkbond!"
},

{
    "location": "design/#VectorMol-1",
    "page": "Design of molecular graph models",
    "title": "VectorMol",
    "category": "section",
    "text": "VectorMol is vector(array)-based molecular model which is specialized for element-wise fast computation of molecular properties.VectorMol can iterate over atom properties faster than MapMoland can store calculated molecular properties and annotation arrays which are suitable for vector computation. On the other hand, VectorMol does not have abilities to modify its graph structure (adding or removing elements). VectorMol can be converted to MapMol but the calculated properties and annotations will be lost."
},

{
    "location": "design/#MapMol-1",
    "page": "Design of molecular graph models",
    "title": "MapMol",
    "category": "section",
    "text": "MapMol is used as a molecular model builder for general purpose.This type inherits AbstractMapMol, a molecular graph model which have map(dict)-based structure. The map-based molecular graph can insert and delete elements (atoms and bonds). This can be easily converted to VectorMol object by using vectormol method"
},

{
    "location": "design/#QueryMol-1",
    "page": "Design of molecular graph models",
    "title": "QueryMol",
    "category": "section",
    "text": "QueryMol consists of QueryAtoms and QueryBonds that represents molecular query (ex. atom symbol is \'O\' and charge is -1, bond order is 1 and not in rings, ...). This type of objects typically built from SMARTS query."
},

{
    "location": "design/#Atom-1",
    "page": "Design of molecular graph models",
    "title": "Atom",
    "category": "section",
    "text": "AbstractNode\nAbstractAtom\nAtom\nSDFileAtom\nSmilesAtom\nQueryAtom\nSmartsAtom"
},

{
    "location": "design/#Bond-1",
    "page": "Design of molecular graph models",
    "title": "Bond",
    "category": "section",
    "text": "AbstractEdge\nAbstractBond\nBond\nSDFileBond\nSmilesBond\nQueryBond\nSmartsBond"
},

]}
