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
    "location": "moleculargraph/io/#MolecularGraph.sdfilereader-Tuple{IO}",
    "page": "I/O",
    "title": "MolecularGraph.sdfilereader",
    "category": "method",
    "text": "sdfilereader(file::IO)\nsdfilereader(path::AbstractString)\n\nRead SDFile data from input stream (or a file path as a string) and return a lazy iterator that yields molecule objects.\n\nsdfilereader does not stop and raise errors when an erroneous or incompatible SDFile block is read but produces an error message and yields an empty molecule. If this behavior is not desirable, you can use the customized supplier function instead of default supplier nohaltsupplier\n\nfunction customsupplier()\n    mol = try\n        parse(SDFile, block)\n    catch e\n        throw(ErrorException(\"incompatible molecule found, aborting...\"))\n    end\nend\n\nfunction sdfilereader(file::IO)\n    return SDFileReader(eachline(file), customsupplier)\nend\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdftomol-Tuple{Any}",
    "page": "I/O",
    "title": "MolecularGraph.sdftomol",
    "category": "method",
    "text": "sdftomol(lines) -> GraphMol{SDFileAtom,SDFileBond}\nsdftomol(file::IO) -> GraphMol{SDFileAtom,SDFileBond}\nsdftomol(path::AbstractString) -> GraphMol{SDFileAtom,SDFileBond}\n\nRead a SDFile(.sdf or .mol) and parse it into a molecule object. The given argument should be a file input stream, a file path as a string or an iterator that yields each sdfile text lines.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#SDFile-reader-1",
    "page": "I/O",
    "title": "SDFile reader",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"sdfilereader.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdfilewriter-Tuple{IO,Any}",
    "page": "I/O",
    "title": "MolecularGraph.sdfilewriter",
    "category": "method",
    "text": "sdfilewriter(io::IO, mols)\nsdfilewriter(filename::AbstractString, mols)\n\nWrite molecule data to the output stream as a SDFile format file.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#SDFile-writer-1",
    "page": "I/O",
    "title": "SDFile writer",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"sdfilewriter.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/io/#MolecularGraph.smartstomol-Tuple{AbstractString}",
    "page": "I/O",
    "title": "MolecularGraph.smartstomol",
    "category": "method",
    "text": "smartstomol(smarts::AbstractString) -> QueryMol{SmartsAtom,SmartsBond}\n\nParse SMARTS string into QueryMol object.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.smilestomol-Tuple{AbstractString}",
    "page": "I/O",
    "title": "MolecularGraph.smilestomol",
    "category": "method",
    "text": "smilestomol(smiles::AbstractString) -> GraphMol{SmilesAtom,SmilesBond}\n\nParse SMILES string into GraphMol object.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#SMILES/SMARTS-1",
    "page": "I/O",
    "title": "SMILES/SMARTS",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"smarts/base.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
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
    "location": "moleculargraph/draw/#MolecularGraph.atomcolor-Tuple{GraphMol}",
    "page": "Structure drawing",
    "title": "MolecularGraph.atomcolor",
    "category": "method",
    "text": "atomcolor(mol::GraphMol; setting=DRAW_SETTING) -> Vector{Color}\n\nReturn atom colors for molecule 2D drawing\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.chargesign-Tuple{Int64}",
    "page": "Structure drawing",
    "title": "MolecularGraph.chargesign",
    "category": "method",
    "text": "chargesign(charge::Int) -> String\n\nGet a charge sign.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.draw2d!-Tuple{Canvas,MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "Structure drawing",
    "title": "MolecularGraph.draw2d!",
    "category": "method",
    "text": "draw2d!(canvas::Canvas, mol::UndirectedGraph;\n        setting=copy(DRAW_SETTING), recalculate=false)\n\nDraw molecular image to the canvas.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.isatomvisible-Tuple{GraphMol}",
    "page": "Structure drawing",
    "title": "MolecularGraph.isatomvisible",
    "category": "method",
    "text": "isatomvisible(mol::GraphMol; setting=DRAW_SETTING) -> Vector{Bool}\n\nReturn whether the atom is visible in the 2D drawing.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.drawsvg-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,Int64,Int64}",
    "page": "Structure drawing",
    "title": "MolecularGraph.drawsvg",
    "category": "method",
    "text": "drawsvg(mol::GraphMol, width::Int, height::Int)\n\nGenerate molecular structure image as a SVG format string.\n\nwidth and height specifies the size of the image (width and height attribute of svg tag).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/draw/#MolecularGraph.initcanvas!-Tuple{Canvas,GraphMol}",
    "page": "Structure drawing",
    "title": "MolecularGraph.initcanvas!",
    "category": "method",
    "text": "initcanvas!(canvas::Canvas, mol::GraphMol)\n\nMove and adjust the size of the molecule for drawing.\n\n\n\n\n\n"
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
    "page": "Chemical properties",
    "title": "Chemical properties",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/properties/#MolecularGraph.hacceptorcount-Tuple{GraphMol}",
    "page": "Chemical properties",
    "title": "MolecularGraph.hacceptorcount",
    "category": "method",
    "text": "hacceptorcount(mol::GraphMol) -> Int\n\nReturn the number of hydrogen bond acceptors (N, O and F).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.hdonorcount-Tuple{GraphMol}",
    "page": "Chemical properties",
    "title": "MolecularGraph.hdonorcount",
    "category": "method",
    "text": "hdonorcount(mol::GraphMol) -> Int\n\nReturn the number of hydrogen bond donors (O and N attached to hydrogens).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.molweight-Tuple{GraphMol}",
    "page": "Chemical properties",
    "title": "MolecularGraph.molweight",
    "category": "method",
    "text": "molweight(mol::GraphMol; digits=2) -> Float64\n\nReturn standard molecular weight.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.rotatablecount-Tuple{GraphMol}",
    "page": "Chemical properties",
    "title": "MolecularGraph.rotatablecount",
    "category": "method",
    "text": "rotatablecount(mol::GraphMol) -> Int\n\nReturn the number of rotatable bonds.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.wclogp-Tuple{GraphMol}",
    "page": "Chemical properties",
    "title": "MolecularGraph.wclogp",
    "category": "method",
    "text": "wclogp(mol::GraphMol) -> Float64\n\nReturn predicted logP value calculated by using Wildman and Crippen method.\n\nReference\n\nWildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical Parameters by Atomic Contributions. Journal of Chemical Information and Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#Molecular-properties-1",
    "page": "Chemical properties",
    "title": "Molecular properties",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"properties.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/preprocess/#",
    "page": "Preprocessing",
    "title": "Preprocessing",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.allhydrogens-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.allhydrogens",
    "category": "method",
    "text": "allhydrogens(mol::GraphMol) -> Set{Int}\n\nReturn a set of hydrogen nodes.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.canonicalize!-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.canonicalize!",
    "category": "method",
    "text": "canonicalize!(mol::GraphMol)\n\nCanonicalize molecule notation and apply the changes to the molecular property vector.\n\nNeutralize oxo acid, 1-3° ammonium and polarized carbonyls except in the case that polarization is required for aromaticity.\nCanonicalize anions next to triple bonds (ex. [C-][N+]#N -> C=[N+]=[N-])\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.depolarize!-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.depolarize!",
    "category": "method",
    "text": "depolarize!(mol::GraphMol)\n\nDepolarize oxo groups except in the case that polarization is required for aromaticity.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.largestcomponentgraph-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.largestcomponentgraph",
    "category": "method",
    "text": "largestcomponentgraph(mol::GraphMol) -> GraphMol\n\nReturn largest connected component of the molecular graph.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.largestcomponentnodes-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.largestcomponentnodes",
    "category": "method",
    "text": "largestcomponentnodes(mol::GraphMol) -> Set{Int}\n\nReturn a set of nodes in the largest connected component.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.makehydrogensexplicit-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.makehydrogensexplicit",
    "category": "method",
    "text": "makehydrogensexplicit(mol::GraphMol) -> GraphMol\n\nReturn molecule whose hydrogens are fully attached. If option all is set to false, only trivial hydrogens are removed (see trivialhydrogens).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.makehydrogensimplicit-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.makehydrogensimplicit",
    "category": "method",
    "text": "makehydrogensimplicit(mol::GraphMol) -> GraphMol\n\nReturn molecule whose hydrogen nodes are removed. If option all is set to false, only trivial hydrogens are removed (see trivialhydrogens).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.neutralizeacids!-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.neutralizeacids!",
    "category": "method",
    "text": "neutralizeacids!(mol::GraphMol)\n\nNeutralize oxo(thio) acids.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.neutralizeoniums!-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.neutralizeoniums!",
    "category": "method",
    "text": "neutralizeoniums!(mol::GraphMol)\n\nNeutralize 1-3° oniums. Permanently charged quart-oniums are not neutralized.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.triplebondanion!-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.triplebondanion!",
    "category": "method",
    "text": "triplebondanion!(mol::GraphMol)\n\nCanonicalize anions next to triple bonds (ex. [C-][N+]#N -> C=[N+]=[N-]).\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.trivialhydrogens-Tuple{GraphMol}",
    "page": "Preprocessing",
    "title": "MolecularGraph.trivialhydrogens",
    "category": "method",
    "text": "trivialhydrogens(mol::GraphMol) -> Set{Int}\n\nReturn a set of trivial hydrogen nodes (light hydrogens which are uncharged, non-radical, non-stereospecific and attached to organic heavy atoms)\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#Preprocessing-1",
    "page": "Preprocessing",
    "title": "Preprocessing",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"preprocess.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/substructure/#",
    "page": "Substructure match",
    "title": "Substructure match",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/substructure/#MolecularGraph.fastquerymatches-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,QueryMol}",
    "page": "Substructure match",
    "title": "MolecularGraph.fastquerymatches",
    "category": "method",
    "text": "fastquerymatches(mol::UndirectedGraph, query::QueryMol; kwargs...\n    ) -> Dict{Int,Int}\n\nGenerate query match mappings between mol and query. If no match found, return nothing.\n\nThe query should not have any component level expression that means it should not have any dots (.). This is intended for use in functional group detection.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/substructure/#MolecularGraph.isquerymatch-Tuple{Any,Any}",
    "page": "Substructure match",
    "title": "MolecularGraph.isquerymatch",
    "category": "method",
    "text": "isquerymatch(mol, query; kwargs...)\n\nReturn whether mol matches with the query.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/substructure/#MolecularGraph.querymatch-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,QueryMol}",
    "page": "Substructure match",
    "title": "MolecularGraph.querymatch",
    "category": "method",
    "text": "querymatch(mol::UndirectedGraph, query::QueryMol; kwargs...) -> Dict{Int,Int}\n\nGenerate substructure match mappings between mol1 and mol2. If no match found, return nothing.\n\nThis accepts also disconnected single atom but returns only the first match. This function is intended for use in SMARTS query search\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/substructure/#MolecularGraph.structmatches-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "Substructure match",
    "title": "MolecularGraph.structmatches",
    "category": "method",
    "text": "structmatches(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...) -> Iterator\n\nGenerate molecular graph match mappings between mol1 and mol2. If no match found, return nothing.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/substructure/#MolecularGraph.substructmatches-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "Substructure match",
    "title": "MolecularGraph.substructmatches",
    "category": "method",
    "text": "substructmatches(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...) -> Iterator\n\nGenerate substructure match mappings between mol1 and mol2. If no match found, return nothing.\n\nThe mapping is based on only edge induced subgraph isomorphism and therefore it does not care disconnected single atom matches. This function is intended for use in substructure search. If you need single atom SMARTS match (ex. [#16;X2;!R]), see querymatch.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/substructure/#Structure-match-1",
    "page": "Substructure match",
    "title": "Structure match",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"substructure.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/mcs/#",
    "page": "MCS",
    "title": "MCS",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/mcs/#MolecularGraph.mcesmol-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "MCS",
    "title": "MolecularGraph.mcesmol",
    "category": "method",
    "text": "mcesmol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...\n    ) -> Tuple{Dict{Int,Int},Symbol}\n\nCompute maximum common edge induced substructure (MCES) of mol1 and mol2.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/mcs/#MolecularGraph.mcismol-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "MCS",
    "title": "MolecularGraph.mcismol",
    "category": "method",
    "text": "mcismol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...\n    ) -> Tuple{Dict{Int,Int},Symbol}\n\nCompute maximum common induced substructure (MCIS) of mol1 and mol2.\n\nKeyword arguments\n\nconnected(Bool): if true, apply connected MCS constraint.\ntopological(Bool): if true, apply topological constraint.\ndiameter(Int): distance cutoff for topological constraint.\ntolerance(Int): distance mismatch tolerance for topological constraint.\ntimeout(Int): abort calculation and return suboptimal results if the execution\n\ntime has reached the given value (default=60, in seconds).\n\ntargetsize(Int): abort calculation and return suboptimal result so far if the\n\ngiven mcs size achieved.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/mcs/#Structure-match-1",
    "page": "MCS",
    "title": "Structure match",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"mcs.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "moleculargraph/funcgroup/#",
    "page": "Functional group detection",
    "title": "Functional group detection",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/funcgroup/#MolecularGraph.functionalgroupgraph-Tuple{GraphMol}",
    "page": "Functional group detection",
    "title": "MolecularGraph.functionalgroupgraph",
    "category": "method",
    "text": "functionalgroupgraph(mol::GraphMol) -> FunctionalGroupClassifier\n\nGenerate functional group graph that is a directed acyclic graph similar to an ontology graph.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/funcgroup/#Functional-group-detection-1",
    "page": "Functional group detection",
    "title": "Functional group detection",
    "category": "section",
    "text": "Modules = [MolecularGraph]\nPages   = [\"funcgroup.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
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
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.addnode!-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,MolecularGraph.MolecularGraphModel.AbstractNode}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.addnode!",
    "category": "method",
    "text": "addnode!(graph) -> Int\naddnode!(graph, attr) -> Int\n\nAdd new node and return the node index. If the node attribute type is required, specify the node attribute object by node keyword.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgeattr-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph,Int64,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgeattr",
    "category": "method",
    "text": "edgeattr(graph::AbstractGraph, u::Int, v::Int\n    ) -> Union{AbstractEdge,Nothing}\n\nReturn the attribute object of an edge that connects u and v. If not found, return nothing.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgeattr-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgeattr",
    "category": "method",
    "text": "edgeattr(graph::AbstractGraph, i::Int) -> AbstractEdge\n\nReturn the attribute object of edge i.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgeattrs-Tuple{Union{OrderedDiGraph, OrderedGraph, OrderedHyperGraph}}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgeattrs",
    "category": "method",
    "text": "edgeattrs(graph::Union{OrderedGraph,OrderedDiGraph}) -> Vector{AbstractEdge}\n\nReturn graph edge attributes.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgeattrtype-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgeattrtype",
    "category": "method",
    "text": "edgeattrtype(graph) -> Type\n\nReturn edge attribute type.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgecount-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgecount",
    "category": "method",
    "text": "edgecount(graph::AbstractGraph) -> Int\n\nReturn the number of graph edges.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.edgeset-Tuple{Union{OrderedDiGraph, OrderedGraph, OrderedHyperGraph}}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.edgeset",
    "category": "method",
    "text": "edgeset(graph::Union{OrderedGraph,OrderedDiGraph}) -> Set{Int}\n\nReturn the set of edge keys.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.getedge-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.getedge",
    "category": "method",
    "text": "getedge(graph::AbstractGraph, i::Int) -> Tuple{Int,Int}\ngetedge(graph::HyperGraph, i::Int) -> Set{Int}\n\nReturn an edge.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.hasedge-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph,Int64,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.hasedge",
    "category": "method",
    "text": "hasedge(graph::UndirectedGraph, u::Int, v::Int) -> Bool\nhasedge(graph::DirectedGraph, source::Int, target::Int) -> Bool\n\nReturn whether the given two nodes are connected by at least one edge.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.indegree-Tuple{MolecularGraph.MolecularGraphModel.DirectedGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.indegree",
    "category": "method",
    "text": "indegree(graph::DirectedGraph, n::Int) -> Int\n\nReturn the number of inneighbors of the node \'n\'.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.inneighbors-Tuple{MolecularGraph.MolecularGraphModel.DirectedGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.inneighbors",
    "category": "method",
    "text": "inneighbors(graph::DirectedGraph, i::Int) -> Dict{Int,Int}\n\nReturn the mapping of predecessor node keys and in edge keys connected to the given node.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.neighborcount-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.neighborcount",
    "category": "method",
    "text": "neighborcount(graph::AbstractGraph, n::Int) -> Int\ndegree(graph::AbstractGraph, n::Int) -> Int\n\nReturn the number of adjacent nodes of the node \'n\'.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.neighbors-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.neighbors",
    "category": "method",
    "text": "neighbors(graph, i) -> Dict{Int,Int}\n\nReturn the mapping of incident edges and adjacent nodes of node i. If the graph is directed graph, both outneighbors and inneighbors are mapped.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodeattr-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodeattr",
    "category": "method",
    "text": "nodeattr(graph::AbstractGraph, i::Int) -> AbstractNode\n\nReturn the attribute object of node i.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodeattrs-Tuple{Union{OrderedDiGraph, OrderedGraph, OrderedHyperGraph}}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodeattrs",
    "category": "method",
    "text": "nodeattrs(graph::Union{OrderedGraph,OrderedDiGraph}) -> Vector{AbstractNode}\n\nReturn graph node attributes.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodeattrtype-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodeattrtype",
    "category": "method",
    "text": "nodeattrtype(graph) -> Type\n\nReturn node attribute type.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodecount-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodecount",
    "category": "method",
    "text": "nodecount(graph::AbstractGraph) -> Int\n\nReturn the number of graph nodes.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.nodeset-Tuple{Union{OrderedDiGraph, OrderedGraph, OrderedHyperGraph}}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.nodeset",
    "category": "method",
    "text": "nodeset(graph::Union{OrderedGraph,OrderedDiGraph}) -> Set{Int}\n\nReturn the set of node keys.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.outdegree-Tuple{MolecularGraph.MolecularGraphModel.DirectedGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.outdegree",
    "category": "method",
    "text": "outdegree(graph::DirectedGraph, n::Int) -> Int\n\nReturn the number of outneighbors of the node \'n\'.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.outneighbors-Tuple{MolecularGraph.MolecularGraphModel.DirectedGraph,Int64}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.outneighbors",
    "category": "method",
    "text": "outneighbors(graph::DirectedGraph, i::Int) -> Dict{Int,Int}\n\nReturn the mapping of successor node keys and out edge keys connected to the given node.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.setedgeattr!-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph,Int64,MolecularGraph.MolecularGraphModel.AbstractEdge}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.setedgeattr!",
    "category": "method",
    "text": "setedgeattr!(graph::AbstractGraph, i::Int, attr::AbstractNode)\n\nUpdate the edge attribute.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.setnodeattr!-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph,Int64,MolecularGraph.MolecularGraphModel.AbstractNode}",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.setnodeattr!",
    "category": "method",
    "text": "setnodeattr!(graph::AbstractGraph, i::Int, attr::AbstractNode)\n\nUpdate the node attribute.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.unlinkedges",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.unlinkedges",
    "category": "function",
    "text": "unlinkedges(graph, edges)\n\nDelete given edges from the graph.\n\n\n\n\n\n"
},

{
    "location": "graph/interface/#MolecularGraph.MolecularGraphModel.unlinknodes",
    "page": "Interface",
    "title": "MolecularGraph.MolecularGraphModel.unlinknodes",
    "category": "function",
    "text": "unlinknodes(graph, nodes)\n\nDelete given nodes and its incident edges from the graph.\n\n\n\n\n\n"
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
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.circularladder-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.circularladder",
    "category": "method",
    "text": "circularladder(n::Int; mutable=false) -> PlainGraph\n\nGenerate circular ladder graph CL_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.completebipartite-Tuple{Int64,Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.completebipartite",
    "category": "method",
    "text": "completebipartite(m::Int,n::Int; mutable=false) -> PlainGraph\n\nGenerate bipartite graph K_mn.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.completegraph-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.completegraph",
    "category": "method",
    "text": "completegraph(length::Int; mutable=false) -> PlainGraph\n\nGenerate complete graph K_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.cyclegraph-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.cyclegraph",
    "category": "method",
    "text": "cyclegraph(length::Int; mutable=false) -> PlainGraph\n\nGenerate cycle graph C_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.generalizedpetersen-Tuple{Int64,Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.generalizedpetersen",
    "category": "method",
    "text": "generalizedpetersen(n::Int, k::Int; mutable=false) -> PlainGraph\n\nGenerate generalized petersen graph G(nk).\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.laddergraph-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.laddergraph",
    "category": "method",
    "text": "laddergraph(n::Int; mutable=false) -> PlainGraph\n\nGenerate ladder graph L_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.moebiusladder-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.moebiusladder",
    "category": "method",
    "text": "moebiusladder(n::Int; mutable=false) -> PlainGraph\n\nGenerate Möbius ladder graph ML_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.pathgraph-Tuple{Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.pathgraph",
    "category": "method",
    "text": "pathgraph(n::Int; mutable=false) -> PlainGraph\n\nGenerate path graph P_n.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#MolecularGraph.MolecularGraphModel.squaregrid-Tuple{Int64,Int64}",
    "page": "Generator",
    "title": "MolecularGraph.MolecularGraphModel.squaregrid",
    "category": "method",
    "text": "squaregrid(m::Int,n::Int; mutable=false) -> PlainGraph\n\nGenerate m times n square grid graph.\n\nUse cartesianproduct for higher dimensional grid graphs.\n\n\n\n\n\n"
},

{
    "location": "graph/generator/#Graph-generator-1",
    "page": "Generator",
    "title": "Graph generator",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/generator.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/traversal/#",
    "page": "Traversal",
    "title": "Traversal",
    "category": "page",
    "text": ""
},

{
    "location": "graph/traversal/#Graph-traversal-1",
    "page": "Traversal",
    "title": "Graph traversal",
    "category": "section",
    "text": ""
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.diameter-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.diameter",
    "category": "method",
    "text": "diameter(graph::AbstractGraph) -> Int\n\nCompute the diameter of the graph (the largest eccentricity of any nodes).\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.distance-Tuple{Function,Any,Any,Any}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.distance",
    "category": "method",
    "text": "distance(graph::AbstractGraph, source::Int, target::Int) -> Int\nreversedistance(graph::DirectedGraph, source::Int, target::Int) -> Int\n\nCompute the distance (shortest path length) from source to target. If the nodes are not reachable each other, the value will be nothing.\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.distancematrix-Tuple{Function,Any}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.distancematrix",
    "category": "method",
    "text": "distancematrix(graph::OrderedGraph) -> Matrix{Float64}\n\nGenerate the distance matrix of the graph.\n\nNote that the type of the generated matrix will be Float64. If the nodes are not reachable each other, the distance value will be Inf.\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.eccentricity-Tuple{Function,Any,Any}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.eccentricity",
    "category": "method",
    "text": "eccentricity(graph::UndirectedGraph, v::Int) -> Int\n\nCompute the eccentricity of the graph (the largest distance between v and any other nodes).\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.isreachable-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,Int64,Int64}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.isreachable",
    "category": "method",
    "text": "reachablenodes(graph::AbstractGraph, u::Int, v::Int) -> Bool\nreversereachablenodes(graph::DirectedGraph, u::Int, v::Int) -> Bool\n\nReturn whether the node v is reachable from u.\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.longestshortestpathnodes-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.longestshortestpathnodes",
    "category": "method",
    "text": "longestshortestpathnodes(graph::UndirectedGraph) -> Vector{Int}\n\nCompute the longest shortest path in the graph (a path between two arbitrary peripheral nodes) as a vector of nodes that starts with one of the peripheral node and ends with the other side.\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.reachablenodes-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,Int64}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.reachablenodes",
    "category": "method",
    "text": "reachablenodes(graph::AbstractGraph, node::Int) -> Set{Int}\nreversereachablenodes(graph::DirectedGraph, node::Int) -> Set{Int}\n\nReturn the set of reachable nodes from node.\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.shortestpathedges-Tuple{Function,Any,Any,Any}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.shortestpathedges",
    "category": "method",
    "text": "shortestpathedges(graph::UndirectedGraph, u::Int, v::Int) -> Vector{Int}\n\nCompute the shortest path between u and v as a vector of the edges that forms the path. Return nothing if not reachable.\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#MolecularGraph.MolecularGraphModel.shortestpathnodes-Tuple{Function,Any,Any,Any}",
    "page": "Traversal",
    "title": "MolecularGraph.MolecularGraphModel.shortestpathnodes",
    "category": "method",
    "text": "shortestpathnodes(graph::UndirectedGraph, u::Int, v::Int) -> Vector{Int}\n\nCompute the shortest path between u and v as a vector of the nodes that forms the path. Return nothing if not reachable.\n\n\n\n\n\n"
},

{
    "location": "graph/traversal/#Shortest-path-1",
    "page": "Traversal",
    "title": "Shortest path",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/shortestpath.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/topology/#",
    "page": "Topology",
    "title": "Topology",
    "category": "page",
    "text": ""
},

{
    "location": "graph/topology/#Graph-topology-1",
    "page": "Topology",
    "title": "Graph topology",
    "category": "section",
    "text": ""
},

{
    "location": "graph/topology/#Cycle-1",
    "page": "Topology",
    "title": "Cycle",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/cycle.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/topology/#Connectivity-1",
    "page": "Topology",
    "title": "Connectivity",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/connectivity.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/topology/#Planarity-1",
    "page": "Topology",
    "title": "Planarity",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/planarity.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/operation/#",
    "page": "Operations",
    "title": "Operations",
    "category": "page",
    "text": ""
},

{
    "location": "graph/operation/#Graph-operations-1",
    "page": "Operations",
    "title": "Graph operations",
    "category": "section",
    "text": ""
},

{
    "location": "graph/operation/#MolecularGraph.MolecularGraphModel.disjointunion!-Union{Tuple{T}, Tuple{T,T,Vararg{T,N} where N}} where T<:MolecularGraph.MolecularGraphModel.OrderedGraph",
    "page": "Operations",
    "title": "MolecularGraph.MolecularGraphModel.disjointunion!",
    "category": "method",
    "text": "disjointunion!(g1::T, g2::T, G::T...) where {T<:OrderedGraph} -> T\n\nGenerate disjoint union graph of given graphs. g1 will be overwritten by the union graph. Unlike non-destructive disjointunion, g1 does not retain any information about other given graphs but a bit faster.\n\n\n\n\n\n"
},

{
    "location": "graph/operation/#MolecularGraph.MolecularGraphModel.disjointunion-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,MolecularGraph.MolecularGraphModel.UndirectedGraph,Vararg{MolecularGraph.MolecularGraphModel.UndirectedGraph,N} where N}",
    "page": "Operations",
    "title": "MolecularGraph.MolecularGraphModel.disjointunion",
    "category": "method",
    "text": "disjointunion(g1::UndirectedGraph, g2::UndirectedGraph,\n    G::UndirectedGraph...) -> DisjointUnionGraph\n\nGenerate disjoint union graph of given graphs. The new graph with type DisjointUnionGraph retains mapping to the original graphs as nodes and edges attributes.\n\n\n\n\n\n"
},

{
    "location": "graph/operation/#Disjoint-union-1",
    "page": "Operations",
    "title": "Disjoint union",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/disjointunion.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/operation/#MolecularGraph.MolecularGraphModel.linegraph-Tuple{MolecularGraph.MolecularGraphModel.AbstractGraph}",
    "page": "Operations",
    "title": "MolecularGraph.MolecularGraphModel.linegraph",
    "category": "method",
    "text": "linegraph(G::AbstractGraph) -> LineGraph\n\nGenerate line graph.\n\n\n\n\n\n"
},

{
    "location": "graph/operation/#Line-graph-1",
    "page": "Operations",
    "title": "Line graph",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/linegraph.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/operation/#MolecularGraph.MolecularGraphModel.cartesianproduct-Tuple{MolecularGraph.MolecularGraphModel.OrderedGraph,MolecularGraph.MolecularGraphModel.OrderedGraph}",
    "page": "Operations",
    "title": "MolecularGraph.MolecularGraphModel.cartesianproduct",
    "category": "method",
    "text": "cartesianproduct(G::OrderedGraph, H::OrderedGraph) -> CartesianProduct\n\nReturn the cartesian product of graphs G and H.\n\n\n\n\n\n"
},

{
    "location": "graph/operation/#MolecularGraph.MolecularGraphModel.modularproduct-Tuple{MolecularGraph.MolecularGraphModel.OrderedGraph,MolecularGraph.MolecularGraphModel.OrderedGraph}",
    "page": "Operations",
    "title": "MolecularGraph.MolecularGraphModel.modularproduct",
    "category": "method",
    "text": "modularproduct(G::OrderedGraph, H::OrderedGraph) -> ModularProduct\n\nReturn the modular product of graphs G and H.\n\n\n\n\n\n"
},

{
    "location": "graph/operation/#Product-of-graphs-1",
    "page": "Operations",
    "title": "Product of graphs",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/product.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
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
    "location": "graph/clique/#MolecularGraph.MolecularGraphModel.maximalcliques-Union{Tuple{T}, Tuple{T}} where T<:MolecularGraph.MolecularGraphModel.UndirectedGraph",
    "page": "Clique",
    "title": "MolecularGraph.MolecularGraphModel.maximalcliques",
    "category": "method",
    "text": "maximalcliques(graph::UndirectedGraph; kwargs...\n    ) -> Tuple{Vector{Set{Int}}, Symbol}\n\nReturn maximal cliques.\n\nReference\n\nTomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time complexity for generating all maximal cliques and computational experiments. Theoretical Computer Science, 363(1), 28–42. https://doi.org/10.1016/J.TCS.2006.06.015\nCazals, F., & Karande, C. (2008). A note on the problem of reporting maximal cliques. Theoretical Computer Science, 407(1–3), 564–568. https://doi.org/10.1016/j.tcs.2008.05.010\n\n\n\n\n\n"
},

{
    "location": "graph/clique/#MolecularGraph.MolecularGraphModel.maximalconncliques-Union{Tuple{T}, Tuple{T}} where T<:MolecularGraph.MolecularGraphModel.UndirectedGraph",
    "page": "Clique",
    "title": "MolecularGraph.MolecularGraphModel.maximalconncliques",
    "category": "method",
    "text": "maximalconncliques(graph::ModularProduct; kwargs...\n    ) -> Tuple{Vector{Set{Int}}, Symbol}\n\nReturn maximal connected cliques.\n\nReference\n\nCazals, F., & Karande, C. (2005). An algorithm for reporting maximal c-cliques. Theoretical Computer Science, 349(3), 484–490. https://doi.org/10.1016/j.tcs.2005.09.038\n\n\n\n\n\n"
},

{
    "location": "graph/clique/#MolecularGraph.MolecularGraphModel.maximumclique-Union{Tuple{T}, Tuple{T}} where T<:MolecularGraph.MolecularGraphModel.UndirectedGraph",
    "page": "Clique",
    "title": "MolecularGraph.MolecularGraphModel.maximumclique",
    "category": "method",
    "text": "maximumclique(graph::UndirectedGraph; kwargs...) -> Tuple{Set{Int}, Symbol}\n\nReturn a maximum clique.\n\n\n\n\n\n"
},

{
    "location": "graph/clique/#MolecularGraph.MolecularGraphModel.maximumconnclique-Union{Tuple{T}, Tuple{T}} where T<:MolecularGraph.MolecularGraphModel.ModularProduct",
    "page": "Clique",
    "title": "MolecularGraph.MolecularGraphModel.maximumconnclique",
    "category": "method",
    "text": "maximumconnclique(graph::ModularProduct; kwargs...\n    ) -> Tuple{Set{Int}, Symbol}\n\nReturn a maximum connected clique.\n\n\n\n\n\n"
},

{
    "location": "graph/clique/#Clique-2",
    "page": "Clique",
    "title": "Clique",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/clique.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/isomorphism/#",
    "page": "Isomorphism",
    "title": "Isomorphism",
    "category": "page",
    "text": ""
},

{
    "location": "graph/isomorphism/#Graph-isomorphism-1",
    "page": "Isomorphism",
    "title": "Graph isomorphism",
    "category": "section",
    "text": ""
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.edgesubgraphmatch-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.edgesubgraphmatch",
    "category": "method",
    "text": "edgesubgraphmatch(\n    G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}\n\nReturn a edge induced subgraph isomorphism mapping between G and H. If no match found, return nothing.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.edgesubgraphmatches-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.edgesubgraphmatches",
    "category": "method",
    "text": "edgesubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator\n\nGenerate edge induced subgraph isomorphism mappings between G and H.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.graphmatch-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.graphmatch",
    "category": "method",
    "text": "graphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}\n\nReturn an isomorphism mapping between G and H. If no match found, return nothing.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.graphmatches-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.graphmatches",
    "category": "method",
    "text": "graphmatches(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator\n\nGenerate isomorphism mappings between G and H. If no match found, return nothing.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.isedgesubgraphmatch-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.isedgesubgraphmatch",
    "category": "method",
    "text": "isedgesubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool\n\nReturn true if a node induced subgraph of G and H are isomorphic.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.isgraphmatch-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.isgraphmatch",
    "category": "method",
    "text": "isgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool\n\nReturn true if G and H are isomorphic.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.issubgraphmatch-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.issubgraphmatch",
    "category": "method",
    "text": "issubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool\n\nReturn true if a node induced subgraph of G and H are isomorphic.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.subgraphmatch-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.subgraphmatch",
    "category": "method",
    "text": "subgraphmatch(\n    G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}\n\nReturn a subgraph isomorphism mapping between G and H. If no match found, return nothing.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.subgraphmatches-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.subgraphmatches",
    "category": "method",
    "text": "subgraphmatches(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator\n\nGenerate subgraph isomorphism mappings between G and H.\n\nKeyword arguments\n\nnodematcher(Function): node matcher function that takes two node indices as\n\narguments.\n\nedgematcher(Function): edge matcher function that takes two edge indices as\n\narguments.\n\nmandatory(Dict{Int,Int}): mandatory node matches (available for only VF2)\nforbidden(Dict{Int,Int}):   forbidden node matches (available for only VF2)\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#Subgraph-match-(VF2-algorithm)-1",
    "page": "Isomorphism",
    "title": "Subgraph match (VF2 algorithm)",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/isomorphism/vf2.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.findmces-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.findmces",
    "category": "method",
    "text": "findmces(G::UndirectedGraph, H::UndirectedGraph; kwargs...\n    ) -> Tuple{Dict{Int,Int},Symbol}\n\nCompute maximum common edge induced subgraph between G and H.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.findmcis-Tuple{MolecularGraph.MolecularGraphModel.UndirectedGraph,MolecularGraph.MolecularGraphModel.UndirectedGraph}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.findmcis",
    "category": "method",
    "text": "findmcis(G::UndirectedGraph, H::UndirectedGraph; kwargs...\n    ) -> Tuple{Dict{Int,Int},Symbol}\n\nCompute maximum common induced subgraph between G and H.\n\nKeyword arguments:\n\nconnected(Bool): if true, apply connected MCS constraint.\ntopological(Bool): if true, apply topological constraint.\ndiameter(Int): distance cutoff for topological constraint.\ntolerance(Int): distance mismatch tolerance for topological constraint.\ntimeout(Int): abort calculation and return suboptimal results so far if the\n\nexecution time has reached the given value (default=60, in seconds).\n\ntargetsize(Int): abort calculation and return suboptimal result so far if the\n\ngiven mcs size achieved.\n\nnodematcher(Function): node matcher function that takes two node indices as\n\narguments.\n\nedgematcher(Function): edge matcher function that takes two edge indices as\n\narguments.\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.mcessize-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.mcessize",
    "category": "method",
    "text": "mcessize(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Int\n\nReturn the maximum common edge induced subgraph size (number of edges).\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#MolecularGraph.MolecularGraphModel.mcissize-Tuple{Any,Any}",
    "page": "Isomorphism",
    "title": "MolecularGraph.MolecularGraphModel.mcissize",
    "category": "method",
    "text": "mcissize(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Int\n\nReturn the maximum common induced subgraph size (number of nodes).\n\n\n\n\n\n"
},

{
    "location": "graph/isomorphism/#Maximum-common-subgraph-(Clique-detection-based-algorithm)-1",
    "page": "Isomorphism",
    "title": "Maximum common subgraph (Clique detection based algorithm)",
    "category": "section",
    "text": "Modules = [MolecularGraphModel]\nPages   = [\"graph/isomorphism/cliquemcs.jl\"]\nPrivate = false\nOrder   = [:constant, :function, :type]"
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
    "location": "design/#Graph-type-hierarchy-1",
    "page": "Design of molecular graph models",
    "title": "Graph type hierarchy",
    "category": "section",
    "text": "(AbstractGraph)\n(UndirectedGraph)\n(OrderedGraph)\nPlainGraph\nImmutablePlainGraph\nGraphMol{Atom,Bond}\nSDFile (Alias of GraphMol{SDFileAtom,SDFileBond})\nSMILES (Alias of GraphMol{SmilesAtom,SmilesBond})\nQueryMol{QueryAtom,QueryBond}\nSMARTS (Alias of QueryMol{SmartsAtom,SmartsBond})\nLineGraph\nCartesianProductGraph\nModularProductGraph\nSubgraphView{UndirectedGraph}\n(DirectedGraph)\n(OrderedDiGraph)\nPlainDiGraph\nFunctionalGroupClassGraph\nDiSubgraphView{DirectedGraph}\nHyperGraph?"
},

{
    "location": "design/#AbstractGraph-methods-1",
    "page": "Design of molecular graph models",
    "title": "AbstractGraph methods",
    "category": "section",
    "text": "getnode, getedge, hasedge\nneighbors and its derivatives\nnodecount\nedgecount\nnodeset\nedgeset"
},

{
    "location": "design/#DirectedGraph-methods-1",
    "page": "Design of molecular graph models",
    "title": "DirectedGraph methods",
    "category": "section",
    "text": "outneighbors and inneighbors"
},

{
    "location": "design/#OrderedGraph-methods-1",
    "page": "Design of molecular graph models",
    "title": "OrderedGraph methods",
    "category": "section",
    "text": "nodesiter\nedgesiter\nnodeattrs\nedgeattrs"
},

{
    "location": "design/#OrderedGraph-1",
    "page": "Design of molecular graph models",
    "title": "OrderedGraph",
    "category": "section",
    "text": "OrderedGraph consists of vectors of neighborhood map (incident edge => adjacent node) and edge (tuple of node index pair) vector."
},

{
    "location": "design/#QueryMol-1",
    "page": "Design of molecular graph models",
    "title": "QueryMol",
    "category": "section",
    "text": "QueryMol consists of QueryAtom and QueryBond that represent molecular query (ex. atom symbol is \'O\' and charge is -1, bond order is 1 and not in rings, ...). This type of objects typically built from SMARTS query."
},

{
    "location": "design/#Node-type-hierarchy-1",
    "page": "Design of molecular graph models",
    "title": "Node type hierarchy",
    "category": "section",
    "text": "(AbstractNode)\n(Atom)\nSDFileAtom\nSmilesAtom\n(QueryAtom)\nSmartsAtom"
},

{
    "location": "design/#Edge-type-hierarchy-1",
    "page": "Design of molecular graph models",
    "title": "Edge type hierarchy",
    "category": "section",
    "text": "(AbstractEdge)\n(UndirectedEdge)\nEdge\n(Bond)\nSDFileBond\nSmilesBond\n(QueryBond)\nSmartsBond\n(DirectedEdge)\nArrow"
},

]}
