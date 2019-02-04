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
    "page": "Molecule I/O",
    "title": "Molecule I/O",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/io/#Molecule-I/O-1",
    "page": "Molecule I/O",
    "title": "Molecule I/O",
    "category": "section",
    "text": ""
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdfilereader-Tuple{IO}",
    "page": "Molecule I/O",
    "title": "MolecularGraph.sdfilereader",
    "category": "method",
    "text": "sdfilereader(file::IO)\n\nRead SDFile data from input stream and return a lazy iterator which yields molecule objects.\n\nsdfilereader does not stop and raise errors when an erroneous or incompatible SDFile block is read but produces an error message and yields an empty molecule. If this behavior is not desirable, you can use the customized supplier function instead of default supplier nohaltsupplier\n\nfunction customsupplier()\n    mol = try\n        parse(SDFile, block)\n    catch e\n        throw(ErrorException(\"incompatible molecule found, aborting...\"))\n    end\n    return defaultpostprocess(mol)\nend\n\nfunction sdfilereader(file::IO)\n    return SDFileReader(eachline(file), customsupplier)\nend\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdfilereader-Tuple{AbstractString}",
    "page": "Molecule I/O",
    "title": "MolecularGraph.sdfilereader",
    "category": "method",
    "text": "sdfilereader(path::AbstractString)\n\nRead a SDFile and return a lazy iterator which yields molecule objects.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdftomol-Tuple{IO}",
    "page": "Molecule I/O",
    "title": "MolecularGraph.sdftomol",
    "category": "method",
    "text": "sdftomol(file::IO)\n\nRead a SDFile mol block from the input stream and parse it into a molecule object.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#MolecularGraph.sdftomol-Tuple{AbstractString}",
    "page": "Molecule I/O",
    "title": "MolecularGraph.sdftomol",
    "category": "method",
    "text": "sdftomol(path::AbstractString)\n\nRead a SDFile and parse it into a molecule object. Single mol block files without optional information are often provided as a .mol file.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#Base.parse-Tuple{Type{GeneralMapMol{SDFileAtom,SDFileBond}},Any}",
    "page": "Molecule I/O",
    "title": "Base.parse",
    "category": "method",
    "text": "parse(::Type{SDFile}, lines)\n\nParse lines of a SDFile mol block data into a molecule object.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/io/#SDFile-1",
    "page": "Molecule I/O",
    "title": "SDFile",
    "category": "section",
    "text": "sdfilereader(file::IO)\nsdfilereader(path::AbstractString)\nsdftomol(file::IO)\nsdftomol(path::AbstractString)\nparse(::Type{SDFile}, lines)"
},

{
    "location": "moleculargraph/properties/#",
    "page": "Basic chemical properties",
    "title": "Basic chemical properties",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/properties/#MolecularGraph.molweight-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.molweight",
    "category": "method",
    "text": "molweight(mol::VectorMol; digits=2) -> Float64\n\nReturn standard molecular weight.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.H_acceptor_count-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.H_acceptor_count",
    "category": "method",
    "text": "H_acceptor_count(mol::VectorMol) -> Int\n\nReturn the number of hydrogen bond acceptors (N, O and F).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.H_donor_count-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.H_donor_count",
    "category": "method",
    "text": "H_donor_count(mol::VectorMol) -> Int\n\nReturn the number of hydrogen bond donors (O and N attached to hydrogens).\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.wclogp-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.wclogp",
    "category": "method",
    "text": "wclogp(mol::VectorMol) -> Float64\n\nReturn predicted logP value calculated by using Wildman and Crippen method.\n\nReference\n\nWildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical Parameters by Atomic Contributions. Journal of Chemical Information and Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#MolecularGraph.rotatable_count-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Basic chemical properties",
    "title": "MolecularGraph.rotatable_count",
    "category": "method",
    "text": "rotatable_count(mol::VectorMol) -> Int\n\nReturn the number of rotatable bonds.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/properties/#Molecular-properties-1",
    "page": "Basic chemical properties",
    "title": "Molecular properties",
    "category": "section",
    "text": "molweight(mol::VectorMol; digits=2)\nH_acceptor_count(mol::VectorMol)\nH_donor_count(mol::VectorMol)\nwclogp(mol::VectorMol; digits=2)\nrotatable_count(mol::VectorMol)"
},

{
    "location": "moleculargraph/preprocess/#",
    "page": "Preprocessing",
    "title": "Preprocessing",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.trivialhydrogens-Tuple{MolGraph}",
    "page": "Preprocessing",
    "title": "MolecularGraph.trivialhydrogens",
    "category": "method",
    "text": "trivialhydrogens(mol::MolGraph) -> Set{Int}\n\nReturn a set of trivial hydrogen nodes (light hydrogen which is uncharged, non-radical, non-stereospecific and attached to organic heavy atoms)\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.allhydrogens-Tuple{MolGraph}",
    "page": "Preprocessing",
    "title": "MolecularGraph.allhydrogens",
    "category": "method",
    "text": "allhydrogens(mol::MolGraph) -> Set{Int}\n\nReturn a set of hydrogen nodes.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.largestcomponent-Tuple{MolGraph}",
    "page": "Preprocessing",
    "title": "MolecularGraph.largestcomponent",
    "category": "method",
    "text": "largestcomponent(mol::MolGraph) -> Set{Int}\n\nReturn a set of nodes in the largest connected component.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.neutralize_acids!-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Preprocessing",
    "title": "MolecularGraph.neutralize_acids!",
    "category": "method",
    "text": "neutralize_acids!(mol::VectorMol)\n\nNeutralize oxo(thio) acids.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.neutralize_oniums!-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Preprocessing",
    "title": "MolecularGraph.neutralize_oniums!",
    "category": "method",
    "text": "neutralize_oniums!(mol::VectorMol)\n\nNeutralize 1-3° oniums. Permanently charged quart-oniums are not neutralized.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.depolarize!-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Preprocessing",
    "title": "MolecularGraph.depolarize!",
    "category": "method",
    "text": "depolarize!(mol::VectorMol)\n\nDepolarize oxo groups except in the case that polarization is required for aromaticity.\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.triplebond_anion!-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Preprocessing",
    "title": "MolecularGraph.triplebond_anion!",
    "category": "method",
    "text": "triplebond_anion!(mol::VectorMol)\n\nCanonicalize anions next to triple bonds (ex. [C-][N+]#N -> C=[N+]=[N-]).\n\nNote that this function edits Atom object fields directly. The molecular property vector needs recalculation to apply the changes. see canonicalize!.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#MolecularGraph.canonicalize!-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Preprocessing",
    "title": "MolecularGraph.canonicalize!",
    "category": "method",
    "text": "canonicalize!(mol::VectorMol)\n\nCanonicalize molecule notation and apply the changes to the molecular property vector.\n\nNeutralize oxo acid, 1-3° ammonium and polarized carbonyls except in the case that polarization is required for aromaticity.\nCanonicalize anions next to triple bonds (ex. [C-][N+]#N -> C=[N+]=[N-])\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/preprocess/#Preprocessing-1",
    "page": "Preprocessing",
    "title": "Preprocessing",
    "category": "section",
    "text": "trivialhydrogens(mol::MolGraph)\nallhydrogens(mol::MolGraph)\nlargestcomponent(mol::MolGraph)\nneutralize_acids!(mol::VectorMol)\nneutralize_oniums!(mol::VectorMol)\ndepolarize!(mol::VectorMol)\ntriplebond_anion!(mol::VectorMol)\ncanonicalize!(mol::VectorMol)"
},

{
    "location": "moleculargraph/structure/#",
    "page": "Structure match",
    "title": "Structure match",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/structure/#Structure-match-1",
    "page": "Structure match",
    "title": "Structure match",
    "category": "section",
    "text": ""
},

{
    "location": "moleculargraph/structure/#MolecularGraph.is_identical-Tuple{Union{VectorMolGraph, VectorMolView},Union{VectorMolGraph, VectorMolView}}",
    "page": "Structure match",
    "title": "MolecularGraph.is_identical",
    "category": "method",
    "text": "is_identical(mol1::VectorMol, mol2::VectorMol)\n\nReturn whether mol1 and mol2 are identical in chemical structure.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/structure/#Identity-1",
    "page": "Structure match",
    "title": "Identity",
    "category": "section",
    "text": "is_identical(mol1::VectorMol, mol2::VectorMol)"
},

{
    "location": "moleculargraph/structure/#MolecularGraph.is_substruct-Tuple{Union{VectorMolGraph, VectorMolView},Union{VectorMolGraph, VectorMolView}}",
    "page": "Structure match",
    "title": "MolecularGraph.is_substruct",
    "category": "method",
    "text": "is_substruct(mol1::VectorMol, mol2::VectorMol)\n\nReturn whether mol1 is a substructure of mol2.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/structure/#MolecularGraph.is_superstruct-Tuple{Any,Any}",
    "page": "Structure match",
    "title": "MolecularGraph.is_superstruct",
    "category": "method",
    "text": "is_superstruct(mol1, mol2)\n\nReturn whether mol1 is a superstructure of mol2.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/structure/#Sub/Superstructure-1",
    "page": "Structure match",
    "title": "Sub/Superstructure",
    "category": "section",
    "text": "is_substruct(mol1::VectorMol, mol2::VectorMol)\nis_superstruct(mol1, mol2)"
},

{
    "location": "moleculargraph/structure/#MolecularGraph.is_querymatch-Tuple{Any,Any}",
    "page": "Structure match",
    "title": "MolecularGraph.is_querymatch",
    "category": "method",
    "text": "is_querymatch(mol, query; kwargs...)\n\nReturn whether mol matches with the query.\n\n\n\n\n\n"
},

{
    "location": "moleculargraph/structure/#Molecular-query-1",
    "page": "Structure match",
    "title": "Molecular query",
    "category": "section",
    "text": "is_querymatch(mol, query; kwargs...)"
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
    "text": "CurrentModule = Graph"
},

{
    "location": "graph/clique/#MolecularGraph.Graph.maxclique-Tuple{Union{UndirectedGraph, UndirectedGraphView}}",
    "page": "Clique",
    "title": "MolecularGraph.Graph.maxclique",
    "category": "method",
    "text": "maxclique(graph::UDGraph; kwargs...) -> Set{Int}\n\nCompute maximum clique of the graph. For details, see maximalcliques.\n\n\n\n\n\n"
},

{
    "location": "graph/clique/#MolecularGraph.Graph.maximalcliques-Tuple{Union{UndirectedGraph, UndirectedGraphView}}",
    "page": "Clique",
    "title": "MolecularGraph.Graph.maximalcliques",
    "category": "method",
    "text": "maximalcliques(graph::UDGraph; kwargs...)\n\nReturn Channel which generates maximal cliques of the graph. Each cliques are represented as a Set of member nodes.\n\nReference\n\nTomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time complexity for generating all maximal cliques and computational experiments. Theoretical Computer Science, 363(1), 28–42. https://doi.org/10.1016/J.TCS.2006.06.015\nCazals, F., & Karande, C. (2005). An algorithm for reporting maximal c-cliques. Theoretical Computer Science, 349(3), 484–490. https://doi.org/10.1016/j.tcs.2005.09.038\n\n\n\n\n\n"
},

{
    "location": "graph/clique/#Clique-2",
    "page": "Clique",
    "title": "Clique",
    "category": "section",
    "text": "maxclique(graph::UDGraph; kwargs...)\nmaximalcliques(graph::UDGraph; kwargs...)"
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
    "location": "moleculargraph/#",
    "page": "Design of molecular graph models",
    "title": "Design of molecular graph models",
    "category": "page",
    "text": ""
},

{
    "location": "moleculargraph/#Design-of-molecular-graph-models-1",
    "page": "Design of molecular graph models",
    "title": "Design of molecular graph models",
    "category": "section",
    "text": ""
},

{
    "location": "moleculargraph/#Molecule-1",
    "page": "Design of molecular graph models",
    "title": "Molecule",
    "category": "section",
    "text": "MolGraph\nMapMolGraph\nGeneralMapMol{A<:Atom, B<:Bond}\nSDFile (Alias of GeneralMapMol{SDFileAtom, SDFileBond})\nSMILES (Alias of GeneralMapMol{SmilesAtom, SmilesBond})\nQueryMolGraph\nConnectedQueryMol{A<:QueryAtom, B<:QueryBond}\nConnectedSMARTS (Alias of ConnectedQueryMol{SmartsAtom, SmartsBond})\nDisconnectedQueryMol{A<:QueryAtom, B<:QueryBond}\nSMARTS (Alias of DisconnectedQueryMol{SmartsAtom, SmartsBond})\nVectorMolGraph\nGeneralVectorMol{A<:Atom, B<:Bond}\nMapMolView\nQueryMolView\nVectorMolView"
},

{
    "location": "moleculargraph/#AbstractMol-1",
    "page": "Design of molecular graph models",
    "title": "AbstractMol",
    "category": "section",
    "text": "MolGraph is the base type of all molecular objects. Fundamental graph element accessor functions and counting functions are available in all molecular model objects."
},

{
    "location": "moleculargraph/#Methods-1",
    "page": "Design of molecular graph models",
    "title": "Methods",
    "category": "section",
    "text": "getatom\ngetbond\nneighbors\nneighborcount (or degree)\natomcount\nbondcount"
},

{
    "location": "moleculargraph/#MapMol-1",
    "page": "Design of molecular graph models",
    "title": "MapMol",
    "category": "section",
    "text": "MapMol is used as a molecular model builder for general purpose.This type inherits AbstractMapMol, a molecular graph model which have map(dict)-based structure. The map-based molecular graph can insert and delete elements (atoms and bonds). This can be easily converted to VectorMol object by using vectormol method"
},

{
    "location": "moleculargraph/#Methods-2",
    "page": "Design of molecular graph models",
    "title": "Methods",
    "category": "section",
    "text": "updateatom!\nupdatebond!\nunlinkatom!\nunlinkbond!\nvectormol"
},

{
    "location": "moleculargraph/#QueryMol-1",
    "page": "Design of molecular graph models",
    "title": "QueryMol",
    "category": "section",
    "text": "QueryMol consists of QueryAtoms and QueryBonds that represents molecular query (ex. atom symbol is \'O\' and charge is -1, bond order is 1 and not in rings, ...). This type of objects typically built from SMARTS query."
},

{
    "location": "moleculargraph/#Methods-3",
    "page": "Design of molecular graph models",
    "title": "Methods",
    "category": "section",
    "text": "updateatom!\nupdatebond!\nunlinkatom!\nunlinkbond!"
},

{
    "location": "moleculargraph/#VectorMol-1",
    "page": "Design of molecular graph models",
    "title": "VectorMol",
    "category": "section",
    "text": "VectorMol is vector(array)-based molecular model which is specialized for element-wise fast computation of molecular properties.VectorMol can iterate over atom properties faster than MapMoland can store calculated molecular properties and annotation arrays which are suitable for vector computation. On the other hand, VectorMol does not have abilities to modify its graph structure (adding or removing elements). VectorMol can be converted to MapMol but the calculated properties and annotations will be lost."
},

{
    "location": "moleculargraph/#Common-methods-1",
    "page": "Design of molecular graph models",
    "title": "Common methods",
    "category": "section",
    "text": "mapmol"
},

{
    "location": "moleculargraph/#Atom-1",
    "page": "Design of molecular graph models",
    "title": "Atom",
    "category": "section",
    "text": "AbstractNode\nAbstractAtom\nAtom\nSDFileAtom\nSmilesAtom\nQueryAtom\nSmartsAtom"
},

{
    "location": "moleculargraph/#Bond-1",
    "page": "Design of molecular graph models",
    "title": "Bond",
    "category": "section",
    "text": "AbstractEdge\nAbstractBond\nBond\nSDFileBond\nSmilesBond\nQueryBond\nSmartsBond"
},

]}
