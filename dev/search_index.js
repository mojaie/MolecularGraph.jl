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
    "text": ""
},

{
    "location": "#Basic-information-1",
    "page": "Home",
    "title": "Basic information",
    "category": "section",
    "text": "README.md on GitHub\nMolecular graph modeling library written in Julia\nAnd some chemoinformatics analysis tools\nMIT License"
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
    "location": "api/properties/#",
    "page": "Molecular properties",
    "title": "Molecular properties",
    "category": "page",
    "text": ""
},

{
    "location": "api/properties/#MolecularGraph.molweight-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Molecular properties",
    "title": "MolecularGraph.molweight",
    "category": "method",
    "text": "molweight(mol::VectorMol; digits=2) -> Float64\n\nReturn standard molecular weight.\n\n\n\n\n\n"
},

{
    "location": "api/properties/#MolecularGraph.H_acceptor_count-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Molecular properties",
    "title": "MolecularGraph.H_acceptor_count",
    "category": "method",
    "text": "H_acceptor_count(mol::VectorMol) -> Int\n\nReturn the number of hydrogen bond acceptors (N, O and F).\n\n\n\n\n\n"
},

{
    "location": "api/properties/#MolecularGraph.H_donor_count-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Molecular properties",
    "title": "MolecularGraph.H_donor_count",
    "category": "method",
    "text": "H_donor_count(mol::VectorMol) -> Int\n\nReturn the number of hydrogen bond donors (O and N attached to hydrogens).\n\n\n\n\n\n"
},

{
    "location": "api/properties/#MolecularGraph.wclogp-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Molecular properties",
    "title": "MolecularGraph.wclogp",
    "category": "method",
    "text": "wclogp(mol::VectorMol) -> Float64\n\nReturn predicted logP value calculated by using Wildman and Crippen method.\n\nReference\n\nWildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical\n\nParameters by Atomic Contributions. Journal of Chemical Information and Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l\n\n\n\n\n\n"
},

{
    "location": "api/properties/#MolecularGraph.rotatable_count-Tuple{Union{VectorMolGraph, VectorMolView}}",
    "page": "Molecular properties",
    "title": "MolecularGraph.rotatable_count",
    "category": "method",
    "text": "rotatable_count(mol::VectorMol) -> Int\n\nReturn the number of rotatable bonds.\n\n\n\n\n\n"
},

{
    "location": "api/properties/#Molecular-properties-1",
    "page": "Molecular properties",
    "title": "Molecular properties",
    "category": "section",
    "text": "molweight(mol::VectorMol; digits=2)\nH_acceptor_count(mol::VectorMol)\nH_donor_count(mol::VectorMol)\nwclogp(mol::VectorMol; digits=2)\nrotatable_count(mol::VectorMol)"
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
