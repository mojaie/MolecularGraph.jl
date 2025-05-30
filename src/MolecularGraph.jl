#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraph

using OrderedCollections
using Printf: @sprintf
import YAML


# Utility

include("./util/meta.jl")
include("./util/iterator.jl")
include("./util/math.jl")


# Geometry
# TODO: migrate to GeometryBasics

using LinearAlgebra: cross, dot, norm, normalize
import LinearAlgebra

include("./geometry/interface.jl")
include("./geometry/cartesian.jl")
# include("./geometry/internal.jl")


# Graph models and algorithms

using DelimitedFiles: readdlm
using Graphs

include("./model/interface.jl")
include("./model/atom.jl")
include("./model/bond.jl")
include("./model/molgraph.jl")
include("./model/query.jl")

include("./graph/operators.jl")
include("./graph/traversals.jl")
include("./graph/cycle.jl")
include("./graph/bipartite.jl")
include("./graph/matching.jl")
include("./graph/clique.jl")
include("./graph/planarity.jl")
include("./graph/isomorphism_edge.jl")
include("./graph/isomorphism_vf2.jl")
include("./graph/isomorphism_clique.jl")

export
    ATOMTABLE, ATOMSYMBOLMAP, ATOM_COVALENT_RADII, ATOM_VANDERWAALS_RADII,
    SDFAtom, SMILESAtom, CommonChemAtom,
    atom_number, atom_symbol, atom_charge, multiplicity, atom_mass,
    SDFBond, SMILESBond, CommonChemBond,
    bond_order,
    AbstractMolGraph, SimpleMolGraph,
    MolGraph, SDFMolGraph, SMILESMolGraph, CommonChemMolGraph,
    SMARTSMolGraph,
    QueryAny, QueryLiteral, QueryOperator, QueryTree, QueryTruthTable,
    AbstractReaction, Reaction,
    vproptype, eproptype,
    props, vprops, eprops,
    get_prop, has_prop,
    set_state!, has_state, get_state,
    set_cache!, has_cache, get_cache, clear_caches!,
    set_prop!, update_edge_rank!,
    edge_neighbors, ordered_edge_neighbors

export
    maxcardmap, maxcard,
    all_maximal_cliques, maximum_clique,
    all_maximal_conn_cliques, maximum_conn_clique,
    approx_maximum_clique,
    mincyclebasis, edgemincyclebasis, disjoint_union,
    ConstraintArrayMCIS, ConstraintArrayMCES,
    maximum_common_subgraph, maximum_common_edge_subgraph,
    mcis_constraints, mces_constraints,
    VF2Matcher,
    isomorphisms, is_isomorphic,
    nodesubgraph_isomorphisms, nodesubgraph_is_isomorphic,
    edgesubgraph_isomorphisms, edgesubgraph_is_isomorphic,
    subgraph_monomorphisms, subgraph_is_monomorphic,
    max_matching, is_perfect_matching,
    modular_product, line_graph,
    planaritytest, outerplanaritytest,
    isplanar, isouterplanar


# Preprocessing

using coordgenlibs_jll: libcoordgen

include("coords.jl")
include("stereo.jl")
include("preprocess.jl")

export
    has_coords, sdf_coords2d, coords2d, coords3d, coordgen, coordgen!,
    set_stereocenter!, set_stereobond!, stereo_hydrogen, safe_stereo_hydrogen!,
    stereocenter_from_smiles!, stereocenter_from_sdf2d!,
    stereobond_from_smiles!, stereobond_from_sdf2d!,
    kekulize, kekulize!,
    removable_hydrogens, all_hydrogens,
    remove_hydrogens!, remove_all_hydrogens!, add_hydrogens!,
    largest_component_nodes, extract_largest_component!,
    protonate_acids, protonate_acids!, deprotonate_oniums, deprotonate_oniums!,
    depolarize, depolarize!, polarize, polarize!,
    to_triple_bond, to_triple_bond!, to_allene_like, to_allene_like!


# Molecule drawing

import Cairo
using Colors: RGB, RGBA, N0f8, hex
using GeometryBasics: mesh, Cylinder, Sphere, Point
using MakieCore: @recipe, Theme, meshscatter!, lines!, mesh!
import MakieCore
import Statistics


include("./draw/color.jl")
include("./draw/interface.jl")
include("./draw/draw2d.jl")
include("./draw/svg.jl")
include("./draw/cairo.jl")
include("./draw/draw3d.jl")

export
    drawsvg, drawpng,
    html_fixed_size, html_grid,
    atom_radius,
    spacefilling, spacefilling!,
    ballstick, ballstick!,
    stick, stick!,
    wire, wire!


# I/O

import Dates
import JSON

include("json.jl")
include("sdfilereader.jl")
include("sdfilewriter.jl")
include("./smarts/base.jl")
include("./smarts/atom.jl")
include("./smarts/bond.jl")
include("./smarts/logicaloperator.jl")
include("./smarts/molecule.jl")

export
    to_dict, to_json,
    SDFileReader,
    sdfilereader, rdfilereader, sdfilescanner,
    printv2mol, printv2sdf, sdfilewriter,
    sdftomol, rxntoreaction,
    SMILESParser, SMARTSParser,
    smilestomol, smartstomol


# Descriptors

using libinchi_jll: libinchi
using RDKitMinimalLib:
    Mol, get_mol, get_smiles,
    get_morgan_fp_as_bytes, get_rdkit_fp_as_bytes,
    get_pattern_fp_as_bytes, get_atom_pair_fp_as_bytes,
    get_topological_torsion_fp_as_bytes

include("properties.jl")
include("mass.jl")
include("wclogp.jl")
include("rdkit.jl")
include("inchi.jl")
include("structurematch.jl")
include("querycontainment.jl")

export
    sssr, sssr!,
    which_ring, edge_which_ring, fused_rings, which_fused_ring,
    smallest_ring, ring_count, is_in_ring, is_edge_in_ring,
    lone_pair, lone_pair!, apparent_valence, apparent_valence!, valence, valence!,
    explicit_hydrogens, implicit_hydrogens, heavy_atoms,
    total_hydrogens, connectivity,
    is_hydrogen_donor, hydrogen_donor_count,
    is_hydrogen_acceptor, hydrogen_acceptor_count,
    is_rotatable, rotatable_count,
    atom_counter, heavy_atom_count, molecular_formula, empirical_formula,
    pi_electron, pi_delocalized, hybridization, hybridization_delocalized,
    is_ring_aromatic, is_ring_aromatic!, is_aromatic, is_edge_aromatic

export
    monoiso_mass_unc, monoiso_mass, nominal_mass,
    exact_mass_unc, exact_mass,
    standard_weight_unc, standard_weight,
    isotopic_composition, massspec_peaks, simulate_massspec

export
    to_rdkdict, to_rdkjson,
    rdkitmol, smiles,
    morgan_fp_vector, rdkit_fp_vector,
    pattern_fp_vector, atom_pair_fp_vector,
    topological_torsion_fp_vector

export
    inchi, inchikey, inchitomol, inchitosdf

export
    vmatchgen, vmatchvecgen, ematchgen, ematchvecgen,
    exact_matches, has_exact_match,
    substruct_matches, has_substruct_match,
    node_substruct_matches, has_node_substruct_match,
    edge_substruct_matches, has_edge_substruct_match,
    disconnected_mcis, disconnected_mces,
    connected_mcis, connected_mces,
    tcmcis, tcmces, tdmcis, tdmces,
    emaptonmap,
    tdmcis_constraints, tdmces_constraints

export
    wclogptype, wclogphydrogentype,
    wclogpcontrib, wclogp

export query_containment_diagram, find_queries


# Package compiler

using Base: unsafe_convert
using Base64: Base64EncodePipe

include("libmoleculargraph.jl")

export
    vertex_count, edge_count,
    molblock, sdfmolblock,
    tdmcis_size, tdmces_size, tdmcis_gls, tdmces_gls

end
