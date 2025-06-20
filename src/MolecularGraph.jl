#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraph

using OrderedCollections
using Printf: @sprintf
import YAML


# Utility

using GeometryBasics:
    mesh, Cylinder, Sphere, Point, Point2d, Point3d, Line
using LinearAlgebra:
    cross, dot, norm, normalize
import LinearAlgebra

include("./util/geometry.jl")
include("./util/meta.jl")
include("./util/iterator.jl")
include("./util/math.jl")


# Interfaces

using Graphs

include("./model/interface.jl")
include("./property/interface.jl")
include("./draw/interface.jl")

export
    Reaction, QueryTree,
    vproptype, eproptype,
    props, get_prop, has_prop, set_prop!,
    u_edge, ordered_neighbors, edge_neighbors, ordered_edge_neighbors


# Graph models and algorithms

using DelimitedFiles: readdlm

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
    MolGraph, SDFMolGraph, SMILESMolGraph, CommonChemMolGraph,
    QueryMolGraph, QueryAtom, QueryBond

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
    induced_subgraph_edges, modular_product, line_graph,
    planaritytest, outerplanaritytest,
    isplanar, isouterplanar


# Basic molecular properties

include("property/topology.jl")
include("property/valence.jl")
include("property/hybridization.jl")
include("property/wclogp.jl")

export
    sssr, sssr!,
    which_ring, edge_which_ring, fused_rings, which_fused_ring,
    smallest_ring, ring_count, is_in_ring, is_edge_in_ring,
    default_atom_charge!, default_bond_order!,
    lone_pair, lone_pair!, apparent_valence, apparent_valence!, valence, valence!,
    explicit_hydrogens, implicit_hydrogens, heavy_atoms,
    total_hydrogens, connectivity,
    is_hydrogen_donor, hydrogen_donor_count,
    is_hydrogen_acceptor, hydrogen_acceptor_count,
    is_rotatable, rotatable_count,
    atom_counter, heavy_atom_count, molecular_formula, empirical_formula,
    pi_electron, pi_delocalized, hybridization, hybridization_delocalized,
    is_ring_aromatic, is_ring_aromatic!, is_aromatic, is_edge_aromatic,
    wclogp


# Preprocessing and molecule manipulation

using coordgenlibs_jll: libcoordgen

include("coords.jl")
include("stereo.jl")
include("preprocess.jl")

export
    coords2d, has_coords2d, coords3d, has_coords3d,
    coordgen, coordgen!, coords_from_sdf!, update_coords!,
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
    sdf_on_init!, sdf_on_update!,
    sdfilereader, rdfilereader, sdfilescanner,
    printv2mol, printv2sdf, sdfilewriter,
    sdftomol, rxntoreaction,
    SMILESParser, SMARTSParser, SMARTSMolGraph,
    smiles_on_init!, smiles_on_update!,
    smilestomol, smartstomol


# Descriptors

using libinchi_jll: libinchi

include("mass.jl")
include("rdkit.jl")
include("inchi.jl")
include("structurematch.jl")
include("querycontainment.jl")

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

export query_containment_diagram, find_queries


# Molecule drawing

import Cairo
using Colors: RGB, RGBA, N0f8, hex, coloralpha
using MakieCore:
    @recipe, Theme, meshscatter!, mesh!
import MakieCore
import Statistics


include("./draw/color.jl")
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


# Package compiler

using Base: unsafe_convert
using Base64: Base64EncodePipe

include("libmoleculargraph.jl")

export
    vertex_count, edge_count,
    molblock, sdfmolblock,
    tdmcis_size, tdmces_size, tdmcis_gls, tdmces_gls

end
