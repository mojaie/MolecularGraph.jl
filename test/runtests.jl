#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraphTest

using Colors: RGB
using GeometryBasics: Point, Point2d, Line
using Graphs
using LinearAlgebra: cross
using Logging
using MolecularGraph

using MolecularGraph:
    combinations, sortstablemax, sortstablemin,
    logfactorial,
    translate, trim_u, trim_v, trim_uv,
    interiorangle, isclockwise, transformmatrix,
    merge_ds!

using MolecularGraph:
    ordered_neighbors, u_edge, edge_rank,
    MolProperty, MolDescriptor, remap!, reconstruct, sdf_on_init!,
    atomsymbol!, atomprop!, atom!,
    bondsymbol!, bond!,
    QueryNode, root, querytree, canonical,
    add_qnode!, add_qedge!, set_qnode!, add_qedge!,
    rem_qnode!, rem_qnodes!, rem_qedge!,
    qeq, qtrue, qand, qor, qnot, qanytrue,
    querypropmap, generate_queryfunc

using MolecularGraph:
    wclogptype, wclogphydrogentype

using MolecularGraph:
    fragment!, specialize_nonaromatic!, resolve_not_hydrogen!, remove_hydrogens!,
    lookahead, forward!, backtrack!,
    lgnot!, lgoperator!, lghighand!, lgor!, lglowand!,
    chain!, fragment!, componentquery!,
    ctab_atom_v2, ctab_bond_v2

using MolecularGraph:
    exact_topology_prefilter, topology_prefilter,
    QueryTruthTable, querymatch,
    resolve_disjoint_not!, resolve_recursive!, generate_truthtable,
    draw2d_bond_style,
    atom_color, is_atom_visible,
    double_bond_style, atomhtml

using MolecularGraph: smiles

using OrderedCollections: OrderedDict
using Test

include("./util/geometry.jl")
include("./util/iterator.jl")
include("./util/math.jl")

include("./graph/operators.jl")
include("./graph/cycle.jl")
include("./graph/bipartite.jl")
include("./graph/matching.jl")
include("./graph/planarity.jl")
include("./graph/clique.jl")
include("./graph/isomorphism_vf2.jl")
include("./graph/isomorphism_clique.jl")

include("./model/atom.jl")
include("./model/bond.jl")
include("./model/molgraph.jl")
include("./model/query.jl")

include("./smarts/base.jl")
include("./smarts/logicaloperator.jl")
include("./smarts/atom.jl")
include("./smarts/bond.jl")
include("./smarts/molecule.jl")
include("./smarts/smiles.jl")
include("./smarts/smarts.jl")
include("sdfilereader.jl")
include("sdfilewriter.jl")

include("properties.jl")
include("wclogp.jl")
include("stereo.jl")
include("preprocess.jl")
include("coords.jl")
include("mass.jl")
include("inchi.jl")

include("structurematch.jl")
include("structurematch_mcs.jl")
include("querycontainment.jl")

include("json.jl")
include("./draw/base.jl")
include("./draw/svg.jl")
include("./draw/3d.jl")

end
