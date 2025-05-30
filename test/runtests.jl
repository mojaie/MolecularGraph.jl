#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraphTest

using Graphs
using Logging
using MolecularGraph
using MolecularGraph:
    combinations, sortstablemax, sortstablemin,
    logfactorial,
    MGPoint, Coordinates, radiantophase,
    Point2D, Point3D, Segment,
    toarray, x_components, y_components, z_components,
    distance, unitvector, midpoint,
    translate, trim_u, trim_v, trim_uv,
    cross2d, interiorangle, isclockwise,
    transformmatrix, rotation,
    Color, atom_color, is_atom_visible,
    sdf_bond_style, double_bond_style, atomhtml,
    merge_ds!,
    Metadata,
    ordered_neighbors, u_edge, edge_rank,
    atomsymbol!, atomprop!, atom!,
    bondsymbol!, bond!,
    querypropmap, generate_queryfunc, querymatch, optimize_query,
    fragment!,specialize_nonaromatic!, remove_hydrogens!,
    lookahead, forward!, backtrack!,
    lgnot!, lghighand!, lgor!, lglowand!,
    chain!, fragment!, componentquery!,
    resolve_disjoint_not, resolve_recursive, generate_truthtable, querymatch, querypropmap,
    ctab_atom_v2, ctab_bond_v2,
    exact_match_prefilter, substruct_match_prefilter
using OrderedCollections: OrderedDict
using Test


include("./util/iterator.jl")
include("./util/math.jl")

include("./geometry/cartesian.jl")
# include("./geometry/internal.jl")

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

include("sdfilereader.jl")
include("sdfilewriter.jl")
include("json.jl")
include("./smarts/base.jl")
include("./smarts/logicaloperator.jl")
include("./smarts/atom.jl")
include("./smarts/bond.jl")
include("./smarts/molecule.jl")
include("./smarts/smiles.jl")
include("./smarts/smarts.jl")

include("stereo.jl")
include("preprocess.jl")
include("coords.jl")
include("properties.jl")
include("mass.jl")
include("wclogp.jl")
include("inchi.jl")
include("structurematch.jl")
include("structurematch_mcs.jl")
include("./model/query.jl")
include("querycontainment.jl")
include("./draw/base.jl")
include("./draw/svg.jl")
include("./draw/3d.jl")

end
