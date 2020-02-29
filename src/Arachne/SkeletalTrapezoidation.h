//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_SKELETAL_TRAPEZOIDATION_H
#define ARACHNE_SKELETAL_TRAPEZOIDATION_H

#include <boost/polygon/voronoi.hpp>

#include <unordered_map>
#include <utility> // pair

#include "utils/HalfEdgeGraph.h"
#include "utils/polygon.h"
#include "utils/PolygonsSegmentIndex.h"
#include "SkeletalTrapezoidationEdge.h"
#include "SkeletalTrapezoidationNode.h"

#include "utils/STLwriter.h"

namespace arachne
{

extern bool generate_MAT_STL;

/*!
 * The input polygon region is decomposed into trapezoids and represented as a half-edge data-structure.
 * 
 * We compute the voronoi diagram on the segments of the input polygon using boost.
 * 
 * Then we copy the internal edges from the VD into the half-edge data-structure,
 * and discretize parabolic edges.
 * 
 * Then we clean up the graph to avoid rounding erros because of integer arithmetic.
 * 
 */
class SkeletalTrapezoidation
{
protected:
    using pos_t = double;
    using vd_t = boost::polygon::voronoi_diagram<pos_t>;
    using graph_t = HalfEdgeGraph<SkeletalTrapezoidationEdge, SkeletalTrapezoidationNode>;
    using edge_t = SkeletalTrapezoidationEdge;
    using node_t = SkeletalTrapezoidationNode;

    const Polygons& polys; //!< input outline boundary shape

    float transitioning_angle; //!< How pointy a region should be before we apply the method. Equals 180* - limit_bisector_angle
    coord_t discretization_step_size; //!< approximate size of segments when parabolic VD edges get discretized (and vertex-vertex edges)
    coord_t snap_dist = 20; //!< generic arithmatic inaccuracy. Only used to determine whether a transition really needs to insert an extra edge.

public:
    using Segment = PolygonsSegmentIndex;

    SkeletalTrapezoidation(const Polygons& polys, float transitioning_angle, coord_t discretization_step_size)
    : polys(polys)
    , transitioning_angle(transitioning_angle)
    , discretization_step_size(discretization_step_size)
    {
    }

    graph_t graph;

    /*!
     * Compute the skeletal trapezoidation decomposition of the input shape.
     * 
     * Compute the Voronoi Diagram (VD) and transfer all inside edges into our half-edge (HE) datastructure.
     * 
     * The algorithm is currently a bit overcomplicated, because the discretization of parabolic edges is performed at the same time as all edges are being transfered,
     * which means that there is no one-to-one mapping from VD edges to HE edges.
     * Instead we map from a VD edge to the last HE edge.
     * This could be cimplified by recording the edges which should be discretized and discretizing the mafterwards.
     * 
     * Another complication arises because the VD uses floating logic, which can result in zero-length segments after rounding to integers.
     * We therefore collapse edges and their whole cells afterwards.
     */
    void initializeGraph();

private:

    /*!
     * mapping each voronoi VD edge to the corresponding halfedge HE edge
     * In case the result segment is discretized, we map the VD edge to the *last* HE edge
     */
    std::unordered_map<vd_t::edge_type*, edge_t*> vd_edge_to_he_edge;
    std::unordered_map<vd_t::vertex_type*, node_t*> vd_node_to_he_node;
    node_t& makeNode(vd_t::vertex_type& vd_node, Point p); //!< Get the node which the VD node maps to, or create a new mapping if there wasn't any yet.
    /*!
     * Transfer an edge vrom the VD to the HE and perform discretization of parabolic edges (and vertex-vertex edges)
     * \p prev_edge serves as input and output. May be null as input.
     */
    void transferEdge(Point from, Point to, vd_t::edge_type& vd_edge, edge_t*& prev_edge, Point& start_source_point, Point& end_source_point, const std::vector<Point>& points, const std::vector<Segment>& segments);

    /*!
     * Make a support edge (a.k.a. rib) from a node to the closest location on the polygon.
     * 
     *    introduced
     * B  rib   .-
     * |  ^    /     }
     * |_ _ _./      }-> parabolic edge which gets discretized,
     * |     |       }   so that a new node is introduced,
     * |     |       }   which is not in the VD.
     * A-----o
     * 
     * \param prev_edge The node, as pointed to by an already created edge.
     * \param start_source_point Start point of polygon segment (A in drawing)
     * \param end_source_point End point of polygon segment (B in drawing)
     */
    void makeRib(edge_t*& prev_edge, Point start_source_point, Point end_source_point);

    /*!
     * Discretize a nonlinear edge into linear segments.
     * 
     * A non-linear edge is defined as an edge where the radial distance increases non-linearly.
     * 1. A parabolic edge, the midpoints between an input polygon segment and an input polygon vertex
     * 2. A vertex-vertex edge, the midpoints between two input polygon vertices
     * 
     * \param segment The nonlinear segment to discretize
     * \param points The input points (always empty)
     * \param segments The input segments of the input polygon
     * \return list of points on the nonlinear edge
     */
    std::vector<Point> discretize(const vd_t::edge_type& segment, const std::vector<Point>& points, const std::vector<Segment>& segments);

    /*!
     * Compute the range of VD edges belonging inside the polygon (if any).
     * 
     * \param cell The VD cell for which to compute the internal range of edges
     * \param[out] start_source_point output the start point of the source segment of this cell
     * \param[out] end_source_point output the end point of the source segment of this cell
     * \param[out] starting_vd_edge output the first edge, departing from the polyon going inward
     * \param[out] ending_vd_edge output the last edge, departing from inside the polyon and going to the polygon
     * \return whether this cell has any associated edges on the inside
     */
    bool computePointCellRange(vd_t::cell_type& cell, Point& start_source_point, Point& end_source_point, vd_t::edge_type*& starting_vd_edge, vd_t::edge_type*& ending_vd_edge, const std::vector<Point>& points, const std::vector<Segment>& segments);

    /*!
     * Compute the range of VD edges belonging inside the polygon (if any).
     * 
     * \param cell The VD cell for which to compute the internal range of edges
     * \param[out] start_source_point output the start point of the source segment of this cell
     * \param[out] end_source_point output the end point of the source segment of this cell
     * \param[out] starting_vd_edge output the first edge, departing from the polyon going inward
     * \param[out] ending_vd_edge output the last edge, departing from inside the polyon and going to the polygon
     * \return whether this cell has any associated edges on the inside
     */
    void computeSegmentCellRange(vd_t::cell_type& cell, Point& start_source_point, Point& end_source_point, vd_t::edge_type*& starting_vd_edge, vd_t::edge_type*& ending_vd_edge, const std::vector<Point>& points, const std::vector<Segment>& segments);

    /*!
     * For VD cells associated with an input polygon vertex, we need to separate the node at the end and start of the cell into two
     * That way we can reach both the quad_start and the quad_end from the [some_edge] of the two new nodes
     * Otherwise if node.some_edge = quad_start you couldnt reach quad_end.twin by normal iteration (i.e. it = it.twin.next)
     */
    void separatePointyQuadEndNodes();

    /*!
     * If an edge is too small, collapse it and its twin and fix the surrounding edges to ensure a consistent graph.
     * 
     * Don't collapse support edges, unless we can collapse the whole quad.
     * 
     * o-,
     * |  "-o
     * |    | > Don't collapse this edge only.
     * o    o
     */
    void collapseSmallEdges(coord_t snap_dist = 5);

    /*!
     * Fix an issue where two nodes are in the graph at the same location, while there should only be a single one.
     */
    void fixNodeDuplication();

public:
    // debugging
    void debugCheckGraphCompleteness(); //!< Checks whether all member fields of edges and nodes are filled. Should be true after initialization.
    void debugCheckEndpointUniqueness(); //!< Checks whether the end points of qauds have unique verts. Should be true after separatePointyQuadEndNodes().
    void debugCheckGraphExistance(); //!< Checks whether all member fields of edges and nodes are existing nodes/edges recorded in graph.nodes and graph.edges. Should be true after any graph update.
    void debugCheckGraphStructure(); //!< Checks whether iterating around a node (using it = it.twin.next) ends up where it started. Should be true after init.
    void debugCheckGraphReachability(); //!< Checks whether an edge is reachable from iterating around its from node. Should be true after init.
    void debugCheckGraphConsistency(bool ignore_duplication = false); //!< Checks whether edge and node relations fit with each other. Should be true after any graph update.
    void debugOutput(SVG& svg, bool draw_arrows, bool draw_dists, bool draw_bead_counts = false, bool draw_locations = false);
    void debugOutputSTL(STLwriter& stl, bool use_bead_count = false);
protected:
    SVG::ColorObject getColor(edge_t& edge);

};




} // namespace arachne
#endif // ARACHNE_SKELETAL_TRAPEZOIDATION_H
