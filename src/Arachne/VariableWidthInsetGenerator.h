//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_VARIABLE_WIDTH_INSET_GENERATOR_H
#define ARACHNE_VARIABLE_WIDTH_INSET_GENERATOR_H

#include <boost/polygon/voronoi.hpp>

#include <unordered_map>
#include <utility> // pair

#include "utils/HalfEdgeGraph.h"
#include "utils/polygon.h"
#include "utils/PolygonsSegmentIndex.h"
#include "utils/ExtrusionJunction.h"
#include "utils/ExtrusionLine.h"
#include "SkeletalTrapezoidation.h"
#include "SkeletalTrapezoidationEdge.h"
#include "SkeletalTrapezoidationJoint.h"
#include "BeadingStrategies/BeadingStrategy.h"

#include "utils/STLwriter.h"

namespace arachne
{

/*!
 * Main class of libArachne.
 * 
 * The input polygon region is decomposed into trapezoids and represented as a half-edge data-structure by \ref SkeletalTrapezoidation.
 * 
 * We determine which edges are 'central' accordinding to the transitioning_angle of the beading strategy,
 * and determine the bead count for these central regions and apply them outward when generating toolpaths. [oversimplified]
 * 
 * The method can be visually explained as generating the 3D union of cones surface on the outline polygons,
 * and changing the heights along central regions of that surface so that they are flat. (See paper)
 * This visual explanation aid explains the use of "upward", "lower" etc,
 * i.e. the radial distance and/or the bead count are used as heights of this visualization, there is no coordinate called 'Z'.
 */
class VariableWidthInsetGenerator
{
    using pos_t = double;
    using vd_t = boost::polygon::voronoi_diagram<pos_t>;
    using graph_t = HalfEdgeGraph<SkeletalTrapezoidationJoint, SkeletalTrapezoidationEdge>;
    using edge_t = HalfEdge<SkeletalTrapezoidationJoint, SkeletalTrapezoidationEdge>;
    using node_t = HalfEdgeNode<SkeletalTrapezoidationJoint, SkeletalTrapezoidationEdge>;
    using Beading = BeadingStrategy::Beading;

public:
    SkeletalTrapezoidation st;

private:
    const Polygons& polys; //!< input outline boundary shape

    float transitioning_angle; //!< How pointy a region should be before we apply the method. Equals 180* - limit_bisector_angle
    coord_t discretization_step_size; //!< approximate size of segments when parabolic VD edges get discretized (and vertex-vertex edges)
    coord_t transition_filter_dist; //!< filter transition mids (i.e. anchors) closer together than this
    coord_t beading_propagation_transition_dist; //!< When there are different beadings propagated from below and from above, use this transitioning distance
    coord_t marking_filter_dist = 20; //!< filter areas marked as 'central' smaller than this
    coord_t snap_dist = 20; //!< generic arithmatic inaccuracy. Only used to determine whether a transition really needs to insert an extra edge.

public:
    VariableWidthInsetGenerator(const Polygons& polys, float transitioning_angle
    , coord_t discretization_step_size = 200
    , coord_t transition_filter_dist = 1000
    , coord_t beading_propagation_transition_dist = 400
    );

    /*!
     * Generate the variable width insets.
     * 
     * \param beading_strategy The strategy to fill space using (variable width0 insets
     * \param filter_outermost_marked_edges Unmark all edges touching the input polygon as non-central. (Used to emulate existing literature)
     */
    std::vector<std::list<ExtrusionLine>> generateToolpaths(const BeadingStrategy& beading_strategy, bool filter_outermost_marked_edges = false);

private:

    /*!
     * Generate ExtrusionLines.
     * 
     * \param[out] segments the generated segments
     */
    void generateSegments(std::vector<std::list<ExtrusionLine>>& result_polylines_per_index, const BeadingStrategy& beading_strategy);

    /*!
     * Get the edge pointing to the node with the maximum distance_to_boundary
     * from the quad/trapezoid for which  \p quad_start_edge is the first edge.
     * 
     * max R
     * ^
     * o.,    } return edge
     * |  "-o
     * |    | } quad_start_edge
     * o====o
     *  poly
     *  segment
     */
    edge_t* getQuadMaxRedgeTo(edge_t* quad_start_edge);

    /*!
     * Helper object to store Beadings on unmarked nodes.
     * These beadings originate from marked nodes, from where they are propagated to unmakred nodes.
     */
    struct BeadingPropagation
    {
        Beading beading;
        coord_t dist_to_bottom_source;
        coord_t dist_from_top_source;
        bool is_upward_propagated_only; //!< Whether the doward propagation phase has not processed this beading yet.
        BeadingPropagation(const Beading& beading)
        : beading(beading)
        , dist_to_bottom_source(0)
        , dist_from_top_source(0)
        , is_upward_propagated_only(false)
        {}
    };

    /*!
     * propagate beading info from lower R nodes to higher R nodes
     * 
     * only propagate from nodes with beading info upward to nodes without beading info
     * 
     * edges are sorted so that we can do a depth-first walk without employing a recursive algorithm
     * 
     * In upward propagated beadings we store the distance traveled, so that we can merge these beadings with the downward propagated beadings in \ref propagateBeadingsDownward(.)
     * 
     * \param upward_quad_mids all upward halfedges of the inner skeletal edges (not directly connected to the outline) sorted on their highest [distance_to_boundary]. Higher dist first.
     * \param node_to_beading mapping each node to a \ref BeadingPropagation
     * \param beading_strategy The beading strategy
     */
    void propagateBeadingsUpward(std::vector<edge_t*>& upward_quad_mids, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, const BeadingStrategy& beading_strategy);

    /*!
     * propagate beading info from higher R nodes to lower R nodes
     * 
     * merge with upward propagated beadings if they are encountered
     * 
     * don't transfer to nodes which lie on the outline polygon
     * 
     * edges are sorted so that we can do a depth-first walk without employing a recursive algorithm
     * 
     * \param upward_quad_mids all upward halfedges of the inner skeletal edges (not directly connected to the outline) sorted on their highest [distance_to_boundary]. Higher dist first.
     * \param node_to_beading mapping each node to a \ref BeadingPropagation
     * \param beading_strategy The beading strategy
     */
    void propagateBeadingsDownward(std::vector<edge_t*>& upward_quad_mids, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, const BeadingStrategy& beading_strategy);

    /*!
     * propagate beading info from higher R nodes to lower R nodes
     * 
     * merge with upward propagated beadings if they are encountered
     * 
     * don't transfer to nodes which lie on the outline polygon
     * 
     * \param upward_quad_mids all upward halfedges of the inner skeletal edges (not directly connected to the outline) sorted on their highest [distance_to_boundary]. Higher dist first.
     * \param node_to_beading mapping each node to a \ref BeadingPropagation
     * \param beading_strategy The beading strategy
     */
    void propagateBeadingsDownward(edge_t* edge_to_peak, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, const BeadingStrategy& beading_strategy);

    /*!
     * Interpolate between the beading propagated from above and the one propagated from below.
     * 
     * \param top The beading propagated from above
     * \param ratio_top_to_whole The ratio of the top beading to appear in the result
     * \param bottom The beading propagated from below
     * \param switching_radius The radius at which we switch from the top beading to a merged one
     * \return the interpolated/merged beading
     */
    Beading interpolate(const Beading& top, float ratio_top_to_whole, const Beading& bottom, coord_t switching_radius) const;

    /*!
     * Interpolate between the beading propagated from above and the one propagated from below.
     * 
     * \param top The beading propagated from above
     * \param ratio_top_to_whole The ratio of the top beading to appear in the result
     * \param bottom The beading propagated from below
     * \return the interpolated/merged beading
     */
    Beading interpolate(const Beading& top, float ratio_top_to_whole, const Beading& bottom) const;

    /*!
     * node_to_beading[node]
     * with extra safety instructions in case that function call should fail somehow.
     * \param node_to_beading mapping each node to a \ref BeadingPropagation
     * \param beading_strategy The beading strategy
     */
    BeadingPropagation& getBeading(node_t* node, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, const BeadingStrategy& beading_strategy);

    /*!
     * In case we cannot find the beading of a node, get a beading from the nearest node
     * 
     * \param node The node for which we weren't able to find a beading
     * \param max_dist The maximum distance from the \p node after which to give up
     * \param node_to_beading mapping each node to a \ref BeadingPropagation
     */
    BeadingPropagation* getNearestBeading(node_t* node, coord_t max_dist, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading);

    /*!
     * generate junctions for each unmarked edge
     * \param node_to_beading mapping each node to a \ref BeadingPropagation
     * \param[out] edge_to_junctions junctions ordered high R to low R
     * \param beading_strategy The beading strategy
     */
    void generateJunctions(std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions, const BeadingStrategy& beading_strategy);

    /*!
     * connect junctions in each quad
     * \param edge_to_junctions junctions ordered high R to low R
     * \param[out] result_polylines_per_index the generated segments
     */
    void connectJunctions(std::unordered_map<edge_t*, std::vector<arachne::ExtrusionJunction>>& edge_to_junctions, std::vector<std::list<ExtrusionLine>>& result_polylines_per_index);

    /*!
     * Whether more than 2 extrusion lines end in this \p node.
     * I.e. whether this node is the end point of a transition and it has an odd bead count.
     */
    bool isMultiIntersection(node_t* node);

    /*!
     * Genrate small extrusion segments for local maxima where the beading would only result in a single bead.
     * \param[out] segments the generated segments
     */
    void generateLocalMaximaSingleBeads(std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, std::vector<std::list<ExtrusionLine>>& result_polylines_per_index);

    /*!
     * edge_to_junctions[edge]
     * 
     * \p edge is assumed to point upward to higher R; otherwise take its twin
     * 
     * \param include_odd_start_junction Whether to leave out the first junction if it coincides with \p edge.from->p
     */
    const std::vector<ExtrusionJunction>& getJunctions(edge_t* edge, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions);
    
    // ^ toolpath generation | v helpers

    // TODO: remove duplicate helper functions!!

    /*!
     * Return the first and last edge of the edges replacing \p edge pointing to the same node
     */
    std::pair<edge_t*, edge_t*> insertRib(edge_t& edge, node_t* mid_node);

    std::pair<Point, Point> getSource(const edge_t& edge);
    bool isEndOfMarking(const edge_t& edge) const;

    /*!
     * Check whether this node has a locally maximal distance_to_boundary
     * 
     * \param strict Whether equidistant edges can count as a local maximum
     */
    bool isLocalMaximum(const node_t& node, bool strict = false) const;

    /*!
     * Check (recursively) whether there is any upward edge from the distance_to_boundary of the from of the \param edge
     * 
     * \param strict Whether equidistant edges can count as a local maximum
     */
    bool canGoUp(const edge_t* edge, bool strict = false) const;

    /*!
     * Calculate the traversed distance until we meet an upward edge.
     * Useful for calling on edges between equidistant points.
     * 
     * If we can go up then the distance includes the length of the \param edge
     */
    std::optional<coord_t> distToGoUp(const edge_t* edge) const;

    /*!
     * Check whether the edge goes from a lower to a higher distance_to_boundary.
     * Effectively deals with equidistant edges by looking beyond this edge.
     * 
     * // TODO: remove this function
     */
    bool isUpward(const edge_t* edge) const;

    /*!
     * Whether the node is connected to any marked edge
     */
    bool isMarked(const node_t* node) const;
public:
    void debugCheckDecorationConsistency(bool transitioned); //!< Check logical relationships relting to distance_to_boundary and is_marked etc. Should be true anywhere after setMarking(.)
    void debugOutput(SVG& svg, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions);
    void debugOutput(STLwriter& stl, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading);
protected:
    SVG::ColorObject getColor(edge_t& edge);

};




} // namespace arachne
#endif // ARACHNE_VARIABLE_WIDTH_INSET_GENERATOR_H
