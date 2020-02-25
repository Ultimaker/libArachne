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
class VariableWidthInsetGenerator : public SkeletalTrapezoidation
{
    using Beading = BeadingStrategy::Beading;

    const Polygons& polys; //!< input outline boundary shape

    coord_t transition_filter_dist; //!< filter transition mids (i.e. anchors) closer together than this
    coord_t beading_propagation_transition_dist; //!< When there are different beadings propagated from below and from above, use this transitioning distance
    coord_t marking_filter_dist = 20; //!< filter areas marked as 'central' smaller than this

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

protected:
    
    void setMarking(const BeadingStrategy& beading_strategy); //!< set the is_marked flag for each edge based on the transitioning_angle

    /*!
     * Filter out small marked areas.
     * 
     * Only used to get rid of small edges which get marked because of rounding errors because hte region is so small.
     */
    void filterMarking(coord_t max_length);

    /*!
     * Filter markings connected to starting_edge recursively.
     * Recursive function.
     * Unmark a region as central if it's too small.
     * 
     * \param starting_edge Edge from which to walk along the marked edges to check if the region is too small
     * \param traveled_dist Total travelled distance before the current recursive call
     * \param max_length Maximum total length of region which is still filtered out.
     * \return Whether we should unmark this marked section on the way back out of the recursion
     */
    bool filterMarking(edge_t* starting_edge, coord_t traveled_dist, coord_t max_length);

    /*!
     * Unmark the outermost edges directly connected to the outline.
     * 
     * Only used to emulate some related literature.
     * 
     * The paper shows that this function is bad for the stability of the framework.
     */
    void filterOuterMarking();

    /*!
     * Set bead count in marked regions based on the optimal_bead_count of the beading strategy.
     */
    void setBeadCount(const BeadingStrategy& beading_strategy);

    /*!
     * Add marked regions and set bead counts
     * where there is an end of marking and when traveling upward we get to another region with the same bead count
     */
    void filterUnmarkedRegions(const BeadingStrategy& beading_strategy);

    /*!
     * Add marked regions and set bead counts
     * where there is an end of marking and when traveling upward we get to another region with the same bead count
     * 
     * Recursive function.
     * 
     * \param to_edge edge pointing to the node from which to continue the search
     * \param bead_count Bead count at the lower end from where we started ascending
     * \param traveled_dist Distance traversed up to the current recursive call
     * \param max_dist Maximum distance along the graph to change from unmarked to marked
     * \param beading_strategy Strategy to determine bead count on nodes along the region in case we decide to mark it
     * \return Whether to set the bead count on the way back
     */
    bool filterUnmarkedRegions(edge_t* to_edge, coord_t bead_count, coord_t traveled_dist, coord_t max_dist, const BeadingStrategy& beading_strategy);

//
// ^^^^^^^^^^^^^^^^^^^^^
//       MARKING
// =====================
//
// =====================
//    TRANSTISIONING
// vvvvvvvvvvvvvvvvvvvvv
//

    
    /*!
     * Representing the location along an edge where the anchor position of a transition should be placed.
     */
    struct TransitionMiddle
    {
        coord_t pos; //! position along edge as measure from edge.from.p
        coord_t lower_bead_count;
        TransitionMiddle(coord_t pos, coord_t lower_bead_count)
        : pos(pos), lower_bead_count(lower_bead_count)
        {}
    };

    /*!
     * Auxiliary for referencing one transition along an edge which may contain multiple transitions
     */
    struct TransitionMidRef
    {
        std::unordered_map<edge_t*, std::list<TransitionMiddle>>::iterator pair_it;
        std::list<TransitionMiddle>::iterator transition_it;
        TransitionMidRef(std::unordered_map<edge_t*, std::list<TransitionMiddle>>::iterator pair_it, std::list<TransitionMiddle>::iterator transition_it)
        : pair_it(pair_it)
        , transition_it(transition_it)
        {}
    };

    /*!
     * Represents the location along an edge where the lower or upper end of a transition should be placed.
     */
    struct TransitionEnd
    {
        coord_t pos; //!< position along edge as measure from edge.from.p, where the edge is always the half edge oriented from lower to higher R
        coord_t lower_bead_count;
        bool is_lower_end; //!< whether this is the ed of the transition with lower bead count
        TransitionEnd(coord_t pos, coord_t lower_bead_count, bool is_lower_end)
        : pos(pos), lower_bead_count(lower_bead_count), is_lower_end(is_lower_end)
        {}
    };

    /*!
     * Generate transitions, update the graph with the addeed ribs aka support edges, set the fractional bead count values.
     */
    void generateTransitioningRibs(const BeadingStrategy& beading_strategy);

    /*!
     * Calculate the anchor locations of transitions.
     * These are the locations where the beading_strategy.optimal_bead_count changes value.
     * 
     * A transitions covers a region which covers the transition anchor.
     * 
     * \param beading_strategy The beading strategy.
     * \param edge_to_transitions Mapping from edges to the transitions on those edges
     */
    void generateTransitionMids(const BeadingStrategy& beading_strategy, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions);

    /*!
     * Filter out transition mids if a region of a bead count value which is locally maximal or locally minimal.
     * Dissolve into the surrounding bead count.
     * \param edge_to_transitions Mapping from edges to the transitions on those edges
     * \param beading_strategy The beading strategy.
     */
    void filterTransitionMids(std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions, const BeadingStrategy& beading_strategy);

    /*!
     * Filter out transition mids if a region of a bead count value which is locally maximal or locally minimal.
     * Dissolve into the surrounding bead count.
     * 
     * Recursive function
     * 
     * \param edge_to_start edge pointing to the node from which to start traveling in all directions except along \p edge_to_start
     * \param origin_transition The transition for which we are checking nearby transitions
     * \param traveled_dist the distance traveled before we came to \p edge_to_start.to
     * \param max_dist Maximum distance of a marked region to consider collapsing
     * \param going_up Whether we are traveling in the upward direction as seen from the \p origin_transition. If this doesn't align with the direction according to the R diff on a consecutive edge we know there was a local optimum
     * \param edge_to_transitions Mapping from edges to the transitions on those edges
     * \param beading_strategy The beading strategy.
     * \return whether the origin transition should be dissolved
     */
    std::list<TransitionMidRef> dissolveNearbyTransitions(edge_t* edge_to_start, TransitionMiddle& origin_transition, coord_t traveled_dist, coord_t max_dist, bool going_up, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions, const BeadingStrategy& beading_strategy);

    /*!
     * Dissolve a small region with a (locally extremal) bead count value into the surrounding bead count value.
     */
    void dissolveBeadCountRegion(edge_t* edge_to_start, coord_t from_bead_count, coord_t to_bead_count);

    /*!
     * Dissolve transition anchors for which the transition regions would partly fall outside of the marked region.
     * 
     * Recursive function
     * 
     * \param edge_to_start edge pointing to the node from which to start traveling in all directions except along \p edge_to_start
     * \param traveled_dist the distance traveled before we came to \p edge_to_start.to
     * \param max_dist Maximum distance of a marked region to consider dissolving
     * \param replacing_bead_count Bead count to set in the region.
     * \param beading_strategy The beading strategy.
     */
    bool filterEndOfMarkingTransition(edge_t* edge_to_start, coord_t traveled_dist, coord_t max_dist, coord_t replacing_bead_count, const BeadingStrategy& beading_strategy);

    /*!
     * Generate the end points of the transitions at some distance from the anchors.
     * This function needs to take care of possible branchings in central regions.
     * 
     * \param beading_strategy The beading strategy.
     * \param edge_to_transitions The transition anchors
     * \param[out] edge_to_transition_ends The generated transition end points
     */
    void generateTransitionEnds(const BeadingStrategy& beading_strategy, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions, std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends);

    /*!
     * Apply transitions and introduce required graph edges.
     * 
     * Also set the rest values (fractional bead count) at nodes in between the transition ends
     * 
     * \param edge_to_transition_ends Locations where transition ends should be introduced in the graph
     */
    void applyTransitions(std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends);

    /*!
     * Insert a node into the graph and connect it to the input polygon using ribs
     * 
     * \return the last edge which replaced [edge], which points to the same [to] node
     */
    edge_t* insertNode(edge_t* edge, Point mid, coord_t mide_node_bead_count);

    /*!
     * \param edge_to_transition_mids From the upward halfedges to their trnsitions mids
     */
    void generateTransition(edge_t& edge, coord_t mid_R, const BeadingStrategy& beading_strategy, coord_t transition_lower_bead_count, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transition_mids, std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends);

    /*!
     * \p start_rest and \p end_rest refer to gap distances at the start and end pos in terms of ratios w.r.t. the inner bead width at the high end of the transition
     * 
     * \p end_pos_along_edge may be beyond this edge!
     * In this case we need to interpolate the rest value at the locations in between
     * 
     * \return whether the subgraph is going downward
     */
    bool generateTransitionEnd(edge_t& edge, coord_t start_pos, coord_t end_pos, coord_t transition_half_length, float start_rest, float end_rest, coord_t transition_lower_bead_count, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transition_mids, std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends);

    bool isGoingDown(edge_t* outgoing, coord_t traveled_dist, coord_t transition_half_length, coord_t lower_bead_count, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transition_mids) const;

    /*!
     * Return the first and last edge of the edges replacing \p edge pointing to the same node
     */
    std::pair<SkeletalTrapezoidation::edge_t*, SkeletalTrapezoidation::edge_t*> insertRib(edge_t& edge, node_t* mid_node);

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
     */
    bool isUpward(const edge_t* edge) const;

    /*!
     * Whether the node is connected to any marked edge
     */
    bool isMarked(const node_t* node) const;


    /*!
     * Generate support edges at locations with a certain distance_to_boundary,
     * as prescribed by \ref beading_strategy.getNonlinearThicknesses()
     */
    void generateExtraRibs(const BeadingStrategy& beading_strategy);

//
// ^^^^^^^^^^^^^^^^^^^^^
//    TRANSTISIONING
// =====================
//
// =====================
//  TOOLPATH GENERATION
// vvvvvvvvvvvvvvvvvvvvv
//

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

public:
    void debugCheckDecorationConsistency(bool transitioned); //!< Check logical relationships relting to distance_to_boundary and is_marked etc. Should be true anywhere after setMarking(.)
    void debugCheckTransitionMids(const std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions) const;
    void debugOutput(SVG& svg, std::unordered_map<edge_t*, std::list<TransitionMiddle>>* edge_to_transition_mids = nullptr, std::unordered_map<edge_t*, std::list<TransitionEnd>>* edge_to_transition_ends = nullptr);
    void debugOutput(SVG& svg, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions);
    void debugOutput(STLwriter& stl, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading);
protected:
    SVG::ColorObject getColor(edge_t& edge);

};




} // namespace arachne
#endif // ARACHNE_VARIABLE_WIDTH_INSET_GENERATOR_H
