//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_TRAPEZOIDATION_QUANTIZER_H
#define ARACHNE_TRAPEZOIDATION_QUANTIZER_H

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
#include "SkeletalTrapezoidationNode.h"
#include "BeadingStrategies/BeadingStrategy.h"

#include "utils/STLwriter.h"

namespace arachne
{

/*!
 * Take a skeletal trapezoidation and 'quantize' the heights:
 * assign integer bead counts, and transition between integer bead counts using fractional bead counts in transition regions.
 */
class TrapezoidationQuantizer
{
    using edge_t = SkeletalTrapezoidationEdge;
    using node_t = SkeletalTrapezoidationNode;
    using Beading = BeadingStrategy::Beading;

    SkeletalTrapezoidation& st;
    const BeadingStrategy& beading_strategy;

private:
    const Polygons& polys; //!< input outline boundary shape

    float transitioning_angle; //!< How pointy a region should be before we apply the method. Equals 180* - limit_bisector_angle
    coord_t discretization_step_size; //!< approximate size of segments when parabolic VD edges get discretized (and vertex-vertex edges)
    coord_t transition_filter_dist; //!< filter transition mids (i.e. anchors) closer together than this
    coord_t marking_filter_dist = 20; //!< filter areas marked as 'central' smaller than this
    coord_t snap_dist = 20; //!< generic arithmatic inaccuracy. Only used to determine whether a transition really needs to insert an extra edge.

public:
    TrapezoidationQuantizer(
    SkeletalTrapezoidation& st
    , const BeadingStrategy& beading_strategy
    , const Polygons& polys
    , float transitioning_angle
    , coord_t discretization_step_size = 200
    , coord_t transition_filter_dist = 1000
    );

    void applyBeadCounts(bool filter_outermost_marked_edges);

private:
    
    void setMarking(); //!< set the is_marked flag for each edge based on the transitioning_angle

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
    void setBeadCount();

    /*!
     * Add marked regions and set bead counts
     * where there is an end of marking and when traveling upward we get to another region with the same bead count
     */
    void filterUnmarkedRegions();

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
    bool filterUnmarkedRegions(edge_t* to_edge, coord_t bead_count, coord_t traveled_dist, coord_t max_dist);

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
     * maps the upward edge to the transitions.
     * We only map the halfedge for which the distance_to_boundary is higher at the end than at the beginning
     */
    std::unordered_map<edge_t*, std::list<TransitionMiddle>> edge_to_transition_mids;

    /*!
     * we only map the half edge in the upward direction. mapped items are not sorted
     */
    std::unordered_map<edge_t*, std::list<TransitionEnd>> edge_to_transition_ends;

    /*!
     * Generate transitions, update the graph with the addeed ribs aka support edges, set the fractional bead count values.
     */
    void generateTransitioningRibs();

    /*!
     * Calculate the anchor locations of transitions.
     * These are the locations where the beading_strategy.optimal_bead_count changes value.
     * 
     * A transitions covers a region which covers the transition anchor.
     * 
     * \param beading_strategy The beading strategy.
     * \param edge_to_transition_mids Mapping from edges to the transitions on those edges
     */
    void generateTransitionMids();

    /*!
     * Filter out transition mids if a region of a bead count value which is locally maximal or locally minimal.
     * Dissolve into the surrounding bead count.
     * \param edge_to_transition_mids Mapping from edges to the transitions on those edges
     * \param beading_strategy The beading strategy.
     */
    void filterTransitionMids(std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transition_mids);

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
     * \param edge_to_transition_mids Mapping from edges to the transitions on those edges
     * \param beading_strategy The beading strategy.
     * \return whether the origin transition should be dissolved
     */
    std::list<TransitionMidRef> dissolveNearbyTransitions(edge_t* edge_to_start, TransitionMiddle& origin_transition, coord_t traveled_dist, coord_t max_dist, bool going_up);

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
    bool filterEndOfMarkingTransition(edge_t* edge_to_start, coord_t traveled_dist, coord_t max_dist, coord_t replacing_bead_count);

    /*!
     * Generate the end points of the transitions at some distance from the anchors.
     * This function needs to take care of possible branchings in central regions.
     * 
     * \param beading_strategy The beading strategy.
     * \param edge_to_transition_mids The transition anchors
     * \param[out] edge_to_transition_ends The generated transition end points
     */
    void generateTransitionEnds();

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
    void generateTransition(edge_t& edge, coord_t mid_R, const BeadingStrategy& beading_strategy, coord_t transition_lower_bead_count);

    /*!
     * \p start_rest and \p end_rest refer to gap distances at the start and end pos in terms of ratios w.r.t. the inner bead width at the high end of the transition
     * 
     * \p end_pos_along_edge may be beyond this edge!
     * In this case we need to interpolate the rest value at the locations in between
     * 
     * \return whether the subgraph is going downward
     */
    bool generateTransitionEnd(edge_t& edge, coord_t start_pos, coord_t end_pos, coord_t transition_half_length, float start_rest, float end_rest, coord_t transition_lower_bead_count);

    /*!
     * Whether the quntized bead count is going down in the direction of the \p outgoing edge.
     * 
     * Note that the dist_to_boundary may be going up, but due to the quantization the bead count may stay the same.
     */
    bool isGoingDown(edge_t* outgoing, coord_t traveled_dist, coord_t transition_half_length, coord_t lower_bead_count) const;

    /*!
     * Return the first and last edge of the edges replacing \p edge pointing to the same node
     */
    std::pair<edge_t*, edge_t*> insertRib(edge_t& edge, node_t* mid_node);

    /*!
     * Get the source polygon segment corresponding to the trapezoid in which the \p edge lies.
     * 
     * In case the source geometry is a concave vertex the two points returned are equal.
     */
    std::pair<Point, Point> getSource(const edge_t& edge);
    
    /*!
     * Check whether this edge is marked and no further edges are marked in the edge.next direction
     */
    bool isEndOfMarking(const edge_t& edge) const;

    /*!
     * Generate support edges at locations with a certain distance_to_boundary,
     * as prescribed by \ref beading_strategy.getNonlinearThicknesses()
     */
    void generateExtraRibs();

    // ^ transitioning | v debug

public:
    void debugCheckDecorationConsistency(bool transitioned); //!< Check logical relationships relting to distance_to_boundary and is_marked etc. Should be true anywhere after setMarking(.)
    void debugCheckTransitionMids() const;
    void debugOutput(SVG& svg);
    void debugOutputJunctions(SVG& svg);
protected:
    SVG::ColorObject getColor(edge_t& edge);

};




} // namespace arachne
#endif // ARACHNE_TRAPEZOIDATION_QUANTIZER_H
