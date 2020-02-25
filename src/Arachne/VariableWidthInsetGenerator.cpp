//Copyright (c) 2019 Ultimaker B.V.
#include "VariableWidthInsetGenerator.h"

#include <stack>
#include <functional>
#include <unordered_set>
#include <sstream>
#include <queue>
#include <functional>

#include "BoostInterface.hpp"

#include "utils/VoronoiUtils.h"

#include "utils/linearAlg2D.h"
#include "utils/IntPoint.h"
#include "utils/polygonUtils.h" // only used in svg output for now

#include "utils/logoutput.h"

#include "utils/macros.h"

namespace arachne
{

VariableWidthInsetGenerator::VariableWidthInsetGenerator(const Polygons& polys, float transitioning_angle
, coord_t discretization_step_size
, coord_t transition_filter_dist
, coord_t beading_propagation_transition_dist
)
: SkeletalTrapezoidation(polys, transitioning_angle, discretization_step_size)
, polys(polys)
, transition_filter_dist(transition_filter_dist)
, beading_propagation_transition_dist(beading_propagation_transition_dist)
{
}

std::vector<std::list<ExtrusionLine>> VariableWidthInsetGenerator::generateToolpaths(const BeadingStrategy& beading_strategy, bool filter_outermost_marked_edges)
{
    initialize_graph(); // generate graph
    
    setMarking(beading_strategy);

    filterMarking(marking_filter_dist);

    if (filter_outermost_marked_edges)
    {
        filterOuterMarking();
    }

        debugCheckGraphCompleteness();
        debugCheckGraphConsistency();

    setBeadCount(beading_strategy);

#ifdef DEBUG
    {
        SVG svg("output/unfiltered.svg", AABB(polys));
        SkeletalTrapezoidation::debugOutput(svg, false, false, true, false);
    }
#endif

    filterUnmarkedRegions(beading_strategy);

    debugCheckDecorationConsistency(false);

#ifdef DEBUG
    {
        SVG svg("output/filtered.svg", AABB(polys));
        SkeletalTrapezoidation::debugOutput(svg, false, false, true, false);
    }
#endif

    generateTransitioningRibs(beading_strategy);

    generateExtraRibs(beading_strategy);

#ifdef DEBUG
    {
        AABB aabb(polys);
        SVG svg("output/radial_dists.svg", aabb);
        SkeletalTrapezoidation::debugOutput(svg, false, true);
    }
    {
        AABB aabb(polys);
        SVG svg("output/bead_counts.svg", aabb);
        SkeletalTrapezoidation::debugOutput(svg, false, false, true);
    }
    {
        AABB aabb(polys);
        SVG svg("output/locations.svg", aabb);
        SkeletalTrapezoidation::debugOutput(svg, false, false, false, true);
    }
#endif // DEBUG

    debugCheckDecorationConsistency(true);

    std::vector<std::list<ExtrusionLine>> result_polylines_per_index;
    generateSegments(result_polylines_per_index, beading_strategy);
    // junctions = generateJunctions

    return result_polylines_per_index;
}

void VariableWidthInsetGenerator::setMarking(const BeadingStrategy& beading_strategy)
{
    //                                            _.-'^`      .
    //                                      _.-'^`            .
    //                                _.-'^` \                .
    //                          _.-'^`        \               .
    //                    _.-'^`               \ R2           .
    //              _.-'^` \              _.-'\`\             .
    //        _.-'^`        \R1     _.-'^`     '`\ dR         .
    //  _.-'^`a/2            \_.-'^`a             \           .
    //  `^'-._````````````````A```````````v````````B```````   .
    //        `^'-._                     dD = |AB|            .
    //              `^'-._                                    .
    //                             sin a = dR / dD            .

    coord_t outer_edge_filter_length = beading_strategy.transition_thickness(0) / 2;

    float cap = sin(beading_strategy.transitioning_angle * 0.5); // = cos(bisector_angle / 2)
    for (edge_t& edge : graph.edges)
    {
        assert(edge.twin);
        if (edge.twin->data.markingIsSet())
        {
            edge.data.setMarked(edge.twin->data.isMarked());
        }
        else if (edge.data.type == SkeletalTrapezoidationEdge::EXTRA_VD)
        {
            edge.data.setMarked(false);
        }
        else if (std::max(edge.from->data.distance_to_boundary, edge.to->data.distance_to_boundary) < outer_edge_filter_length)
        {
            edge.data.setMarked(false);
        }
        else
        {
            Point a = edge.from->p;
            Point b = edge.to->p;
            Point ab = b - a;
            coord_t dR = std::abs(edge.to->data.distance_to_boundary - edge.from->data.distance_to_boundary);
            coord_t dD = vSize(ab);
            edge.data.setMarked(dR < dD * cap);
        }
    }
}


void VariableWidthInsetGenerator::filterMarking(coord_t max_length)
{
    for (edge_t& edge : graph.edges)
    {
        if (isEndOfMarking(edge) && !isLocalMaximum(*edge.to) && !isLocalMaximum(*edge.to))
        {
            filterMarking(edge.twin, 0, max_length);
        }
    }
}


bool VariableWidthInsetGenerator::filterMarking(edge_t* starting_edge, coord_t traveled_dist, coord_t max_length)
{
    coord_t length = vSize(starting_edge->from->p - starting_edge->to->p);
    if (traveled_dist + length > max_length)
    {
        return false;
    }
    bool should_dissolve = true;
    for (edge_t* next_edge = starting_edge->next; next_edge && next_edge != starting_edge->twin; next_edge = next_edge->twin->next)
    {
        if (next_edge->data.isMarked())
        {
            should_dissolve &= filterMarking(next_edge, traveled_dist + length, max_length);
        }
    }
    should_dissolve &= !isLocalMaximum(*starting_edge->to); // don't filter marked regions with a local maximum!
    if (should_dissolve)
    {
        starting_edge->data.setMarked(false);
        starting_edge->twin->data.setMarked(false);
    }
    return should_dissolve;
}

void VariableWidthInsetGenerator::filterOuterMarking()
{
    for (edge_t& edge : graph.edges)
    {
        if (!edge.prev)
        {
            edge.data.setMarked(false);
            edge.twin->data.setMarked(false);
        }
    }
}

void VariableWidthInsetGenerator::setBeadCount(const BeadingStrategy& beading_strategy)
{
    for (edge_t& edge : graph.edges)
    {
        if (edge.data.isMarked())
        {
            edge.to->data.bead_count = beading_strategy.optimal_bead_count(edge.to->data.distance_to_boundary * 2);
        }
    }

    // fix bead count at locally maximal R
    // also for marked regions!! See TODO s in generateTransitionEnd(.)
    for (node_t& node : graph.nodes)
    {
        if (isLocalMaximum(node))
        {
            if (node.data.distance_to_boundary < 0)
            {
                RUN_ONCE(logWarning("Distance to boundary not yet computed for local maximum!\n"));
                node.data.distance_to_boundary = std::numeric_limits<coord_t>::max();
                bool first = true;
                for (edge_t* edge = node.some_edge; first || edge != node.some_edge; edge = edge->twin->next)
                {
                    node.data.distance_to_boundary = std::min(node.data.distance_to_boundary, edge->to->data.distance_to_boundary + vSize(edge->from->p - edge->to->p));
                }
            }
            coord_t bead_count = beading_strategy.optimal_bead_count(node.data.distance_to_boundary * 2);
            node.data.bead_count = bead_count;
        }
    }
}

void VariableWidthInsetGenerator:: filterUnmarkedRegions(const BeadingStrategy& beading_strategy)
{
    for (edge_t& edge : graph.edges)
    {
        if (!isEndOfMarking(edge))
        {
            continue;
        }
        assert(edge.to->data.bead_count >= 0 || edge.to->data.distance_to_boundary == 0);
        coord_t max_dist = 400; // beading_strategy.getTransitioningLength(edge.to->data.bead_count)
        filterUnmarkedRegions(&edge, edge.to->data.bead_count, 0, max_dist, beading_strategy);
    }
}

bool VariableWidthInsetGenerator::filterUnmarkedRegions(edge_t* to_edge, coord_t bead_count, coord_t traveled_dist, coord_t max_dist, const BeadingStrategy& beading_strategy)
{
    coord_t r = to_edge->to->data.distance_to_boundary;
    bool dissolve = false;
    for (edge_t* next_edge = to_edge->next; next_edge && next_edge != to_edge->twin; next_edge = next_edge->twin->next)
    {
        coord_t length = vSize(next_edge->to->p - next_edge->from->p);
        if (next_edge->to->data.distance_to_boundary < r && !shorterThen(next_edge->to->p - next_edge->from->p, 10))
        { // only walk upward
            continue;
        }
        if (next_edge->to->data.bead_count == bead_count)
        {
            dissolve = true;
        }
        else if (next_edge->to->data.bead_count < 0)
        {
            dissolve = filterUnmarkedRegions(next_edge, bead_count, traveled_dist + length, max_dist, beading_strategy);
        }
        else // upward bead count is different
        {
            // dissolve if two marked regions with different bead count are closer together than the max_dist (= transition distance)
            dissolve = (traveled_dist + length < max_dist) && std::abs(next_edge->to->data.bead_count - bead_count) == 1;
        }
        if (dissolve)
        {
            next_edge->data.setMarked(true);
            next_edge->twin->data.setMarked(true);
            next_edge->to->data.bead_count = beading_strategy.optimal_bead_count(next_edge->to->data.distance_to_boundary * 2);
            next_edge->to->data.transition_ratio = 0;
        }
        return dissolve; // dissolving only depend on the one edge going upward. There cannot be multiple edges going upward.
    }
    return dissolve;
}

void VariableWidthInsetGenerator::generateTransitioningRibs(const BeadingStrategy& beading_strategy)
{
        debugCheckGraphCompleteness();

    std::unordered_map<edge_t*, std::list<TransitionMiddle>> edge_to_transitions; // maps the upward edge to the transitions. WE only map the halfedge for which the distance_to_boundary is higher at the end than at the beginning
    generateTransitionMids(beading_strategy, edge_to_transitions);

    for (edge_t& edge : graph.edges)
    { // check if there is a transition in between nodes with different bead counts
        if (edge.data.isMarked() && edge.from->data.bead_count != edge.to->data.bead_count)
            assert(edge_to_transitions.find(&edge) != edge_to_transitions.end()
                || edge_to_transitions.find(edge.twin) != edge_to_transitions.end() );
    }
    
        debugCheckGraphCompleteness();
        debugCheckGraphConsistency();

#ifdef DEBUG
    {
        SVG svg("output/transition_mids_unfiltered.svg", AABB(polys));
        SkeletalTrapezoidation::debugOutput(svg, false, false, true, false);
        debugOutput(svg, &edge_to_transitions);
    }
#endif

    filterTransitionMids(edge_to_transitions, beading_strategy);

#ifdef DEBUG
    {
        SVG svg("output/transition_mids.svg", AABB(polys));
        SkeletalTrapezoidation::debugOutput(svg, false, false, true, false);
        debugOutput(svg, &edge_to_transitions);
    }
#endif

    debugCheckTransitionMids(edge_to_transitions);

    std::unordered_map<edge_t*, std::list<TransitionEnd>> edge_to_transition_ends; // we only map the half edge in the upward direction. mapped items are not sorted
    generateTransitionEnds(beading_strategy, edge_to_transitions, edge_to_transition_ends);

#ifdef DEBUG
    {
        SVG svg("output/transition_ends.svg", AABB(polys));
        SkeletalTrapezoidation::debugOutput(svg, false, false, true, false);
        debugOutput(svg, &edge_to_transitions, &edge_to_transition_ends);
    }
#endif

    applyTransitions(edge_to_transition_ends);
}


void VariableWidthInsetGenerator::generateTransitionMids(const BeadingStrategy& beading_strategy, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions)
{
    for (edge_t& edge : graph.edges)
    {
        assert(edge.data.markingIsSet());
        if (!edge.data.isMarked())
        { // only marked regions introduce transitions
            continue;
        }
        coord_t start_R = edge.from->data.distance_to_boundary;
        coord_t end_R = edge.to->data.distance_to_boundary;
        coord_t start_bead_count = edge.from->data.bead_count;
        coord_t end_bead_count = edge.to->data.bead_count;

        if (start_R == end_R)
        { // no transitions occur when both end points have the same distance_to_boundary
            assert(edge.from->data.bead_count == edge.to->data.bead_count);// TODO: what to do in this case?
            continue;
        }
        else if (start_R > end_R)
        { // only consider those half-edges which are going from a lower to a higher distance_to_boundary
            continue;
        }

        if (edge.from->data.bead_count == edge.to->data.bead_count)
        { // no transitions should accur according to the enforced bead counts
            continue;
        }

        if (start_bead_count > beading_strategy.optimal_bead_count(start_R * 2)
            || end_bead_count > beading_strategy.optimal_bead_count(end_R * 2))
        { // wasn't the case earlier in this function because of already introduced transitions
            RUN_ONCE(logError("transitioning segment overlap! (?)\n"));
        }
        assert(start_R < end_R);
        coord_t edge_size = vSize(edge.from->p - edge.to->p);
        for (coord_t transition_lower_bead_count = start_bead_count; transition_lower_bead_count < end_bead_count; transition_lower_bead_count++)
        {
            coord_t mid_R = beading_strategy.transition_thickness(transition_lower_bead_count) / 2;
            if (mid_R > end_R)
            {
                RUN_ONCE(logError("transition on segment lies outside of segment!\n"));
                mid_R = end_R;
            }
            if (mid_R < start_R)
            {
                RUN_ONCE(logError("transition on segment lies outside of segment!\n"));
                mid_R = start_R;
            }
            coord_t mid_pos = edge_size * (mid_R - start_R) / (end_R - start_R);
            assert(mid_pos >= 0);
            assert(mid_pos <= edge_size);
            assert(edge_to_transitions[&edge].empty() || mid_pos >= edge_to_transitions[&edge].back().pos);
            edge_to_transitions[&edge].emplace_back(mid_pos, transition_lower_bead_count);
        }
        if (edge.from->data.bead_count != edge.to->data.bead_count)
        {
            assert(edge_to_transitions[&edge].size() >= 1);
        }
    }
}

void VariableWidthInsetGenerator::filterTransitionMids(std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions, const BeadingStrategy& beading_strategy)
{
    for (auto pair_it = edge_to_transitions.begin(); pair_it != edge_to_transitions.end();)
    {
        std::pair<edge_t* const, std::list<TransitionMiddle>>& pair = *pair_it;
        edge_t* edge = pair.first;
        std::list<TransitionMiddle>& transitions = pair.second;
        if (transitions.empty())
        {
            pair_it = edge_to_transitions.erase(pair_it);
            continue;
        }
        assert(transitions.front().lower_bead_count <= transitions.back().lower_bead_count); // this is how stuff should be stored in edge_to_transitions
        assert(edge->from->data.distance_to_boundary <= edge->to->data.distance_to_boundary); // this is how stuff should be stored in edge_to_transitions
        Point a = edge->from->p;
        Point b = edge->to->p;
        Point ab = b - a;
        coord_t ab_size = vSize(ab);

        bool going_up = true;
        std::list<TransitionMidRef> to_be_dissolved_back = dissolveNearbyTransitions(edge, transitions.back(), ab_size - transitions.back().pos, transition_filter_dist, going_up, edge_to_transitions, beading_strategy);
        bool should_dissolve_back = !to_be_dissolved_back.empty();
        for (TransitionMidRef& ref : to_be_dissolved_back)
        {
            dissolveBeadCountRegion(edge, transitions.back().lower_bead_count + 1, transitions.back().lower_bead_count);
            if (ref.pair_it->second.size() <= 1)
            {
                edge_to_transitions.erase(ref.pair_it);
            }
            else
            {
                ref.pair_it->second.erase(ref.transition_it);
            }
        }

        {
            coord_t trans_bead_count = transitions.back().lower_bead_count;
            coord_t upper_transition_half_length = (1.0 - beading_strategy.getTransitionAnchorPos(trans_bead_count)) * beading_strategy.getTransitioningLength(trans_bead_count);
            should_dissolve_back |= filterEndOfMarkingTransition(edge, ab_size - transitions.back().pos, upper_transition_half_length, trans_bead_count, beading_strategy);
        }
        if (should_dissolve_back)
        {
            transitions.pop_back();
        }
        if (transitions.empty())
        { // filterEndOfMarkingTransition gives inconsistent new bead count when executing for the same transition in two directions.
            pair_it = edge_to_transitions.erase(pair_it);
            continue;
        }

        going_up = false;
        std::list<TransitionMidRef> to_be_dissolved_front = dissolveNearbyTransitions(edge->twin, transitions.front(), transitions.front().pos, transition_filter_dist, going_up, edge_to_transitions, beading_strategy);
        bool should_dissolve_front = !to_be_dissolved_front.empty();
        for (TransitionMidRef& ref : to_be_dissolved_front)
        {
            dissolveBeadCountRegion(edge->twin, transitions.front().lower_bead_count, transitions.front().lower_bead_count + 1);
            if (ref.pair_it->second.size() <= 1)
            {
                edge_to_transitions.erase(ref.pair_it);
            }
            else
            {
                ref.pair_it->second.erase(ref.transition_it);
            }
        }

        {
            coord_t trans_bead_count = transitions.front().lower_bead_count;
            coord_t lower_transition_half_length = beading_strategy.getTransitionAnchorPos(trans_bead_count) * beading_strategy.getTransitioningLength(trans_bead_count);
            should_dissolve_front |= filterEndOfMarkingTransition(edge->twin, transitions.front().pos, lower_transition_half_length, trans_bead_count + 1, beading_strategy);
        }
        if (should_dissolve_front)
        {
            transitions.pop_front();
        }
        if (transitions.empty())
        { // filterEndOfMarkingTransition gives inconsistent new bead count when executing for the same transition in two directions.
            pair_it = edge_to_transitions.erase(pair_it);
            continue;
        }
        ++pair_it; // normal update of loop
    }
}

std::list<VariableWidthInsetGenerator::TransitionMidRef> VariableWidthInsetGenerator::dissolveNearbyTransitions(edge_t* edge_to_start, TransitionMiddle& origin_transition, coord_t traveled_dist, coord_t max_dist, bool going_up, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions, const BeadingStrategy& beading_strategy)
{
    std::list<TransitionMidRef> to_be_dissolved;
    if (traveled_dist > max_dist)
    {
        return to_be_dissolved;
    }
    bool should_dissolve = true;
    for (edge_t* edge = edge_to_start->next; edge && edge != edge_to_start->twin; edge = edge->twin->next)
    {
        if (!edge->data.isMarked())
        {
            continue;
        }
        Point a = edge->from->p;
        Point b = edge->to->p;
        Point ab = b - a;
        coord_t ab_size = vSize(ab);
        bool is_aligned = isUpward(edge);
        edge_t* aligned_edge = is_aligned? edge : edge->twin;
        bool seen_transition_on_this_edge = false;
        auto edge_transitions_it = edge_to_transitions.find(aligned_edge);
        if (edge_transitions_it != edge_to_transitions.end())
        {
            std::list<TransitionMiddle>& transitions = edge_transitions_it->second;
            for (auto transition_it = transitions.begin(); transition_it != transitions.end(); ++ transition_it)
            { // note: this is not neccesarily iterating in the traveling direction!
                // check whether we should dissolve
                coord_t pos = is_aligned? transition_it->pos : ab_size - transition_it->pos;
                if (traveled_dist + pos < max_dist
                    && transition_it->lower_bead_count == origin_transition.lower_bead_count) // only dissolve local optima
                {
                    if (traveled_dist + pos < beading_strategy.getTransitioningLength(transition_it->lower_bead_count))
                    {
                        assert(going_up != is_aligned || transition_it->lower_bead_count == 0); // consecutive transitions both in/decreasing in bead count should never be closer together than the transition distance
                    }
                    to_be_dissolved.emplace_back(edge_transitions_it, transition_it);
                    seen_transition_on_this_edge = true;
                }
            }
        }
        if (!seen_transition_on_this_edge)
        {
            std::list<VariableWidthInsetGenerator::TransitionMidRef> to_be_dissolved_here = dissolveNearbyTransitions(edge, origin_transition, traveled_dist + ab_size, max_dist, going_up, edge_to_transitions, beading_strategy);
            if (to_be_dissolved_here.empty())
            { // the region is too long to be dissolved in this direction, so it cannot be dissolved in any direction.
                to_be_dissolved.clear();
                return to_be_dissolved;
            }
            to_be_dissolved.splice(to_be_dissolved.end(), to_be_dissolved_here); // transfer to_be_dissolved_here into to_be_dissolved
            should_dissolve = should_dissolve && !to_be_dissolved.empty();
        }
    }
    if (!should_dissolve)
    {
        to_be_dissolved.clear();
    }
    return to_be_dissolved;
}


void VariableWidthInsetGenerator::dissolveBeadCountRegion(edge_t* edge_to_start, coord_t from_bead_count, coord_t to_bead_count)
{
    assert(from_bead_count != to_bead_count);
    if (edge_to_start->to->data.bead_count != from_bead_count)
    {
        return;
    }
    edge_to_start->to->data.bead_count = to_bead_count;
    for (edge_t* edge = edge_to_start->next; edge && edge != edge_to_start->twin; edge = edge->twin->next)
    {
        if (!edge->data.isMarked())
        {
            continue;
        }
        dissolveBeadCountRegion(edge, from_bead_count, to_bead_count);
    }
}

bool VariableWidthInsetGenerator::filterEndOfMarkingTransition(edge_t* edge_to_start, coord_t traveled_dist, coord_t max_dist, coord_t replacing_bead_count, const BeadingStrategy& beading_strategy)
{
    if (traveled_dist > max_dist)
    {
        return false;
    }
    bool is_end_of_marking = true;
    bool should_dissolve = false;
    for (edge_t* next_edge = edge_to_start->next; next_edge && next_edge != edge_to_start->twin; next_edge = next_edge->twin->next)
    {
        if (next_edge->data.isMarked())
        {
            coord_t length = vSize(next_edge->to->p - next_edge->from->p);
            should_dissolve |= filterEndOfMarkingTransition(next_edge, traveled_dist + length, max_dist, replacing_bead_count, beading_strategy);
            is_end_of_marking = false;
        }
    }
    if (is_end_of_marking && traveled_dist < max_dist)
    {
        should_dissolve = true;
    }
    if (should_dissolve)
    {
        edge_to_start->to->data.bead_count = replacing_bead_count;
    }
    return should_dissolve;
}

void VariableWidthInsetGenerator::generateTransitionEnds(const BeadingStrategy& beading_strategy, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions, std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends)
{
    for (std::pair<edge_t*, std::list<TransitionMiddle>> pair : edge_to_transitions)
    {
        edge_t* edge = pair.first;
        std::list<TransitionMiddle>& transition_positions = pair.second;

        assert(edge->from->data.distance_to_boundary <= edge->to->data.distance_to_boundary);
        for (TransitionMiddle& transition_middle : transition_positions)
        {
            assert(transition_positions.front().pos <= transition_middle.pos);
            assert(transition_middle.pos <= transition_positions.back().pos);
            generateTransition(*edge, transition_middle.pos, beading_strategy, transition_middle.lower_bead_count, edge_to_transitions, edge_to_transition_ends);
        }
    }
}

void VariableWidthInsetGenerator::generateTransition(edge_t& edge, coord_t mid_pos, const BeadingStrategy& beading_strategy, coord_t lower_bead_count, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transition_mids, std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends)
{
    Point a = edge.from->p;
    Point b = edge.to->p;
    Point ab = b - a;
    coord_t ab_size = vSize(ab);

    coord_t transition_length = beading_strategy.getTransitioningLength(lower_bead_count);
    float transition_mid_position = beading_strategy.getTransitionAnchorPos(lower_bead_count);
    float inner_bead_width_ratio_after_transition = 1.0;

    coord_t start_rest = 0;
    float mid_rest = transition_mid_position * inner_bead_width_ratio_after_transition;
    float end_rest = inner_bead_width_ratio_after_transition;

    
        debugCheckGraphCompleteness();
        debugCheckGraphConsistency();

    { // lower bead count transition end
        coord_t start_pos = ab_size - mid_pos;
        coord_t transition_half_length = transition_mid_position * transition_length;
        coord_t end_pos = start_pos + transition_half_length;
        generateTransitionEnd(*edge.twin, start_pos, end_pos, transition_half_length, mid_rest, start_rest, lower_bead_count, edge_to_transition_mids, edge_to_transition_ends);
    }
    debugCheckGraphConsistency();
    { // upper bead count transition end
        coord_t start_pos = mid_pos;
        coord_t transition_half_length = (1.0 - transition_mid_position) * transition_length;
        coord_t end_pos = mid_pos +  transition_half_length;
        bool is_going_down_everywhere = generateTransitionEnd(edge, start_pos, end_pos, transition_half_length, mid_rest, end_rest, lower_bead_count, edge_to_transition_mids, edge_to_transition_ends);
        assert(!is_going_down_everywhere && "There must have been at least one direction in which the bead count is increasing enough for the transition to happen!");
    }

        debugCheckGraphCompleteness();
        debugCheckGraphConsistency();
}

bool VariableWidthInsetGenerator::generateTransitionEnd(edge_t& edge, coord_t start_pos, coord_t end_pos, coord_t transition_half_length, float start_rest, float end_rest, coord_t lower_bead_count, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transition_mids, std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends)
{
    Point a = edge.from->p;
    Point b = edge.to->p;
    Point ab = b - a;
    coord_t ab_size = vSize(ab); // TODO: prevent recalculation of these values

    assert(start_pos <= ab_size);

    bool going_up = end_rest > start_rest;

    assert(edge.data.isMarked());
    if (!edge.data.isMarked())
    { // This function shouldn't generate ends in or beyond unmarked regions
        return false;
    }

    if (end_pos > ab_size)
    { // recurse on all further edges
        coord_t R = edge.to->data.distance_to_boundary;
        float rest = end_rest - (start_rest - end_rest) * (end_pos - ab_size) / (start_pos - end_pos);
        assert(rest >= 0);
        assert(rest <= std::max(end_rest, start_rest));
        assert(rest >= std::min(end_rest, start_rest));

        coord_t marked_edge_count = 0;
        for (edge_t* outgoing = edge.next; outgoing && outgoing != edge.twin; outgoing = outgoing->twin->next)
        {
            if (!outgoing->data.isMarked()) continue;
            marked_edge_count++;
        }

        bool is_only_going_down = true;
        bool has_recursed = false;
        for (edge_t* outgoing = edge.next; outgoing && outgoing != edge.twin;)
        {
            edge_t* next = outgoing->twin->next; // before we change the outgoing edge itself
            if (!outgoing->data.isMarked())
            {
                outgoing = next;
                continue; // don't put transition ends in non-marked regions
            }
            if (marked_edge_count > 1 && going_up && isGoingDown(outgoing, 0, end_pos - ab_size + transition_half_length, lower_bead_count, edge_to_transition_mids))
            { // we're after a 3-way_all-marked_junction-node and going in the direction of lower bead count
                // don't introduce a transition end along this marked direction, because this direction is the downward direction
                // while we are supposed to be [going_up]
                outgoing = next;
                continue;
            }
            bool is_going_down = generateTransitionEnd(*outgoing, 0, end_pos - ab_size, transition_half_length, rest, end_rest, lower_bead_count, edge_to_transition_mids, edge_to_transition_ends);
            is_only_going_down &= is_going_down;
            outgoing = next;
            has_recursed = true;
        }
        if (!going_up || (has_recursed && !is_only_going_down))
        {
            edge.to->data.transition_ratio = rest;
            edge.to->data.bead_count = lower_bead_count;
        }
        return is_only_going_down;
    }
    else // end_pos < ab_size
    { // add transition end point here
//         assert(edge.data.isMarked() && "we should only be adding transition ends in marked regions");
        
        bool is_lower_end = end_rest == 0; // TODO collapse this parameter into the bool for which it is used here!
        std::list<TransitionEnd>* transitions = nullptr;
        coord_t pos = -1;
        if (isUpward(&edge))
        {
            transitions = &edge_to_transition_ends[&edge];
            pos = end_pos;
        }
        else
        {
            transitions = &edge_to_transition_ends[edge.twin];
            pos = ab_size - end_pos;
        }
        assert(ab_size == vSize(edge.twin->from->p - edge.twin->to->p));
        assert(pos <= ab_size);
        if (transitions->empty() || pos < transitions->front().pos)
        { // preorder so that sorting later on is faster
            transitions->emplace_front(pos, lower_bead_count, is_lower_end);
        }
        else
        {
            transitions->emplace_back(pos, lower_bead_count, is_lower_end);
        }
        return false;
    }
}


bool VariableWidthInsetGenerator::isGoingDown(edge_t* outgoing, coord_t traveled_dist, coord_t max_dist, coord_t lower_bead_count, std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transition_mids) const
{
    // NOTE: the logic below is not fully thought through.
    // TODO: take transition mids into account
    if (outgoing->to->data.distance_to_boundary == 0)
    {
        return true;
    }
    bool is_upward = outgoing->to->data.distance_to_boundary >= outgoing->from->data.distance_to_boundary;
    edge_t* upward_edge = is_upward? outgoing : outgoing->twin;
    if (outgoing->to->data.bead_count > lower_bead_count + 1)
    {
        assert(edge_to_transition_mids.find(upward_edge) != edge_to_transition_mids.end() && "If the bead count is going down there has to be a transition mid!");
        return false;
    }
    coord_t length = vSize(outgoing->to->p - outgoing->from->p);
    auto transition_mids_it = edge_to_transition_mids.find(upward_edge);
    if (transition_mids_it != edge_to_transition_mids.end())
    {
        std::list<TransitionMiddle>& transition_mids = transition_mids_it->second;
        if (!transition_mids.empty())
        {
            TransitionMiddle& mid = is_upward? transition_mids.front() : transition_mids.back();
            if (
                mid.lower_bead_count == lower_bead_count &&
                ((is_upward && mid.pos + traveled_dist < max_dist)
                 || (!is_upward && length - mid.pos + traveled_dist < max_dist))
            )
            {
                return true;
            }
        }
    }
    if (traveled_dist + length > max_dist)
    {
        return false;
    }
    if (outgoing->to->data.bead_count <= lower_bead_count
        && !(outgoing->to->data.bead_count == lower_bead_count && outgoing->to->data.transition_ratio > 0.0))
    {
        return true;
    }
    
    bool is_only_going_down = true;
    bool has_recursed = false;
    for (edge_t* next = outgoing->next; next && next != outgoing->twin; next = next->twin->next)
    {
        if (!next->data.isMarked())
        {
            continue;
        }
        bool is_going_down = isGoingDown(next, traveled_dist + length, max_dist, lower_bead_count, edge_to_transition_mids);
        is_only_going_down &= is_going_down;
        has_recursed = true;
    }
    return has_recursed && is_only_going_down;
}

void VariableWidthInsetGenerator::applyTransitions(std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends)
{
    for (std::pair<edge_t* const, std::list<TransitionEnd>>& pair : edge_to_transition_ends)
    {
        edge_t* edge = pair.first;
        auto twin_ends_it = edge_to_transition_ends.find(edge->twin);
        if (twin_ends_it != edge_to_transition_ends.end())
        {
            coord_t length = vSize(edge->from->p - edge->to->p);
            for (TransitionEnd& end : twin_ends_it->second)
            {
                pair.second.emplace_back(length - end.pos, end.lower_bead_count, end.is_lower_end);
            }
            edge_to_transition_ends.erase(twin_ends_it);
        }
    }
    for (std::pair<edge_t* const, std::list<TransitionEnd>>& pair : edge_to_transition_ends)
    {
        edge_t* edge = pair.first;
        assert(edge->data.isMarked());

        std::list<TransitionEnd>& transitions = pair.second;

        transitions.sort([](const TransitionEnd& a, const TransitionEnd& b) { return a.pos < b.pos; } );

        node_t* from = edge->from;
        node_t* to = edge->to;
        Point a = from->p;
        Point b = to->p;
        Point ab = b - a;
        coord_t ab_size = vSize(ab);

        edge_t* last_edge_replacing_input = edge;
        for (TransitionEnd& transition_end : transitions)
        {
            coord_t new_node_bead_count = transition_end.is_lower_end? transition_end.lower_bead_count : transition_end.lower_bead_count + 1;
            coord_t end_pos = transition_end.pos;
            node_t* close_node = (end_pos < ab_size / 2)? from : to;
            if ((end_pos < snap_dist || end_pos > ab_size - snap_dist)
                && close_node->data.bead_count == new_node_bead_count
            )
            {
                assert(end_pos <= ab_size);
//                 close_node->data.bead_count = new_node_bead_count;
                close_node->data.transition_ratio = 0;
                continue;
            }
            Point mid = a + normal(ab, end_pos);
            
            debugCheckGraphCompleteness();
            debugCheckGraphConsistency();

            debugCheckDecorationConsistency(false);

            assert(last_edge_replacing_input->data.isMarked());
            assert(last_edge_replacing_input->data.type != SkeletalTrapezoidationEdge::EXTRA_VD);
            last_edge_replacing_input = insertNode(last_edge_replacing_input, mid, new_node_bead_count);
            assert(last_edge_replacing_input->data.type != SkeletalTrapezoidationEdge::EXTRA_VD);
            assert(last_edge_replacing_input->data.isMarked());


            debugCheckGraphCompleteness();
            debugCheckGraphConsistency();
        }
    }
}

VariableWidthInsetGenerator::edge_t* VariableWidthInsetGenerator::insertNode(edge_t* edge, Point mid, coord_t mide_node_bead_count)
{
    edge_t* last_edge_replacing_input = edge;

    graph.nodes.emplace_back(SkeletalTrapezoidationJoint(), mid);
    node_t* mid_node = &graph.nodes.back();

    edge_t* twin = last_edge_replacing_input->twin;
    last_edge_replacing_input->twin = nullptr;
    twin->twin = nullptr;
    std::pair<VariableWidthInsetGenerator::edge_t*, VariableWidthInsetGenerator::edge_t*> left_pair
        = insertRib(*last_edge_replacing_input, mid_node);
    std::pair<VariableWidthInsetGenerator::edge_t*, VariableWidthInsetGenerator::edge_t*> right_pair
        = insertRib(*twin, mid_node);
    edge_t* first_edge_replacing_input = left_pair.first;
    last_edge_replacing_input = left_pair.second;
    edge_t* first_edge_replacing_twin = right_pair.first;
    edge_t* last_edge_replacing_twin = right_pair.second;

    first_edge_replacing_input->twin = last_edge_replacing_twin;
    last_edge_replacing_twin->twin = first_edge_replacing_input;
    last_edge_replacing_input->twin = first_edge_replacing_twin;
    first_edge_replacing_twin->twin = last_edge_replacing_input;

    mid_node->data.bead_count = mide_node_bead_count;

    return last_edge_replacing_input;
}

std::pair<VariableWidthInsetGenerator::edge_t*, VariableWidthInsetGenerator::edge_t*> VariableWidthInsetGenerator::insertRib(edge_t& edge, node_t* mid_node)
{
    debugCheckGraphConsistency();
    edge_t* edge_before = edge.prev;
    edge_t* edge_after = edge.next;
    node_t* node_before = edge.from;
    node_t* node_after = edge.to;
    
    Point p = mid_node->p;

    std::pair<Point, Point> source_segment = getSource(edge);
    Point px = LinearAlg2D::getClosestOnLineSegment(p, source_segment.first, source_segment.second);
    coord_t dist = vSize(p - px);
    assert(dist > 0);
    mid_node->data.distance_to_boundary = dist;
    mid_node->data.transition_ratio = 0; // both transition end should have rest = 0, because at the ends a whole number of beads fits without rest

    graph.nodes.emplace_back(SkeletalTrapezoidationJoint(), px);
    node_t* source_node = &graph.nodes.back();
    source_node->data.distance_to_boundary = 0;

    edge_t* first = &edge;
    graph.edges.emplace_back(SkeletalTrapezoidationEdge());
    edge_t* second = &graph.edges.back();
    graph.edges.emplace_back(SkeletalTrapezoidationEdge(SkeletalTrapezoidationEdge::TRANSITION_END));
    edge_t* outward_edge = &graph.edges.back();
    graph.edges.emplace_back(SkeletalTrapezoidationEdge(SkeletalTrapezoidationEdge::TRANSITION_END));
    edge_t* inward_edge = &graph.edges.back();

    if (edge_before) edge_before->next = first;
    first->next = outward_edge;
    outward_edge->next = nullptr;
    inward_edge->next = second;
    second->next = edge_after;

    if (edge_after) edge_after->prev = second;
    second->prev = inward_edge;
    inward_edge->prev = nullptr;
    outward_edge->prev = first;
    first->prev = edge_before;

    first->to = mid_node;
    outward_edge->to = source_node;
    inward_edge->to = mid_node;
    second->to = node_after;

    first->from = node_before;
    outward_edge->from = mid_node;
    inward_edge->from = source_node;
    second->from = mid_node;

    node_before->some_edge = first;
    mid_node->some_edge = outward_edge;
    source_node->some_edge = inward_edge;
    if (edge_after) node_after->some_edge = edge_after;

    first->data.setMarked(true);
    outward_edge->data.setMarked(false); // TODO verify this is always the case.
    inward_edge->data.setMarked(false);
    second->data.setMarked(true);

    outward_edge->twin = inward_edge;
    inward_edge->twin = outward_edge;

    first->twin = nullptr; // we don't know these yet!
    second->twin = nullptr;

    assert(second->prev->from->data.distance_to_boundary == 0);

    debugCheckGraphConsistency();

    return std::make_pair(first, second);
}
std::pair<Point, Point> VariableWidthInsetGenerator::getSource(const edge_t& edge)
{
    const edge_t* from_edge;
    for (from_edge = &edge; from_edge->prev; from_edge = from_edge->prev) {}
    const edge_t* to_edge;
    for (to_edge = &edge; to_edge->next; to_edge = to_edge->next) {}
    return std::make_pair(from_edge->from->p, to_edge->to->p);
}

bool VariableWidthInsetGenerator::isEndOfMarking(const edge_t& edge_to) const
{
    if (!edge_to.data.isMarked())
    {
        return false;
    }
    if (!edge_to.next)
    {
        return true;
    }
    for (const edge_t* edge = edge_to.next; edge && edge != edge_to.twin; edge = edge->twin->next)
    {
        if (edge->data.isMarked())
        {
            return false;
        }
        assert(edge->twin);
    }
    return true;
}

bool VariableWidthInsetGenerator::isLocalMaximum(const node_t& node, bool strict) const
{
    if (node.data.distance_to_boundary == 0)
    {
        return false;
    }
    bool first = true;
    for (edge_t* edge = node.some_edge; first || edge != node.some_edge; edge = edge->twin->next)
    {
        if (canGoUp(edge, strict))
        {
            return false;
        }
        first = false;
        assert(edge->twin); if (!edge->twin) return false;
        if (!edge->twin->next)
        { // This point is on the boundary
            return false;
        }
    }
    return true;
}

bool VariableWidthInsetGenerator::canGoUp(const edge_t* edge, bool strict) const
{
    if (edge->to->data.distance_to_boundary > edge->from->data.distance_to_boundary)
    {
        return true;
    }
    if (edge->to->data.distance_to_boundary < edge->from->data.distance_to_boundary
        || strict
    )
    {
        return false;
    }
    // edge is between equidistqant verts; recurse!
    for (edge_t* outgoing = edge->next; outgoing != edge->twin; outgoing = outgoing->twin->next)
    {
        if (canGoUp(outgoing))
        {
            return true;
        }
        assert(outgoing->twin); if (!outgoing->twin) return false;
        assert(outgoing->twin->next); if (!outgoing->twin->next) return true; // This point is on the boundary?! Should never occur
    }
    return false;
}

std::optional<coord_t> VariableWidthInsetGenerator::distToGoUp(const edge_t* edge) const
{
    if (edge->to->data.distance_to_boundary > edge->from->data.distance_to_boundary)
    {
        return 0;
    }
    if (edge->to->data.distance_to_boundary < edge->from->data.distance_to_boundary)
    {
        return std::optional<coord_t>();
    }
    // edge is between equidistqant verts; recurse!
    std::optional<coord_t> ret;
    for (edge_t* outgoing = edge->next; outgoing != edge->twin; outgoing = outgoing->twin->next)
    {
        std::optional<coord_t> dist_to_up = distToGoUp(outgoing);
        if (dist_to_up)
        {
            if (ret)
            {
                ret = std::min(*ret, *dist_to_up);
            }
            else
            {
                ret = dist_to_up;
            }
        }
        assert(outgoing->twin); if (!outgoing->twin) return std::optional<coord_t>();
        assert(outgoing->twin->next); if (!outgoing->twin->next) return 0; // This point is on the boundary?! Should never occur
    }
    if (ret)
    {
        ret =  *ret + vSize(edge->to->p - edge->from->p);
    }
    return ret;
}

bool VariableWidthInsetGenerator::isUpward(const edge_t* edge) const
{
    if (edge->to->data.distance_to_boundary > edge->from->data.distance_to_boundary)
    {
        return true;
    }
    if (edge->to->data.distance_to_boundary < edge->from->data.distance_to_boundary)
    {
        return false;
    }
    // equidistant edge case:
    std::optional<coord_t> forward_up_dist = distToGoUp(edge);
    std::optional<coord_t> backward_up_dist = distToGoUp(edge->twin);
    if (forward_up_dist && backward_up_dist)
    {
        return forward_up_dist < backward_up_dist;
    }
    if (forward_up_dist) return true;
    if (backward_up_dist) return false;
    return edge->to->p < edge->from->p; // arbitrary ordering, which returns the opposite for the twin edge
}

bool VariableWidthInsetGenerator::isMarked(const node_t* node) const
{
    bool first = true;
    for (edge_t* edge = node->some_edge; first || edge != node->some_edge; edge = edge->twin->next)
    {
        if (edge->data.isMarked())
        {
            return true;
        }
        first = false;
        assert(edge->twin); if (!edge->twin) return false;
    }
    return false;
}







void VariableWidthInsetGenerator::generateExtraRibs(const BeadingStrategy& beading_strategy)
{
    auto end_edge_it = --graph.edges.end(); // don't check newly introduced edges
    for (auto edge_it = graph.edges.begin(); std::prev(edge_it) != end_edge_it; ++edge_it)
    {
        edge_t& edge = *edge_it;
        if ( ! edge.data.isMarked()) continue;
        if (shorterThen(edge.to->p - edge.from->p, discretization_step_size)) continue;
        if (edge.from->data.distance_to_boundary >= edge.to->data.distance_to_boundary) continue;

        std::vector<coord_t> rib_thicknesses = beading_strategy.getNonlinearThicknesses(edge.from->data.bead_count);

        if (rib_thicknesses.empty()) continue;

        // preload some variables before [edge] gets changed
        node_t* from = edge.from;
        node_t* to = edge.to;
        Point a = from->p;
        Point b = to->p;
        Point ab = b - a;
        coord_t ab_size = vSize(ab);
        coord_t a_R = edge.from->data.distance_to_boundary;
        coord_t b_R = edge.to->data.distance_to_boundary;
        
        edge_t* last_edge_replacing_input = &edge;
        for (coord_t rib_thickness : rib_thicknesses)
        {
            if (rib_thickness / 2 <= a_R) continue;
            if (rib_thickness / 2 >= b_R) break;
            coord_t new_node_bead_count = std::min(edge.from->data.bead_count, edge.to->data.bead_count);
            coord_t end_pos = ab_size * (rib_thickness / 2 - a_R) / (b_R - a_R);
            assert(end_pos > 0);
            assert(end_pos < ab_size);
            node_t* close_node = (end_pos < ab_size / 2)? from : to;
            if ((end_pos < snap_dist || end_pos > ab_size - snap_dist)
                && close_node->data.bead_count == new_node_bead_count
            )
            {
                assert(end_pos <= ab_size);
//                 close_node->data.bead_count = new_node_bead_count;
                close_node->data.transition_ratio = 0;
                continue;
            }
            Point mid = a + normal(ab, end_pos);
            
            debugCheckGraphCompleteness();
            debugCheckGraphConsistency();

            debugCheckDecorationConsistency(false);

            assert(last_edge_replacing_input->data.isMarked());
            assert(last_edge_replacing_input->data.type != SkeletalTrapezoidationEdge::EXTRA_VD);
            last_edge_replacing_input = insertNode(last_edge_replacing_input, mid, new_node_bead_count);
            assert(last_edge_replacing_input->data.type != SkeletalTrapezoidationEdge::EXTRA_VD);
            assert(last_edge_replacing_input->data.isMarked());


            debugCheckGraphCompleteness();
            debugCheckGraphConsistency();
        }
    }
}

//
// ^^^^^^^^^^^^^^^^^^^^^
//    TRANSTISIONING
// =====================
//
// =====================
//  TOOLPATH GENERATION
// vvvvvvvvvvvvvvvvvvvvv
//

void VariableWidthInsetGenerator::generateSegments(std::vector<std::list<ExtrusionLine>>& result_polylines_per_index, const BeadingStrategy& beading_strategy)
{
    std::vector<edge_t*> upward_quad_mids;
    for (edge_t& edge : graph.edges)
    {
        if (edge.prev && edge.next && isUpward(&edge))
        {
            upward_quad_mids.emplace_back(&edge);
        }
    }
    std::sort(upward_quad_mids.begin(), upward_quad_mids.end(), [this](edge_t* a, edge_t* b)
        {
            if (a->to->data.distance_to_boundary == b->to->data.distance_to_boundary)
            { // ordering between two 'upward' edges of the same distance is important when one of the edges is flat and connected to the other
                if (a->from->data.distance_to_boundary == a->to->data.distance_to_boundary
                    && b->from->data.distance_to_boundary == b->to->data.distance_to_boundary)
                {
                    coord_t max = std::numeric_limits<coord_t>::max();
                    coord_t a_dist_from_up = std::min(distToGoUp(a).value_or(max), distToGoUp(a->twin).value_or(max)) - vSize(a->to->p - a->from->p);
                    coord_t b_dist_from_up = std::min(distToGoUp(b).value_or(max), distToGoUp(b->twin).value_or(max)) - vSize(b->to->p - b->from->p);
                    return a_dist_from_up < b_dist_from_up;
                }
                else if (a->from->data.distance_to_boundary == a->to->data.distance_to_boundary)
                {
                    return true; // edge a might be 'above' edge b
                }
                else if (b->from->data.distance_to_boundary == b->to->data.distance_to_boundary)
                {
                    return false; // edge b might be 'above' edge a
                }
                else
                {
                    // ordering is not important
                }
            }
            return a->to->data.distance_to_boundary > b->to->data.distance_to_boundary;
        });
    
    std::unordered_map<node_t*, BeadingPropagation> node_to_beading;
    { // store beading
        for (node_t& node : graph.nodes)
        {
            if (node.data.bead_count <= 0)
            {
                continue;
            }
            if (node.data.transition_ratio == 0)
            {
                auto pair = node_to_beading.emplace(&node, beading_strategy.compute(node.data.distance_to_boundary * 2, node.data.bead_count));
                assert(pair.first->second.beading.total_thickness == node.data.distance_to_boundary * 2);
            }
            else
            {
                Beading low_count_beading = beading_strategy.compute(node.data.distance_to_boundary * 2, node.data.bead_count);
                Beading high_count_beading = beading_strategy.compute(node.data.distance_to_boundary * 2, node.data.bead_count + 1);
                Beading merged = interpolate(low_count_beading, 1.0 - node.data.transition_ratio, high_count_beading);
                node_to_beading.emplace(&node, merged);
                assert(merged.total_thickness == node.data.distance_to_boundary * 2);
            }
        }
    }
    
    propagateBeadingsUpward(upward_quad_mids, node_to_beading, beading_strategy);

    propagateBeadingsDownward(upward_quad_mids, node_to_beading, beading_strategy);
    
    std::unordered_map<edge_t*, std::vector<ExtrusionJunction>> edge_to_junctions; // junctions ordered high R to low R
    generateJunctions(node_to_beading, edge_to_junctions, beading_strategy);


#ifdef DEBUG
    {
        SVG svg("output/junctions.svg", AABB(polys));
        SkeletalTrapezoidation::debugOutput(svg, false, false, true, false);
        debugOutput(svg, edge_to_junctions);
    }
    if (generate_MAT_STL)
    {
        STLwriter stl("output/vq.stl");
        debugOutput(stl, edge_to_junctions, node_to_beading);
    }
#endif

    connectJunctions(edge_to_junctions, result_polylines_per_index);
    
    generateLocalMaximaSingleBeads(node_to_beading, result_polylines_per_index);
}

VariableWidthInsetGenerator::edge_t* VariableWidthInsetGenerator::getQuadMaxRedgeTo(edge_t* quad_start_edge)
{
    assert(quad_start_edge->prev == nullptr);
    assert(quad_start_edge->from->data.distance_to_boundary == 0);
    coord_t max_R = -1;
    edge_t* ret = nullptr;
    for (edge_t* edge = quad_start_edge; edge; edge = edge->next)
    {
        coord_t r = edge->to->data.distance_to_boundary;
        if (r > max_R)
        {
            max_R = r;
            ret = edge;
        }
    }
    if (!ret->next && ret->to->data.distance_to_boundary - 5 < ret->from->data.distance_to_boundary)
    {
        ret = ret->prev;
    }
    assert(ret);
    assert(ret->next);
    return ret;
}

void VariableWidthInsetGenerator::propagateBeadingsUpward(std::vector<edge_t*>& upward_quad_mids, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, const BeadingStrategy& beading_strategy)
{
    for (auto upward_quad_mids_it = upward_quad_mids.rbegin(); upward_quad_mids_it != upward_quad_mids.rend(); ++upward_quad_mids_it)
    {
        edge_t* upward_edge = *upward_quad_mids_it;
        if (upward_edge->to->data.bead_count >= 0)
        { // don't override local beading
            continue;
        }
        auto lower_beading_it = node_to_beading.find(upward_edge->from);
        if (lower_beading_it == node_to_beading.end())
        { // only propagate if we have something to propagate
            continue;
        }
        auto upper_beading_it = node_to_beading.find(upward_edge->to);
        if (upper_beading_it != node_to_beading.end())
        { // only propagate to places where there is place
            continue;
        }
        assert((upward_edge->from->data.distance_to_boundary != upward_edge->to->data.distance_to_boundary || shorterThen(upward_edge->to->p - upward_edge->from->p, marking_filter_dist)) && "zero difference R edges should always be marked");
        BeadingPropagation& lower_beading = lower_beading_it->second;
        coord_t length = vSize(upward_edge->to->p - upward_edge->from->p);
        BeadingPropagation upper_beading = lower_beading;
        upper_beading.dist_to_bottom_source += length;
        upper_beading.is_upward_propagated_only = true;
        auto pair = node_to_beading.emplace(upward_edge->to, upper_beading);
        assert(upper_beading.beading.total_thickness <= upward_edge->to->data.distance_to_boundary * 2);
    }
}

void VariableWidthInsetGenerator::propagateBeadingsDownward(std::vector<edge_t*>& upward_quad_mids, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, const BeadingStrategy& beading_strategy)
{
    for (edge_t* upward_quad_mid : upward_quad_mids)
    {
        // transfer beading information to lower nodes
        if (!upward_quad_mid->data.isMarked())
        {
            // for equidistant edge: propagate from known beading to node with unknown beading
            if (upward_quad_mid->from->data.distance_to_boundary == upward_quad_mid->to->data.distance_to_boundary
                && node_to_beading.find(upward_quad_mid->from) != node_to_beading.end()
                && node_to_beading.find(upward_quad_mid->to) == node_to_beading.end()
            )
            {
                propagateBeadingsDownward(upward_quad_mid->twin, node_to_beading, beading_strategy);
            }
            else
            {
                propagateBeadingsDownward(upward_quad_mid, node_to_beading, beading_strategy);
            }
        }
    }
}

void VariableWidthInsetGenerator::propagateBeadingsDownward(edge_t* edge_to_peak, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, const BeadingStrategy& beading_strategy)
{
    coord_t length = vSize(edge_to_peak->to->p - edge_to_peak->from->p);
    BeadingPropagation& top_beading = getBeading(edge_to_peak->to, node_to_beading, beading_strategy);
    assert(top_beading.beading.total_thickness >= edge_to_peak->to->data.distance_to_boundary * 2);
    assert( ! top_beading.is_upward_propagated_only);

    auto it = node_to_beading.find(edge_to_peak->from);
    if (it == node_to_beading.end())
    { // set new beading if there is no beading associatied with the node yet
        BeadingPropagation propagated_beading = top_beading;
        propagated_beading.dist_from_top_source += length;
        auto pair = node_to_beading.emplace(edge_to_peak->from, propagated_beading);
        assert(propagated_beading.beading.total_thickness >= edge_to_peak->from->data.distance_to_boundary * 2);
        assert(pair.second && "we emplaced something");
#ifdef DEBUG
    {
        auto it = node_to_beading.find(edge_to_peak->from);
        assert(it != node_to_beading.end());
        assert(it->second.beading.total_thickness >= edge_to_peak->to->data.distance_to_boundary * 2);
    }
#endif
    }
    else // if (!it->second.is_finished)
    {
        BeadingPropagation& bottom_beading = it->second;
        coord_t total_dist = top_beading.dist_from_top_source + length + bottom_beading.dist_to_bottom_source;
        float ratio_of_top = static_cast<float>(bottom_beading.dist_to_bottom_source) / std::min(total_dist, beading_propagation_transition_dist);
        ratio_of_top = std::max(0.0f, ratio_of_top);
        if (ratio_of_top >= 1.0)
        {
            bottom_beading = top_beading;
            bottom_beading.dist_from_top_source += length;
        }
        else
        {
            Beading merged_beading = interpolate(top_beading.beading, ratio_of_top, bottom_beading.beading, edge_to_peak->from->data.distance_to_boundary);
            bottom_beading = BeadingPropagation(merged_beading);
            bottom_beading.is_upward_propagated_only = false;
            assert(merged_beading.total_thickness >= edge_to_peak->from->data.distance_to_boundary * 2);
        }
#ifdef DEBUG
    {
        auto it = node_to_beading.find(edge_to_peak->from);
        assert(it != node_to_beading.end());
        assert(it->second.beading.total_thickness >= edge_to_peak->from->data.distance_to_boundary * 2);
    }
#endif
    }
}


VariableWidthInsetGenerator::Beading VariableWidthInsetGenerator::interpolate(const Beading& left, float ratio_left_to_whole, const Beading& right, coord_t switching_radius) const
{
    assert(ratio_left_to_whole >= 0.0 && ratio_left_to_whole <= 1.0);
    Beading ret = interpolate(left, ratio_left_to_whole, right);

    // TODO: don't use toolpath locations past the middle!
    // TODO: stretch bead widths and locations of the higher bead count beading to fit in the left over space
    coord_t next_inset_idx;
    for (next_inset_idx = left.toolpath_locations.size() - 1; next_inset_idx >= 0; next_inset_idx--)
    {
        if (switching_radius > left.toolpath_locations[next_inset_idx])
        {
            break;
        }
    }
    if (next_inset_idx < 0)
    { // there is no next inset, because there is only one
        assert(left.toolpath_locations.empty() || left.toolpath_locations.front() >= switching_radius);
        return ret;
    }
    if (next_inset_idx + 1 == coord_t(left.toolpath_locations.size()))
    { // we cant adjust to fit the next edge because there is no previous one?!
        return ret;
    }
    assert(next_inset_idx < coord_t(left.toolpath_locations.size()));
    assert(left.toolpath_locations[next_inset_idx] <= switching_radius);
    assert(left.toolpath_locations[next_inset_idx + 1] >= switching_radius);
    if (ret.toolpath_locations[next_inset_idx] > switching_radius)
    { // one inset disappeared between left and the merged one
        // solve for ratio f:
        // f*l + (1-f)*r = s
        // f*l + r - f*r = s
        // f*(l-r) + r = s
        // f*(l-r) = s - r
        // f = (s-r) / (l-r)
        float new_ratio = static_cast<float>(switching_radius - right.toolpath_locations[next_inset_idx]) / static_cast<float>(left.toolpath_locations[next_inset_idx] - right.toolpath_locations[next_inset_idx]);
        new_ratio = std::min(1.0, new_ratio + 0.1);
        return interpolate(left, new_ratio, right);
    }
    return ret;
}


VariableWidthInsetGenerator::Beading VariableWidthInsetGenerator::interpolate(const Beading& left, float ratio_left_to_whole, const Beading& right) const
{
    assert(ratio_left_to_whole >= 0.0 && ratio_left_to_whole <= 1.0);
    float ratio_right_to_whole = 1.0 - ratio_left_to_whole;

    Beading ret = (left.total_thickness > right.total_thickness)? left : right;
    for (size_t inset_idx = 0; inset_idx < std::min(left.bead_widths.size(), right.bead_widths.size()); inset_idx++)
    {
        ret.bead_widths[inset_idx] = ratio_left_to_whole * left.bead_widths[inset_idx] + ratio_right_to_whole * right.bead_widths[inset_idx];
        ret.toolpath_locations[inset_idx] = ratio_left_to_whole * left.toolpath_locations[inset_idx] + ratio_right_to_whole * right.toolpath_locations[inset_idx];
    }
    return ret;
}

void VariableWidthInsetGenerator::generateJunctions(std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions, const BeadingStrategy& beading_strategy)
{
    for (edge_t& edge_ : graph.edges)
    {
        edge_t* edge = &edge_;
        if (edge->from->data.distance_to_boundary > edge->to->data.distance_to_boundary)
        { // only consider the upward half-edges
            continue;
        }

        coord_t start_R = edge->to->data.distance_to_boundary; // higher R
        coord_t end_R = edge->from->data.distance_to_boundary; // lower R

        Beading* beading = &getBeading(edge->to, node_to_beading, beading_strategy).beading;
        std::vector<ExtrusionJunction>& ret = edge_to_junctions[edge]; // emplace a new vector
        if ((edge->from->data.bead_count == edge->to->data.bead_count && edge->from->data.bead_count >= 0)
            || end_R >= start_R)
        { // no beads to generate
            continue;
        }

        assert(beading->total_thickness >= edge->to->data.distance_to_boundary * 2);

        Point a = edge->to->p;
        Point b = edge->from->p;
        Point ab = b - a;

        coord_t junction_idx;
        // compute starting junction_idx for this segment
        for (junction_idx = (std::max(size_t(1), beading->toolpath_locations.size()) - 1) / 2; junction_idx >= 0 && junction_idx < coord_t(beading->toolpath_locations.size()); junction_idx--)
        {
            coord_t bead_R = beading->toolpath_locations[junction_idx];
            if (bead_R <= start_R)
            { // junction coinciding with start node is used in this function call
                break;
            }
        }

        // rebustness against odd segments which might lie just slightly outside of the range due to rounding errors
        // not sure if this is really needed (TODO)
        if (junction_idx + 1 < coord_t(beading->toolpath_locations.size())
            && beading->toolpath_locations[junction_idx + 1] <= start_R + 5
            && beading->total_thickness < start_R + 5
        )
        {
            junction_idx++;
        }

        for (; junction_idx >= 0 && junction_idx < coord_t(beading->toolpath_locations.size()); junction_idx--)
        {
            coord_t bead_R = beading->toolpath_locations[junction_idx];
            assert(bead_R >= 0);
            if (bead_R < end_R)
            { // junction coinciding with a node is handled by the next segment
                break;
            }
            Point junction(a + ab * (bead_R - start_R) / (end_R - start_R));
            if (bead_R > start_R - 5)
            { // snap to start node if it is really close, in order to be able to see 3-way intersection later on more robustly
                junction = a;
            }
            ret.emplace_back(junction, beading->bead_widths[junction_idx], junction_idx);
        }
    }
}

const std::vector<ExtrusionJunction>& VariableWidthInsetGenerator::getJunctions(edge_t* edge, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions)
{
    assert(edge->to->data.distance_to_boundary >= edge->from->data.distance_to_boundary);
    auto ret_it = edge_to_junctions.find(edge);
    assert(ret_it != edge_to_junctions.end());
    return ret_it->second;
}


VariableWidthInsetGenerator::BeadingPropagation& VariableWidthInsetGenerator::getBeading(node_t* node, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, const BeadingStrategy& beading_strategy)
{
    auto beading_it = node_to_beading.find(node);
    if (beading_it == node_to_beading.end())
    {
        if (node->data.bead_count == -1)
        { // This bug is due to too small marked edges
            constexpr coord_t nearby_dist = 100; // TODO
            BeadingPropagation* nearest_beading = getNearestBeading(node, nearby_dist, node_to_beading);
            if (nearest_beading) return *nearest_beading;
            // else make a new beading:
            bool has_marked_edge = false;
            bool first = true;
            coord_t dist = std::numeric_limits<coord_t>::max();
            for (edge_t* edge = node->some_edge; edge && (first || edge != node->some_edge); edge = edge->twin->next)
            {
                if (edge->data.isMarked())
                {
                    has_marked_edge = true;
                }
                assert(edge->to->data.distance_to_boundary >= 0);
                dist = std::min(dist, edge->to->data.distance_to_boundary + vSize(edge->to->p - edge->from->p));
                first = false;
            }
            RUN_ONCE(logError("Unknown beading for unmarked node!\n"));
//             assert(false);
            assert(dist != std::numeric_limits<coord_t>::max());
            node->data.bead_count = beading_strategy.optimal_bead_count(dist * 2);
        }
        assert(node->data.bead_count != -1);
        beading_it = node_to_beading.emplace(node, beading_strategy.compute(node->data.distance_to_boundary * 2, node->data.bead_count)).first;
    }
    assert(beading_it != node_to_beading.end());
    return beading_it->second;
}

VariableWidthInsetGenerator::BeadingPropagation* VariableWidthInsetGenerator::getNearestBeading(node_t* node, coord_t max_dist, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading)
{
    struct DistEdge
    {
        edge_t* edge_to;
        coord_t dist;
        DistEdge(edge_t* edge_to, coord_t dist)
        : edge_to(edge_to), dist(dist)
        {}
    };

    auto compare = [](const DistEdge& l, const DistEdge& r) -> bool { return l.dist > r.dist; };
    std::priority_queue<DistEdge, std::vector<DistEdge>, decltype(compare)> further_edges(compare);
    bool first = true;
    for (edge_t* outgoing = node->some_edge; outgoing && (first || outgoing != node->some_edge); outgoing = outgoing->twin->next)
    {
        further_edges.emplace(outgoing, vSize(outgoing->to->p - outgoing->from->p));
        first = false;
    }

    for (coord_t counter = 0; counter < 1000; counter++)
    { // prevent endless recursion
        if (further_edges.empty()) return nullptr;
        DistEdge here = further_edges.top();
        further_edges.pop();
        if (here.dist > max_dist) return nullptr;
        auto it = node_to_beading.find(here.edge_to->to);
        if (it != node_to_beading.end())
        {
            return &it->second;
        }
        else
        { // recurse
            for (edge_t* further_edge = here.edge_to->next; further_edge && further_edge != here.edge_to->twin; further_edge = further_edge->twin->next)
            {
                further_edges.emplace(further_edge, here.dist + vSize(further_edge->to->p - further_edge->from->p));
            }
        }
    }
    return nullptr;
}

void VariableWidthInsetGenerator::connectJunctions(std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions, std::vector<std::list<ExtrusionLine>>& result_polylines_per_index)
{   
    auto getNextQuad = [](edge_t* quad_start)
    {
        edge_t* quad_end = quad_start;
        while (quad_end->next) quad_end = quad_end->next;
        return quad_end->twin;
    };
    
    auto addSegment = [&result_polylines_per_index](ExtrusionJunction& from, ExtrusionJunction& to, bool is_odd, bool force_new_path)
    {
        if (from == to) return;

        size_t inset_idx = from.perimeter_index;
        if (inset_idx >= result_polylines_per_index.size()) result_polylines_per_index.resize(inset_idx + 1);
        assert((result_polylines_per_index[inset_idx].empty() || !result_polylines_per_index[inset_idx].back().junctions.empty()) && "empty extrusion lines should never have been generated");
        if ( ! force_new_path
            && ! result_polylines_per_index[inset_idx].empty()
            && result_polylines_per_index[inset_idx].back().is_odd == is_odd
            && shorterThen(result_polylines_per_index[inset_idx].back().junctions.back().p - to.p, 10)
            && std::abs(result_polylines_per_index[inset_idx].back().junctions.back().w - to.w) < 10
            )
        {
            result_polylines_per_index[inset_idx].back().junctions.push_back(from);
        }
        else if ( ! force_new_path
            && ! result_polylines_per_index[inset_idx].empty()
            && result_polylines_per_index[inset_idx].back().is_odd == is_odd
            && shorterThen(result_polylines_per_index[inset_idx].back().junctions.back().p - from.p, 10)
            && std::abs(result_polylines_per_index[inset_idx].back().junctions.back().w - from.w) < 10
            )
        {
            result_polylines_per_index[inset_idx].back().junctions.push_back(to);
        }
        else
        {
            result_polylines_per_index[inset_idx].emplace_back(inset_idx, is_odd);
            result_polylines_per_index[inset_idx].back().junctions.push_back(from);
            result_polylines_per_index[inset_idx].back().junctions.push_back(to);
        }
    };
    
    std::unordered_set<edge_t*> unprocessed_quad_starts(graph.edges.size() * 5 / 2);
    for (edge_t& edge : graph.edges)
    {
        if (!edge.prev)
        {
            unprocessed_quad_starts.insert(&edge);
        }
    }
    
    std::unordered_set<edge_t*> passed_odd_edges;
    
    while ( ! unprocessed_quad_starts.empty())
    {
        edge_t* poly_domain_start = *unprocessed_quad_starts.begin();
        bool first = true;
        for (edge_t* quad_start = poly_domain_start; first || quad_start != poly_domain_start; quad_start = getNextQuad(quad_start))
        {
            first = false;
            edge_t* quad_end = quad_start; while (quad_end->next) quad_end = quad_end->next;
            edge_t* edge_to_peak = getQuadMaxRedgeTo(quad_start);
            // walk down on both sides and connect junctions
            edge_t* edge_from_peak = edge_to_peak->next; assert(edge_from_peak);
            
            unprocessed_quad_starts.erase(quad_start);
            
            
            
            std::vector<ExtrusionJunction> from_junctions = getJunctions(edge_to_peak, edge_to_junctions);
            std::vector<ExtrusionJunction> to_junctions = getJunctions(edge_from_peak->twin, edge_to_junctions);
            if (edge_to_peak->prev)
            {
                std::vector<ExtrusionJunction> from_prev_junctions = getJunctions(edge_to_peak->prev, edge_to_junctions);
                if (!from_junctions.empty() && !from_prev_junctions.empty() && from_junctions.back().perimeter_index == from_prev_junctions.front().perimeter_index)
                {
                    from_junctions.pop_back();
                }
                from_junctions.reserve(from_junctions.size() + from_prev_junctions.size());
                from_junctions.insert(from_junctions.end(), from_prev_junctions.begin(), from_prev_junctions.end());
                assert(!edge_to_peak->prev->prev);
            }
            if (edge_from_peak->next)
            {
                std::vector<ExtrusionJunction> to_next_junctions = getJunctions(edge_from_peak->next->twin, edge_to_junctions);
                if (!to_junctions.empty() && !to_next_junctions.empty() && to_junctions.back().perimeter_index == to_next_junctions.front().perimeter_index)
                {
                    to_junctions.pop_back();
                }
                to_junctions.reserve(to_junctions.size() + to_next_junctions.size());
                to_junctions.insert(to_junctions.end(), to_next_junctions.begin(), to_next_junctions.end());
                assert(!edge_from_peak->next->next);
            }
            assert(std::abs(int(from_junctions.size()) - int(to_junctions.size())) <= 1); // at transitions one end has more beads


            size_t segment_count = std::min(from_junctions.size(), to_junctions.size());
            for (size_t junction_rev_idx = 0; junction_rev_idx < segment_count; junction_rev_idx++)
            {
                ExtrusionJunction& from = from_junctions[from_junctions.size() - 1 - junction_rev_idx];
                ExtrusionJunction& to = to_junctions[to_junctions.size() - 1 - junction_rev_idx];
                assert(from.perimeter_index == to.perimeter_index);
                bool is_odd_segment = edge_to_peak->to->data.bead_count > 0 && edge_to_peak->to->data.bead_count % 2 == 1 // quad contains single bead segment
                    && edge_to_peak->to->data.transition_ratio == 0 && edge_to_peak->from->data.transition_ratio == 0 && edge_from_peak->to->data.transition_ratio == 0 // we're not in a transition
                    && junction_rev_idx == segment_count - 1 // is single bead segment
                    && shorterThen(from.p - quad_start->to->p, 5) && shorterThen(to.p - quad_end->from->p, 5);
                if (is_odd_segment
                    && passed_odd_edges.count(quad_start->next->twin) > 0) // only generate toolpath for odd segments once
                {
                    continue; // prevent duplication of single bead segments
                }
                passed_odd_edges.emplace(quad_start->next);
                bool force_new_path = is_odd_segment && isMultiIntersection(quad_start->to);
                addSegment(from, to, is_odd_segment, force_new_path);
            }
        }
    }
}

bool VariableWidthInsetGenerator::isMultiIntersection(node_t* node)
{
    int odd_path_count = 0;
    bool first = true;
    for (edge_t* outgoing = node->some_edge; first || outgoing != node->some_edge; outgoing = outgoing->twin->next)
    {
        first = false;
        if (outgoing->data.isMarked())
            odd_path_count++;
    }
    return odd_path_count > 2;
}

void VariableWidthInsetGenerator::generateLocalMaximaSingleBeads(std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, std::vector<std::list<ExtrusionLine>>& result_polylines_per_index)
{
    for (auto pair : node_to_beading)
    {
        node_t* node = pair.first;
        Beading& beading = pair.second.beading;
        if (beading.bead_widths.size() % 2 == 1
            && isLocalMaximum(*node, true)
            && !isMarked(node)
        )
        {
            size_t inset_index = beading.bead_widths.size() / 2;
            bool is_odd = true;
            if (inset_index >= result_polylines_per_index.size()) result_polylines_per_index.resize(inset_index + 1);
            result_polylines_per_index[inset_index].emplace_back(inset_index, is_odd);
            ExtrusionLine& line = result_polylines_per_index[inset_index].back();
            line.junctions.emplace_back(node->p, beading.bead_widths[inset_index], inset_index);
            line.junctions.emplace_back(node->p + Point(50, 0), beading.bead_widths[inset_index], inset_index);
        }
    }
}

//
// ^^^^^^^^^^^^^^^^^^^^^
//  TOOLPATH GENERATION
// =====================
//
// =====================
//       HELPERS
// vvvvvvvvvvvvvvvvvvvvv
//

void VariableWidthInsetGenerator::debugCheckDecorationConsistency(bool transitioned)
{
#ifdef DEBUG
    for (const edge_t& edge : graph.edges)
    {
        const edge_t* edge_p = &edge;
        assert(edge.data.type >= SkeletalTrapezoidationEdge::NORMAL);
        if (edge.data.type != SkeletalTrapezoidationEdge::NORMAL)
        {
            if (edge.from->data.distance_to_boundary != -1 && edge.to->data.distance_to_boundary != -1)
            {
                assert(edge.from->data.distance_to_boundary == 0 || edge.to->data.distance_to_boundary == 0);
            }
            assert(!edge.data.isMarked());
        }
        assert(edge.data.isMarked() == edge.twin->data.isMarked());
        if (edge.data.isMarked())
        {
            if (transitioned && edge.from->data.bead_count != -1 && edge.to->data.bead_count != -1)
            {
                assert(!edge.data.isMarked() || std::abs(edge.from->data.bead_count - edge.to->data.bead_count) <= 1);
            }
        }
    }
#endif // DEBUG
}

void VariableWidthInsetGenerator::debugCheckTransitionMids(const std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transitions) const
{
#ifdef DEBUG
    for (std::pair<edge_t*, std::list<TransitionMiddle>> pair : edge_to_transitions)
    {
        const edge_t* edge = pair.first;
        const std::list<TransitionMiddle>& transition_positions = pair.second;
        
        assert(edge->from->data.distance_to_boundary <= edge->to->data.distance_to_boundary);
        
        const TransitionMiddle* prev = nullptr;
        for (const TransitionMiddle& here : transition_positions)
        {
            if (prev)
            {
                assert(here.pos > prev->pos);
                assert(here.lower_bead_count > prev->lower_bead_count);
                assert(std::abs(here.lower_bead_count - prev->lower_bead_count) == 1);
            }
            prev = &here;
        }
        
    }
#endif // DEBUG
}

SVG::ColorObject VariableWidthInsetGenerator::getColor(edge_t& edge)
{
    switch (edge.data.type)
    {
        case SkeletalTrapezoidationEdge::TRANSITION_END:
            return SVG::Color::MAGENTA;
        case SkeletalTrapezoidationEdge::NORMAL:
        case SkeletalTrapezoidationEdge::EXTRA_VD:
        default:
            return SVG::ColorObject(100, 100, 100);
    }
}

void VariableWidthInsetGenerator::debugOutput(SVG& svg, std::unordered_map<edge_t*, std::list<TransitionMiddle>>* edge_to_transition_mids, std::unordered_map<edge_t*, std::list<TransitionEnd>>* edge_to_transition_ends)
{
    coord_t font_size = 20;
    SVG::ColorObject up_clr(255, 0, 150);
    SVG::ColorObject mid_clr(150, 0, 150);
    SVG::ColorObject down_clr(150, 0, 255);

    if (edge_to_transition_mids)
    {
        for (auto& pair : *edge_to_transition_mids)
        {
            edge_t* edge = pair.first;
            Point a = edge->from->p;
            Point b = edge->to->p;
            Point ab = b - a;
            coord_t ab_length = vSize(ab);
            for (TransitionMiddle& transition : pair.second)
            {
                Point p = a + ab * transition.pos / ab_length;
                svg.writePoint(p, false, 3, mid_clr);
                std::ostringstream ss;
                ss << transition.lower_bead_count << ">" << (transition.lower_bead_count + 1);
                svg.writeText(p, ss.str(), mid_clr, font_size);
            }
        }
    }
    if (edge_to_transition_ends)
    {
        for (auto& pair : *edge_to_transition_ends)
        {
            edge_t* edge = pair.first;
            Point a = edge->from->p;
            Point b = edge->to->p;
            Point ab = b - a;
            coord_t ab_length = vSize(ab);
            for (TransitionEnd& transition : pair.second)
            {
                Point p = a + ab * transition.pos / ab_length;
                SVG::ColorObject clr = transition.is_lower_end? down_clr : up_clr;
                svg.writePoint(p, false, 3, clr);
                coord_t end_bead_count = transition.lower_bead_count + (transition.is_lower_end? 0 : 1);
                std::ostringstream ss;
                ss << end_bead_count;
                svg.writeText(p, ss.str(), clr, font_size);
            }
        }
    }
}

void VariableWidthInsetGenerator::debugOutput(SVG& svg, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions)
{
    for (auto& pair : edge_to_junctions)
        for (ExtrusionJunction& junction : pair.second)
        {
            svg.writePoint(junction.p, false, svg.getScale() * junction.w / 2, SVG::Color::BLACK);
            
            for (double w = .95; w > .25; w = 1.0 - (1.0 - w) * 1.2)
            {
                int c = std::min(255.0, 255 - 300 * w);
                SVG::ColorObject clr(c, c, c);
                svg.writePoint(junction.p, false, svg.getScale() * junction.w / 2 * w, clr);
            }
        }
    for (auto& pair : edge_to_junctions)
        for (ExtrusionJunction& junction : pair.second)
            svg.writePoint(junction.p, false, 2, SVG::Color::YELLOW);
}

void VariableWidthInsetGenerator::debugOutput(STLwriter& stl, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions, std::unordered_map<node_t*, BeadingPropagation>& node_to_beading)
{
    auto toPoint3 = [](Point p, float h)
    {
        return Point3(p.X, p.Y, h * 200);
    };
    auto getHeight = [&node_to_beading](node_t* node)
    {
        auto found = node_to_beading.find(node);
        if (found == node_to_beading.end())
            return 0.0f;
        BeadingPropagation& beading = found->second;
        coord_t r = node->data.distance_to_boundary;
        size_t upper_inset_idx = 0;
        while (upper_inset_idx < beading.beading.toolpath_locations.size() && beading.beading.toolpath_locations[upper_inset_idx] < r)
            upper_inset_idx++;
        if (upper_inset_idx >= beading.beading.toolpath_locations.size())
        {
            return static_cast<float>(std::max(size_t(1), beading.beading.toolpath_locations.size()));
        }
        coord_t upper_inset_r = beading.beading.toolpath_locations[upper_inset_idx];
        if (upper_inset_idx <= 0)
        {
            return static_cast<float>(r) / upper_inset_r ;
        }
        coord_t lower_inset_r = beading.beading.toolpath_locations[upper_inset_idx - 1];
        return (upper_inset_idx - 1) * 2 + 1 + 2 * static_cast<float>(r - lower_inset_r) / (upper_inset_r - lower_inset_r);
    };
    auto getPoint = [](float h, std::vector<ExtrusionJunction>& junctions)
    {
        float inset_ifx_f = (h - 1) / 2;
        size_t lower_inset_idx = std::floor(inset_ifx_f);
        float rest = inset_ifx_f - lower_inset_idx;
        if (lower_inset_idx >= junctions.size()
            || lower_inset_idx + 1 >= junctions.size())
            return junctions.back().p;
        Point lower = junctions[lower_inset_idx].p;
        Point upper = junctions[lower_inset_idx + 1].p;
        coord_t scaler = 1000;
        return (lower * scaler * rest + upper * scaler * (1.0 - rest)) / scaler;
    };
    for (edge_t& edge : graph.edges)
    {
        if (edge.prev)
        {
            continue;
        }
        edge_t* quad_start = &edge;
        edge_t* quad_end = quad_start; while (quad_end->next) quad_end = quad_end->next;
        std::vector<ExtrusionJunction>* start_junctions = &edge_to_junctions[quad_start];
        std::vector<ExtrusionJunction>* end_junctions = &edge_to_junctions[quad_end->twin];
        coord_t start_junction_index = start_junctions->size() - 1;
        coord_t end_junction_index = end_junctions->size() - 1;
        
        Point3 start_prev = toPoint3(quad_start->from->p, getHeight(quad_start->from));
        Point3 end_prev = toPoint3(quad_end->to->p, getHeight(quad_end->to));
        while (true)
        {
            std::vector<ExtrusionJunction>& ss = *start_junctions;
            std::vector<ExtrusionJunction>& ee = *end_junctions;
            if (end_junction_index < 0 && start_junction_index < 0)
            {
                break;
            }
            else if (end_junction_index < 0)
            {
                if (quad_end->prev == quad_start)
                {
                    break;
                }
                Point3 end = toPoint3(quad_end->from->p, getHeight(quad_end->from));
                Point3 start_focus_point = toPoint3(quad_start->to->p, getHeight(quad_start->to));
                Point3 start = start_prev + (start_focus_point - start_prev) * (end.z - end_prev.z) / (start_focus_point.z - start_prev.z);
                stl.writeQuad(start_prev, start, end_prev, end);
                start_prev = start;
                end_prev = end;
                quad_end = quad_end->prev;
                end_junctions = &edge_to_junctions[quad_end->twin];
                if (end_junctions->empty()) break;
                end_junction_index = end_junctions->size() - 1;
            }
            else if (start_junction_index < 0)
            {
                if (quad_end == quad_start->next)
                {
                    break;
                }
                Point3 start = toPoint3(quad_start->to->p, getHeight(quad_start->to));
                Point3 end_focus_point = toPoint3(quad_end->from->p, getHeight(quad_end->from));
                Point3 end = end_prev + (end_focus_point - end_prev) * (start.z - start_prev.z) / (end_focus_point.z - end_prev.z);
                stl.writeQuad(start_prev, start, end_prev, end);
                start_prev = start;
                end_prev = end;
                quad_start = quad_start->next;
                start_junctions = &edge_to_junctions[quad_start];
                if (start_junctions->empty()) break;
                start_junction_index = start_junctions->size() - 1;
            }
            Point3 start = toPoint3((*start_junctions)[start_junction_index].p, (*start_junctions)[start_junction_index].perimeter_index * 2 + 1);
            Point3 end = toPoint3((*end_junctions)[end_junction_index].p, (*end_junctions)[end_junction_index].perimeter_index * 2 + 1);
            
            stl.writeQuad(start_prev, start, end_prev, end);
            start_prev = start;
            end_prev = end;
            start_junction_index--;
            end_junction_index--;
        }

        Point3 start = toPoint3(quad_start->to->p, getHeight(quad_start->to));
        Point3 end = toPoint3(quad_end->from->p, getHeight(quad_end->from));
        stl.writeQuad(start_prev, start, end_prev, end);
    }
}

} // namespace arachne
