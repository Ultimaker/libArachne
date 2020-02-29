//Copyright (c) 2019 Ultimaker B.V.
#include "TrapezoidationQuantizer.h"

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

TrapezoidationQuantizer::TrapezoidationQuantizer(
SkeletalTrapezoidation& st
, const BeadingStrategy& beading_strategy
, const Polygons& polys
, float transitioning_angle
, coord_t discretization_step_size
, coord_t transition_filter_dist
)
: st(st)
, beading_strategy(beading_strategy)
, polys(polys)
, transitioning_angle(transitioning_angle)
, discretization_step_size(discretization_step_size)
, transition_filter_dist(transition_filter_dist)
{
}

void TrapezoidationQuantizer::applyBeadCounts(bool filter_outermost_marked_edges)
{
    setMarking();

    filterMarking(marking_filter_dist);

    if (filter_outermost_marked_edges)
    {
        filterOuterMarking();
    }

        st.debugCheckGraphCompleteness();
        st.debugCheckGraphConsistency();

    setBeadCount();

#ifdef DEBUG
    {
        SVG svg("output/unfiltered.svg", AABB(polys));
        st.debugOutput(svg, false, false, true, false);
    }
#endif

    filterUnmarkedRegions();

    debugCheckDecorationConsistency(false);

#ifdef DEBUG
    {
        SVG svg("output/filtered.svg", AABB(polys));
        st.debugOutput(svg, false, false, true, false);
    }
#endif

    generateTransitioningRibs();

    generateExtraRibs();

#ifdef DEBUG
    {
        AABB aabb(polys);
        SVG svg("output/radial_dists.svg", aabb);
        st.debugOutput(svg, false, true);
    }
    {
        AABB aabb(polys);
        SVG svg("output/bead_counts.svg", aabb);
        st.debugOutput(svg, false, false, true);
    }
    {
        AABB aabb(polys);
        SVG svg("output/locations.svg", aabb);
        st.debugOutput(svg, false, false, false, true);
    }
#endif // DEBUG

    debugCheckDecorationConsistency(true);
}

void TrapezoidationQuantizer::setMarking()
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

    coord_t outer_edge_filter_length = beading_strategy.transitionThickness(0) / 2;

    float cap = sin(beading_strategy.transitioning_angle * 0.5); // = cos(bisector_angle / 2)
    for (edge_t& edge : st.graph.edges)
    {
        assert(edge.twin);
        if (edge.twin->markingIsSet())
        {
            edge.setMarked(edge.twin->isMarked());
        }
        else if (edge.type == SkeletalTrapezoidationEdge::EXTRA_VD)
        {
            edge.setMarked(false);
        }
        else if (std::max(edge.from->distance_to_boundary, edge.to->distance_to_boundary) < outer_edge_filter_length)
        {
            edge.setMarked(false);
        }
        else
        {
            Point a = edge.from->p;
            Point b = edge.to->p;
            Point ab = b - a;
            coord_t dR = std::abs(edge.to->distance_to_boundary - edge.from->distance_to_boundary);
            coord_t dD = vSize(ab);
            edge.setMarked(dR < dD * cap);
        }
    }
}


void TrapezoidationQuantizer::filterMarking(coord_t max_length)
{
    for (edge_t& edge : st.graph.edges)
    {
        if (isEndOfMarking(edge) && ! edge.to->isLocalMaximum() && ! edge.from->isLocalMaximum())
        {
            filterMarking(edge.twin, 0, max_length);
        }
    }
}


bool TrapezoidationQuantizer::filterMarking(edge_t* starting_edge, coord_t traveled_dist, coord_t max_length)
{
    coord_t length = vSize(starting_edge->from->p - starting_edge->to->p);
    if (traveled_dist + length > max_length)
    {
        return false;
    }
    bool should_dissolve = true;
    for (edge_t* next_edge = starting_edge->next; next_edge && next_edge != starting_edge->twin; next_edge = next_edge->twin->next)
    {
        if (next_edge->isMarked())
        {
            should_dissolve &= filterMarking(next_edge, traveled_dist + length, max_length);
        }
    }
    should_dissolve &= ! starting_edge->to->isLocalMaximum(); // don't filter marked regions with a local maximum!
    if (should_dissolve)
    {
        starting_edge->setMarked(false);
        starting_edge->twin->setMarked(false);
    }
    return should_dissolve;
}

void TrapezoidationQuantizer::filterOuterMarking()
{
    for (edge_t& edge : st.graph.edges)
    {
        if (!edge.prev)
        {
            edge.setMarked(false);
            edge.twin->setMarked(false);
        }
    }
}

void TrapezoidationQuantizer::setBeadCount()
{
    for (edge_t& edge : st.graph.edges)
    {
        if (edge.isMarked())
        {
            edge.to->bead_count = beading_strategy.optimalBeadCount(edge.to->distance_to_boundary * 2);
        }
    }

    // fix bead count at locally maximal R
    // also for marked regions!! See TODO s in generateTransitionEnd(.)
    for (node_t& node : st.graph.nodes)
    {
        if (node.isLocalMaximum())
        {
            if (node.distance_to_boundary < 0)
            {
                RUN_ONCE(logWarning("Distance to boundary not yet computed for local maximum!\n"));
                node.distance_to_boundary = std::numeric_limits<coord_t>::max();
                bool first = true;
                for (edge_t* edge = node.some_edge; first || edge != node.some_edge; edge = edge->twin->next)
                {
                    node.distance_to_boundary = std::min(node.distance_to_boundary, edge->to->distance_to_boundary + vSize(edge->from->p - edge->to->p));
                }
            }
            coord_t bead_count = beading_strategy.optimalBeadCount(node.distance_to_boundary * 2);
            node.bead_count = bead_count;
        }
    }
}

void TrapezoidationQuantizer:: filterUnmarkedRegions()
{
    for (edge_t& edge : st.graph.edges)
    {
        if (!isEndOfMarking(edge))
        {
            continue;
        }
        assert(edge.to->bead_count >= 0 || edge.to->distance_to_boundary == 0);
        coord_t max_dist = 400; // beading_strategy.getTransitioningLength(edge.to->bead_count)
        filterUnmarkedRegions(&edge, edge.to->bead_count, 0, max_dist);
    }
}

bool TrapezoidationQuantizer::filterUnmarkedRegions(edge_t* to_edge, coord_t bead_count, coord_t traveled_dist, coord_t max_dist)
{
    coord_t r = to_edge->to->distance_to_boundary;
    bool dissolve = false;
    for (edge_t* next_edge = to_edge->next; next_edge && next_edge != to_edge->twin; next_edge = next_edge->twin->next)
    {
        coord_t length = vSize(next_edge->to->p - next_edge->from->p);
        if (next_edge->to->distance_to_boundary < r && !shorterThen(next_edge->to->p - next_edge->from->p, 10))
        { // only walk upward
            continue;
        }
        if (next_edge->to->bead_count == bead_count)
        {
            dissolve = true;
        }
        else if (next_edge->to->bead_count < 0)
        {
            dissolve = filterUnmarkedRegions(next_edge, bead_count, traveled_dist + length, max_dist);
        }
        else // upward bead count is different
        {
            // dissolve if two marked regions with different bead count are closer together than the max_dist (= transition distance)
            dissolve = (traveled_dist + length < max_dist) && std::abs(next_edge->to->bead_count - bead_count) == 1;
        }
        if (dissolve)
        {
            next_edge->setMarked(true);
            next_edge->twin->setMarked(true);
            next_edge->to->bead_count = beading_strategy.optimalBeadCount(next_edge->to->distance_to_boundary * 2);
            next_edge->to->transition_ratio = 0;
        }
        return dissolve; // dissolving only depend on the one edge going upward. There cannot be multiple edges going upward.
    }
    return dissolve;
}

//
// ^^^^^^^^^^^^^^^^^^^^^
//       MARKING
// =====================
//
// =====================
//    TRANSTISIONING
// vvvvvvvvvvvvvvvvvvvvv
//

void TrapezoidationQuantizer::generateTransitioningRibs()
{
        st.debugCheckGraphCompleteness();

    generateTransitionMids();

    for (edge_t& edge : st.graph.edges)
    { // check if there is a transition in between nodes with different bead counts
        if (edge.isMarked() && edge.from->bead_count != edge.to->bead_count)
            assert(edge_to_transition_mids.find(&edge) != edge_to_transition_mids.end()
                || edge_to_transition_mids.find(edge.twin) != edge_to_transition_mids.end() );
    }
    
        st.debugCheckGraphCompleteness();
        st.debugCheckGraphConsistency();

#ifdef DEBUG
    {
        SVG svg("output/transition_mids_unfiltered.svg", AABB(polys));
        st.debugOutput(svg, false, false, true, false);
        debugOutput(svg);
    }
#endif

    filterTransitionMids(edge_to_transition_mids);

#ifdef DEBUG
    {
        SVG svg("output/transition_mids.svg", AABB(polys));
        st.debugOutput(svg, false, false, true, false);
        debugOutput(svg);
    }
#endif

    debugCheckTransitionMids();

    generateTransitionEnds();

#ifdef DEBUG
    {
        SVG svg("output/transition_ends.svg", AABB(polys));
        st.debugOutput(svg, false, false, true, false);
        debugOutput(svg);
    }
#endif

    applyTransitions(edge_to_transition_ends);
}


void TrapezoidationQuantizer::generateTransitionMids()
{
    for (edge_t& edge : st.graph.edges)
    {
        assert(edge.markingIsSet());
        if (!edge.isMarked())
        { // only marked regions introduce transitions
            continue;
        }
        coord_t start_R = edge.from->distance_to_boundary;
        coord_t end_R = edge.to->distance_to_boundary;
        coord_t start_bead_count = edge.from->bead_count;
        coord_t end_bead_count = edge.to->bead_count;

        if (start_R == end_R)
        { // no transitions occur when both end points have the same distance_to_boundary
            assert(edge.from->bead_count == edge.to->bead_count);// TODO: what to do in this case?
            continue;
        }
        else if (start_R > end_R)
        { // only consider those half-edges which are going from a lower to a higher distance_to_boundary
            continue;
        }

        if (edge.from->bead_count == edge.to->bead_count)
        { // no transitions should accur according to the enforced bead counts
            continue;
        }

        if (start_bead_count > beading_strategy.optimalBeadCount(start_R * 2)
            || end_bead_count > beading_strategy.optimalBeadCount(end_R * 2))
        { // wasn't the case earlier in this function because of already introduced transitions
            RUN_ONCE(logError("transitioning segment overlap! (?)\n"));
        }
        assert(start_R < end_R);
        coord_t edge_size = vSize(edge.from->p - edge.to->p);
        for (coord_t transition_lower_bead_count = start_bead_count; transition_lower_bead_count < end_bead_count; transition_lower_bead_count++)
        {
            coord_t mid_R = beading_strategy.transitionThickness(transition_lower_bead_count) / 2;
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
            assert(edge_to_transition_mids[&edge].empty() || mid_pos >= edge_to_transition_mids[&edge].back().pos);
            edge_to_transition_mids[&edge].emplace_back(mid_pos, transition_lower_bead_count);
        }
        if (edge.from->bead_count != edge.to->bead_count)
        {
            assert(edge_to_transition_mids[&edge].size() >= 1);
        }
    }
}

void TrapezoidationQuantizer::filterTransitionMids(std::unordered_map<edge_t*, std::list<TransitionMiddle>>& edge_to_transition_mids)
{
    for (auto pair_it = edge_to_transition_mids.begin(); pair_it != edge_to_transition_mids.end();)
    {
        std::pair<edge_t* const, std::list<TransitionMiddle>>& pair = *pair_it;
        edge_t* edge = pair.first;
        std::list<TransitionMiddle>& transitions = pair.second;
        if (transitions.empty())
        {
            pair_it = edge_to_transition_mids.erase(pair_it);
            continue;
        }
        assert(transitions.front().lower_bead_count <= transitions.back().lower_bead_count); // this is how stuff should be stored in edge_to_transition_mids
        assert(edge->from->distance_to_boundary <= edge->to->distance_to_boundary); // this is how stuff should be stored in edge_to_transition_mids
        Point a = edge->from->p;
        Point b = edge->to->p;
        Point ab = b - a;
        coord_t ab_size = vSize(ab);

        bool going_up = true;
        std::list<TransitionMidRef> to_be_dissolved_back = dissolveNearbyTransitions(edge, transitions.back(), ab_size - transitions.back().pos, transition_filter_dist, going_up);
        bool should_dissolve_back = !to_be_dissolved_back.empty();
        for (TransitionMidRef& ref : to_be_dissolved_back)
        {
            dissolveBeadCountRegion(edge, transitions.back().lower_bead_count + 1, transitions.back().lower_bead_count);
            if (ref.pair_it->second.size() <= 1)
            {
                edge_to_transition_mids.erase(ref.pair_it);
            }
            else
            {
                ref.pair_it->second.erase(ref.transition_it);
            }
        }

        {
            coord_t trans_bead_count = transitions.back().lower_bead_count;
            coord_t upper_transition_half_length = (1.0 - beading_strategy.getTransitionAnchorPos(trans_bead_count)) * beading_strategy.getTransitioningLength(trans_bead_count);
            should_dissolve_back |= filterEndOfMarkingTransition(edge, ab_size - transitions.back().pos, upper_transition_half_length, trans_bead_count);
        }
        if (should_dissolve_back)
        {
            transitions.pop_back();
        }
        if (transitions.empty())
        { // filterEndOfMarkingTransition gives inconsistent new bead count when executing for the same transition in two directions.
            pair_it = edge_to_transition_mids.erase(pair_it);
            continue;
        }

        going_up = false;
        std::list<TransitionMidRef> to_be_dissolved_front = dissolveNearbyTransitions(edge->twin, transitions.front(), transitions.front().pos, transition_filter_dist, going_up);
        bool should_dissolve_front = !to_be_dissolved_front.empty();
        for (TransitionMidRef& ref : to_be_dissolved_front)
        {
            dissolveBeadCountRegion(edge->twin, transitions.front().lower_bead_count, transitions.front().lower_bead_count + 1);
            if (ref.pair_it->second.size() <= 1)
            {
                edge_to_transition_mids.erase(ref.pair_it);
            }
            else
            {
                ref.pair_it->second.erase(ref.transition_it);
            }
        }

        {
            coord_t trans_bead_count = transitions.front().lower_bead_count;
            coord_t lower_transition_half_length = beading_strategy.getTransitionAnchorPos(trans_bead_count) * beading_strategy.getTransitioningLength(trans_bead_count);
            should_dissolve_front |= filterEndOfMarkingTransition(edge->twin, transitions.front().pos, lower_transition_half_length, trans_bead_count + 1);
        }
        if (should_dissolve_front)
        {
            transitions.pop_front();
        }
        if (transitions.empty())
        { // filterEndOfMarkingTransition gives inconsistent new bead count when executing for the same transition in two directions.
            pair_it = edge_to_transition_mids.erase(pair_it);
            continue;
        }
        ++pair_it; // normal update of loop
    }
}

std::list<TrapezoidationQuantizer::TransitionMidRef> TrapezoidationQuantizer::dissolveNearbyTransitions(edge_t* edge_to_start, TransitionMiddle& origin_transition, coord_t traveled_dist, coord_t max_dist, bool going_up)
{
    std::list<TransitionMidRef> to_be_dissolved;
    if (traveled_dist > max_dist)
    {
        return to_be_dissolved;
    }
    bool should_dissolve = true;
    for (edge_t* edge = edge_to_start->next; edge && edge != edge_to_start->twin; edge = edge->twin->next)
    {
        if (!edge->isMarked())
        {
            continue;
        }
        Point a = edge->from->p;
        Point b = edge->to->p;
        Point ab = b - a;
        coord_t ab_size = vSize(ab);
        bool is_aligned = edge->isUpward();
        edge_t* aligned_edge = is_aligned? edge : edge->twin;
        bool seen_transition_on_this_edge = false;
        auto edge_transitions_it = edge_to_transition_mids.find(aligned_edge);
        if (edge_transitions_it != edge_to_transition_mids.end())
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
            std::list<TrapezoidationQuantizer::TransitionMidRef> to_be_dissolved_here = dissolveNearbyTransitions(edge, origin_transition, traveled_dist + ab_size, max_dist, going_up);
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


void TrapezoidationQuantizer::dissolveBeadCountRegion(edge_t* edge_to_start, coord_t from_bead_count, coord_t to_bead_count)
{
    assert(from_bead_count != to_bead_count);
    if (edge_to_start->to->bead_count != from_bead_count)
    {
        return;
    }
    edge_to_start->to->bead_count = to_bead_count;
    for (edge_t* edge = edge_to_start->next; edge && edge != edge_to_start->twin; edge = edge->twin->next)
    {
        if (!edge->isMarked())
        {
            continue;
        }
        dissolveBeadCountRegion(edge, from_bead_count, to_bead_count);
    }
}

bool TrapezoidationQuantizer::filterEndOfMarkingTransition(edge_t* edge_to_start, coord_t traveled_dist, coord_t max_dist, coord_t replacing_bead_count)
{
    if (traveled_dist > max_dist)
    {
        return false;
    }
    bool is_end_of_marking = true;
    bool should_dissolve = false;
    for (edge_t* next_edge = edge_to_start->next; next_edge && next_edge != edge_to_start->twin; next_edge = next_edge->twin->next)
    {
        if (next_edge->isMarked())
        {
            coord_t length = vSize(next_edge->to->p - next_edge->from->p);
            should_dissolve |= filterEndOfMarkingTransition(next_edge, traveled_dist + length, max_dist, replacing_bead_count);
            is_end_of_marking = false;
        }
    }
    if (is_end_of_marking && traveled_dist < max_dist)
    {
        should_dissolve = true;
    }
    if (should_dissolve)
    {
        edge_to_start->to->bead_count = replacing_bead_count;
    }
    return should_dissolve;
}

void TrapezoidationQuantizer::generateTransitionEnds()
{
    for (std::pair<edge_t*, std::list<TransitionMiddle>> pair : edge_to_transition_mids)
    {
        edge_t* edge = pair.first;
        std::list<TransitionMiddle>& transition_positions = pair.second;

        assert(edge->from->distance_to_boundary <= edge->to->distance_to_boundary);
        for (TransitionMiddle& transition_middle : transition_positions)
        {
            assert(transition_positions.front().pos <= transition_middle.pos);
            assert(transition_middle.pos <= transition_positions.back().pos);
            generateTransition(*edge, transition_middle.pos, beading_strategy, transition_middle.lower_bead_count);
        }
    }
}

void TrapezoidationQuantizer::generateTransition(edge_t& edge, coord_t mid_pos, const BeadingStrategy& beading_strategy, coord_t lower_bead_count)
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

    
        st.debugCheckGraphCompleteness();
        st.debugCheckGraphConsistency();

    { // lower bead count transition end
        coord_t start_pos = ab_size - mid_pos;
        coord_t transition_half_length = transition_mid_position * transition_length;
        coord_t end_pos = start_pos + transition_half_length;
        generateTransitionEnd(*edge.twin, start_pos, end_pos, transition_half_length, mid_rest, start_rest, lower_bead_count);
    }
    st.debugCheckGraphConsistency();
    { // upper bead count transition end
        coord_t start_pos = mid_pos;
        coord_t transition_half_length = (1.0 - transition_mid_position) * transition_length;
        coord_t end_pos = mid_pos +  transition_half_length;
        bool is_going_down_everywhere = generateTransitionEnd(edge, start_pos, end_pos, transition_half_length, mid_rest, end_rest, lower_bead_count);
        assert(!is_going_down_everywhere && "There must have been at least one direction in which the bead count is increasing enough for the transition to happen!");
    }

        st.debugCheckGraphCompleteness();
        st.debugCheckGraphConsistency();
}

bool TrapezoidationQuantizer::generateTransitionEnd(edge_t& edge, coord_t start_pos, coord_t end_pos, coord_t transition_half_length, float start_rest, float end_rest, coord_t lower_bead_count)
{
    Point a = edge.from->p;
    Point b = edge.to->p;
    Point ab = b - a;
    coord_t ab_size = vSize(ab); // TODO: prevent recalculation of these values

    assert(start_pos <= ab_size);

    bool going_up = end_rest > start_rest;

    assert(edge.isMarked());
    if (!edge.isMarked())
    { // This function shouldn't generate ends in or beyond unmarked regions
        return false;
    }

    if (end_pos > ab_size)
    { // recurse on all further edges
        coord_t R = edge.to->distance_to_boundary;
        float rest = end_rest - (start_rest - end_rest) * (end_pos - ab_size) / (start_pos - end_pos);
        assert(rest >= 0);
        assert(rest <= std::max(end_rest, start_rest));
        assert(rest >= std::min(end_rest, start_rest));

        coord_t marked_edge_count = 0;
        for (edge_t* outgoing = edge.next; outgoing && outgoing != edge.twin; outgoing = outgoing->twin->next)
        {
            if (!outgoing->isMarked()) continue;
            marked_edge_count++;
        }

        bool is_only_going_down = true;
        bool has_recursed = false;
        for (edge_t* outgoing = edge.next; outgoing && outgoing != edge.twin;)
        {
            edge_t* next = outgoing->twin->next; // before we change the outgoing edge itself
            if (!outgoing->isMarked())
            {
                outgoing = next;
                continue; // don't put transition ends in non-marked regions
            }
            if (marked_edge_count > 1 && going_up && isGoingDown(outgoing, 0, end_pos - ab_size + transition_half_length, lower_bead_count))
            { // we're after a 3-way_all-marked_junction-node and going in the direction of lower bead count
                // don't introduce a transition end along this marked direction, because this direction is the downward direction
                // while we are supposed to be [going_up]
                outgoing = next;
                continue;
            }
            bool is_going_down = generateTransitionEnd(*outgoing, 0, end_pos - ab_size, transition_half_length, rest, end_rest, lower_bead_count);
            is_only_going_down &= is_going_down;
            outgoing = next;
            has_recursed = true;
        }
        if (!going_up || (has_recursed && !is_only_going_down))
        {
            edge.to->transition_ratio = rest;
            edge.to->bead_count = lower_bead_count;
        }
        return is_only_going_down;
    }
    else // end_pos < ab_size
    { // add transition end point here
//         assert(edge.isMarked() && "we should only be adding transition ends in marked regions");
        
        bool is_lower_end = end_rest == 0; // TODO collapse this parameter into the bool for which it is used here!
        std::list<TransitionEnd>* transitions = nullptr;
        coord_t pos = -1;
        if (edge.isUpward())
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


bool TrapezoidationQuantizer::isGoingDown(edge_t* outgoing, coord_t traveled_dist, coord_t max_dist, coord_t lower_bead_count) const
{
    // NOTE: the logic below is not fully thought through.
    // TODO: take transition mids into account
    if (outgoing->to->distance_to_boundary == 0)
    {
        return true;
    }
    bool is_upward = outgoing->to->distance_to_boundary >= outgoing->from->distance_to_boundary;
    edge_t* upward_edge = is_upward? outgoing : outgoing->twin;
    if (outgoing->to->bead_count > lower_bead_count + 1)
    {
        assert(edge_to_transition_mids.find(upward_edge) != edge_to_transition_mids.end() && "If the bead count is going down there has to be a transition mid!");
        return false;
    }
    coord_t length = vSize(outgoing->to->p - outgoing->from->p);
    auto transition_mids_it = edge_to_transition_mids.find(upward_edge);
    if (transition_mids_it != edge_to_transition_mids.end())
    {
        const std::list<TransitionMiddle>& transition_mids = transition_mids_it->second;
        if (!transition_mids.empty())
        {
            const TransitionMiddle& mid = is_upward? transition_mids.front() : transition_mids.back();
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
    if (outgoing->to->bead_count <= lower_bead_count
        && !(outgoing->to->bead_count == lower_bead_count && outgoing->to->transition_ratio > 0.0))
    {
        return true;
    }
    
    bool is_only_going_down = true;
    bool has_recursed = false;
    for (edge_t* next = outgoing->next; next && next != outgoing->twin; next = next->twin->next)
    {
        if (!next->isMarked())
        {
            continue;
        }
        bool is_going_down = isGoingDown(next, traveled_dist + length, max_dist, lower_bead_count);
        is_only_going_down &= is_going_down;
        has_recursed = true;
    }
    return has_recursed && is_only_going_down;
}

void TrapezoidationQuantizer::applyTransitions(std::unordered_map<edge_t*, std::list<TransitionEnd>>& edge_to_transition_ends)
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
        assert(edge->isMarked());

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
                && close_node->bead_count == new_node_bead_count
            )
            {
                assert(end_pos <= ab_size);
//                 close_node->bead_count = new_node_bead_count;
                close_node->transition_ratio = 0;
                continue;
            }
            Point mid = a + normal(ab, end_pos);
            
            st.debugCheckGraphCompleteness();
            st.debugCheckGraphConsistency();

            debugCheckDecorationConsistency(false);

            assert(last_edge_replacing_input->isMarked());
            assert(last_edge_replacing_input->type != SkeletalTrapezoidationEdge::EXTRA_VD);
            last_edge_replacing_input = insertNode(last_edge_replacing_input, mid, new_node_bead_count);
            assert(last_edge_replacing_input->type != SkeletalTrapezoidationEdge::EXTRA_VD);
            assert(last_edge_replacing_input->isMarked());


            st.debugCheckGraphCompleteness();
            st.debugCheckGraphConsistency();
        }
    }
}

TrapezoidationQuantizer::edge_t* TrapezoidationQuantizer::insertNode(edge_t* edge, Point mid, coord_t mide_node_bead_count)
{
    edge_t* last_edge_replacing_input = edge;

    st.graph.nodes.emplace_back(mid);
    node_t* mid_node = &st.graph.nodes.back();

    edge_t* twin = last_edge_replacing_input->twin;
    last_edge_replacing_input->twin = nullptr;
    twin->twin = nullptr;
    std::pair<TrapezoidationQuantizer::edge_t*, TrapezoidationQuantizer::edge_t*> left_pair
        = insertRib(*last_edge_replacing_input, mid_node);
    std::pair<TrapezoidationQuantizer::edge_t*, TrapezoidationQuantizer::edge_t*> right_pair
        = insertRib(*twin, mid_node);
    edge_t* first_edge_replacing_input = left_pair.first;
    last_edge_replacing_input = left_pair.second;
    edge_t* first_edge_replacing_twin = right_pair.first;
    edge_t* last_edge_replacing_twin = right_pair.second;

    first_edge_replacing_input->twin = last_edge_replacing_twin;
    last_edge_replacing_twin->twin = first_edge_replacing_input;
    last_edge_replacing_input->twin = first_edge_replacing_twin;
    first_edge_replacing_twin->twin = last_edge_replacing_input;

    mid_node->bead_count = mide_node_bead_count;

    return last_edge_replacing_input;
}

std::pair<TrapezoidationQuantizer::edge_t*, TrapezoidationQuantizer::edge_t*> TrapezoidationQuantizer::insertRib(edge_t& edge, node_t* mid_node)
{
    st.debugCheckGraphConsistency();
    edge_t* edge_before = edge.prev;
    edge_t* edge_after = edge.next;
    node_t* node_before = edge.from;
    node_t* node_after = edge.to;
    
    Point p = mid_node->p;

    std::pair<Point, Point> source_segment = getSource(edge);
    Point px = LinearAlg2D::getClosestOnLineSegment(p, source_segment.first, source_segment.second);
    coord_t dist = vSize(p - px);
    assert(dist > 0);
    mid_node->distance_to_boundary = dist;
    mid_node->transition_ratio = 0; // both transition end should have rest = 0, because at the ends a whole number of beads fits without rest

    st.graph.nodes.emplace_back(px);
    node_t* source_node = &st.graph.nodes.back();
    source_node->distance_to_boundary = 0;

    edge_t* first = &edge;
    st.graph.edges.emplace_back(SkeletalTrapezoidationEdge());
    edge_t* second = &st.graph.edges.back();
    st.graph.edges.emplace_back(SkeletalTrapezoidationEdge(SkeletalTrapezoidationEdge::TRANSITION_END));
    edge_t* outward_edge = &st.graph.edges.back();
    st.graph.edges.emplace_back(SkeletalTrapezoidationEdge(SkeletalTrapezoidationEdge::TRANSITION_END));
    edge_t* inward_edge = &st.graph.edges.back();

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

    first->setMarked(true);
    outward_edge->setMarked(false); // TODO verify this is always the case.
    inward_edge->setMarked(false);
    second->setMarked(true);

    outward_edge->twin = inward_edge;
    inward_edge->twin = outward_edge;

    first->twin = nullptr; // we don't know these yet!
    second->twin = nullptr;

    assert(second->prev->from->distance_to_boundary == 0);

    st.debugCheckGraphConsistency();

    return std::make_pair(first, second);
}
std::pair<Point, Point> TrapezoidationQuantizer::getSource(const edge_t& edge)
{
    const edge_t* from_edge;
    for (from_edge = &edge; from_edge->prev; from_edge = from_edge->prev) {}
    const edge_t* to_edge;
    for (to_edge = &edge; to_edge->next; to_edge = to_edge->next) {}
    return std::make_pair(from_edge->from->p, to_edge->to->p);
}

bool TrapezoidationQuantizer::isEndOfMarking(const edge_t& edge_to) const
{
    if (!edge_to.isMarked())
    {
        return false;
    }
    if (!edge_to.next)
    {
        return true;
    }
    for (const edge_t* edge = edge_to.next; edge && edge != edge_to.twin; edge = edge->twin->next)
    {
        if (edge->isMarked())
        {
            return false;
        }
        assert(edge->twin);
    }
    return true;
}




void TrapezoidationQuantizer::generateExtraRibs()
{
    auto end_edge_it = --st.graph.edges.end(); // don't check newly introduced edges
    for (auto edge_it = st.graph.edges.begin(); std::prev(edge_it) != end_edge_it; ++edge_it)
    {
        edge_t& edge = *edge_it;
        if ( ! edge.isMarked()) continue;
        if (shorterThen(edge.to->p - edge.from->p, discretization_step_size)) continue;
        if (edge.from->distance_to_boundary >= edge.to->distance_to_boundary) continue;

        std::vector<coord_t> rib_thicknesses = beading_strategy.getNonlinearThicknesses(edge.from->bead_count);

        if (rib_thicknesses.empty()) continue;

        // preload some variables before [edge] gets changed
        node_t* from = edge.from;
        node_t* to = edge.to;
        Point a = from->p;
        Point b = to->p;
        Point ab = b - a;
        coord_t ab_size = vSize(ab);
        coord_t a_R = edge.from->distance_to_boundary;
        coord_t b_R = edge.to->distance_to_boundary;
        
        edge_t* last_edge_replacing_input = &edge;
        for (coord_t rib_thickness : rib_thicknesses)
        {
            if (rib_thickness / 2 <= a_R) continue;
            if (rib_thickness / 2 >= b_R) break;
            coord_t new_node_bead_count = std::min(edge.from->bead_count, edge.to->bead_count);
            coord_t end_pos = ab_size * (rib_thickness / 2 - a_R) / (b_R - a_R);
            assert(end_pos > 0);
            assert(end_pos < ab_size);
            node_t* close_node = (end_pos < ab_size / 2)? from : to;
            if ((end_pos < snap_dist || end_pos > ab_size - snap_dist)
                && close_node->bead_count == new_node_bead_count
            )
            {
                assert(end_pos <= ab_size);
//                 close_node->bead_count = new_node_bead_count;
                close_node->transition_ratio = 0;
                continue;
            }
            Point mid = a + normal(ab, end_pos);
            
            st.debugCheckGraphCompleteness();
            st.debugCheckGraphConsistency();

            debugCheckDecorationConsistency(false);

            assert(last_edge_replacing_input->isMarked());
            assert(last_edge_replacing_input->type != SkeletalTrapezoidationEdge::EXTRA_VD);
            last_edge_replacing_input = insertNode(last_edge_replacing_input, mid, new_node_bead_count);
            assert(last_edge_replacing_input->type != SkeletalTrapezoidationEdge::EXTRA_VD);
            assert(last_edge_replacing_input->isMarked());


            st.debugCheckGraphCompleteness();
            st.debugCheckGraphConsistency();
        }
    }
}

//
// ^^^^^^^^^^^^^^^^^^^^^
//    TRANSTISIONING
// =====================
//
// =====================
//       HELPERS
// vvvvvvvvvvvvvvvvvvvvv
//

void TrapezoidationQuantizer::debugCheckDecorationConsistency(bool transitioned)
{
#ifdef DEBUG
    for (const edge_t& edge : st.graph.edges)
    {
        const edge_t* edge_p = &edge;
        assert(edge.type >= SkeletalTrapezoidationEdge::NORMAL);
        if (edge.type != SkeletalTrapezoidationEdge::NORMAL)
        {
            if (edge.from->distance_to_boundary != -1 && edge.to->distance_to_boundary != -1)
            {
                assert(edge.from->distance_to_boundary == 0 || edge.to->distance_to_boundary == 0);
            }
            assert(!edge.isMarked());
        }
        assert(edge.isMarked() == edge.twin->isMarked());
        if (edge.isMarked())
        {
            if (transitioned && edge.from->bead_count != -1 && edge.to->bead_count != -1)
            {
                assert(!edge.isMarked() || std::abs(edge.from->bead_count - edge.to->bead_count) <= 1);
            }
        }
    }
#endif // DEBUG
}

void TrapezoidationQuantizer::debugCheckTransitionMids() const
{
#ifdef DEBUG
    for (std::pair<edge_t*, std::list<TransitionMiddle>> pair : edge_to_transition_mids)
    {
        const edge_t* edge = pair.first;
        const std::list<TransitionMiddle>& transition_positions = pair.second;
        
        assert(edge->from->distance_to_boundary <= edge->to->distance_to_boundary);
        
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

SVG::ColorObject TrapezoidationQuantizer::getColor(edge_t& edge)
{
    switch (edge.type)
    {
        case SkeletalTrapezoidationEdge::TRANSITION_END:
            return SVG::Color::MAGENTA;
        case SkeletalTrapezoidationEdge::NORMAL:
        case SkeletalTrapezoidationEdge::EXTRA_VD:
        default:
            return SVG::ColorObject(100, 100, 100);
    }
}

void TrapezoidationQuantizer::debugOutput(SVG& svg)
{
    coord_t font_size = 20;
    SVG::ColorObject up_clr(255, 0, 150);
    SVG::ColorObject mid_clr(150, 0, 150);
    SVG::ColorObject down_clr(150, 0, 255);

    if ( ! edge_to_transition_mids.empty())
    {
        for (auto& pair : edge_to_transition_mids)
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
    if ( ! edge_to_transition_ends.empty())
    {
        for (auto& pair : edge_to_transition_ends)
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

} // namespace arachne
