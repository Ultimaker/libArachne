//Copyright (c) 2019 Ultimaker B.V.
#include "VariableWidthInsetGenerator.h"

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

VariableWidthInsetGenerator::VariableWidthInsetGenerator(const Polygons& polys, float transitioning_angle
, coord_t discretization_step_size
, coord_t transition_filter_dist
, coord_t beading_propagation_transition_dist
)
: st(polys, transitioning_angle, discretization_step_size)
, polys(polys)
, transitioning_angle(transitioning_angle)
, discretization_step_size()
, transition_filter_dist(transition_filter_dist)
, beading_propagation_transition_dist(beading_propagation_transition_dist)
{
}

std::vector<std::list<ExtrusionLine>> VariableWidthInsetGenerator::generateToolpaths(const BeadingStrategy& beading_strategy, bool filter_outermost_marked_edges)
{
    st.initialize_graph(); // generate st.graph
    
    TrapezoidationQuantizer quantizer(st, beading_strategy, polys, transitioning_angle, discretization_step_size, transition_filter_dist);
    quantizer.applyBeadCounts(filter_outermost_marked_edges);

    std::vector<std::list<ExtrusionLine>> result_polylines_per_index;
    generateSegments(result_polylines_per_index, beading_strategy);
    // junctions = generateJunctions

    return result_polylines_per_index;
}


void VariableWidthInsetGenerator::generateSegments(std::vector<std::list<ExtrusionLine>>& result_polylines_per_index, const BeadingStrategy& beading_strategy)
{
    std::vector<edge_t*> upward_quad_mids;
    for (edge_t& edge : st.graph.edges)
    {
        if (edge.prev && edge.next && st.isUpward(&edge))
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
                    coord_t a_dist_from_up = std::min(st.distToGoUp(a).value_or(max), st.distToGoUp(a->twin).value_or(max)) - vSize(a->to->p - a->from->p);
                    coord_t b_dist_from_up = std::min(st.distToGoUp(b).value_or(max), st.distToGoUp(b->twin).value_or(max)) - vSize(b->to->p - b->from->p);
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
        for (node_t& node : st.graph.nodes)
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
        st.debugOutput(svg, false, false, true, false);
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


VariableWidthInsetGenerator::Beading VariableWidthInsetGenerator::interpolate(const Beading& top, float ratio_top_to_whole, const Beading& bottom, coord_t switching_radius) const
{
    assert( ratio_top_to_whole >= 0.0 && ratio_top_to_whole <= 1.0);
    Beading ret = interpolate(top, ratio_top_to_whole, bottom);

    // TODO: don't use toolpath locations past the middle!
    // TODO: stretch bead widths and locations of the higher bead count beading to fit in the left over space
    coord_t next_inset_idx;
    for (next_inset_idx = top.toolpath_locations.size() - 1; next_inset_idx >= 0; next_inset_idx--)
    {
        if (switching_radius > top.toolpath_locations[next_inset_idx])
        {
            break;
        }
    }
    if (next_inset_idx < 0)
    { // there is no next inset, because there is only one
        assert(top.toolpath_locations.empty() || top.toolpath_locations.front() >= switching_radius);
        return ret;
    }
    if (next_inset_idx + 1 == coord_t( top.toolpath_locations.size()))
    { // we cant adjust to fit the next edge because there is no previous one?!
        return ret;
    }
    assert(next_inset_idx < coord_t( top.toolpath_locations.size()));
    assert(top.toolpath_locations[next_inset_idx] <= switching_radius);
    assert(top.toolpath_locations[next_inset_idx + 1] >= switching_radius);
    if (ret.toolpath_locations[next_inset_idx] > switching_radius)
    { // one inset disappeared between left and the merged one
        // solve for ratio f:
        // f*l + (1-f)*r = s
        // f*l + r - f*r = s
        // f*(l-r) + r = s
        // f*(l-r) = s - r
        // f = (s-r) / (l-r)
        float new_ratio = static_cast<float>(switching_radius - bottom.toolpath_locations[next_inset_idx]) / static_cast<float>( top.toolpath_locations[next_inset_idx] - bottom.toolpath_locations[next_inset_idx]);
        new_ratio = std::min(1.0, new_ratio + 0.1);
        return interpolate( top, new_ratio, bottom);
    }
    return ret;
}


VariableWidthInsetGenerator::Beading VariableWidthInsetGenerator::interpolate(const Beading& top, float ratio_top_to_whole, const Beading& bottom) const
{
    assert( ratio_top_to_whole >= 0.0 && ratio_top_to_whole <= 1.0);
    float ratio_right_to_whole = 1.0 - ratio_top_to_whole;

    Beading ret = (top.total_thickness > bottom.total_thickness)? top : bottom;
    for (size_t inset_idx = 0; inset_idx < std::min( top.bead_widths.size(), bottom.bead_widths.size()); inset_idx++)
    {
        ret.bead_widths[inset_idx] = ratio_top_to_whole * top.bead_widths[inset_idx] + ratio_right_to_whole * bottom.bead_widths[inset_idx];
        ret.toolpath_locations[inset_idx] = ratio_top_to_whole * top.toolpath_locations[inset_idx] + ratio_right_to_whole * bottom.toolpath_locations[inset_idx];
    }
    return ret;
}

void VariableWidthInsetGenerator::generateJunctions(std::unordered_map<node_t*, BeadingPropagation>& node_to_beading, std::unordered_map<edge_t*, std::vector<ExtrusionJunction>>& edge_to_junctions, const BeadingStrategy& beading_strategy)
{
    for (edge_t& edge_ : st.graph.edges)
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
    
    std::unordered_set<edge_t*> unprocessed_quad_starts(st.graph.edges.size() * 5 / 2);
    for (edge_t& edge : st.graph.edges)
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
            && st.isLocalMaximum(*node, true)
            && !st.isMarked(node)
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
    for (const edge_t& edge : st.graph.edges)
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
    for (edge_t& edge : st.graph.edges)
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
