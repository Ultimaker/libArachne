//Copyright (c) 2019 Ultimaker B.V.
#include "SkeletalTrapezoidation.h"

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

bool generate_MAT_STL = false;

SkeletalTrapezoidation::node_t& SkeletalTrapezoidation::make_node(vd_t::vertex_type& vd_node, Point p)
{
    auto he_node_it = vd_node_to_he_node.find(&vd_node);
    if (he_node_it == vd_node_to_he_node.end())
    {
        graph.nodes.emplace_front(SkeletalTrapezoidationJoint(), p);
        node_t& node = graph.nodes.front();
        vd_node_to_he_node.emplace(&vd_node, &node);
        return node;
    }
    else
    {
        return *he_node_it->second;
    }
}

void SkeletalTrapezoidation::transfer_edge(Point from, Point to, vd_t::edge_type& vd_edge, edge_t*& prev_edge, Point& start_source_point, Point& end_source_point, const std::vector<Point>& points, const std::vector<Segment>& segments)
{
    auto he_edge_it = vd_edge_to_he_edge.find(vd_edge.twin());
    if (he_edge_it != vd_edge_to_he_edge.end())
    { // twin segment(s) already made
        edge_t* source_twin = he_edge_it->second;
        assert(source_twin);
        auto end_node_it = vd_node_to_he_node.find(vd_edge.vertex1());
        assert(end_node_it != vd_node_to_he_node.end());
        node_t* end_node = end_node_it->second;
        for (edge_t* twin = source_twin;
            ; //  !twin;
            twin = twin->prev->twin->prev)
        {
            assert(twin);
            graph.edges.emplace_front(SkeletalTrapezoidationEdge());
            edge_t* edge = &graph.edges.front();
            edge->from = twin->to;
            edge->to = twin->from;
            edge->twin = twin;
            twin->twin = edge;
            edge->from->some_edge = edge;
            
            if (prev_edge)
            {
                edge->prev = prev_edge;
                prev_edge->next = edge;
            }

            prev_edge = edge;

            if (prev_edge->to == end_node)
            {
                return;
            }
            
            if (!twin->prev || !twin->prev->twin || !twin->prev->twin->prev)
            {
                RUN_ONCE(logError("Discretized segment behaves oddly!\n"));
                return;
            }
            assert(twin->prev); // forth rib
            assert(twin->prev->twin); // back rib
            assert(twin->prev->twin->prev); // prev segment along parabola
            make_rib(prev_edge, start_source_point, end_source_point);
        }
        assert(prev_edge);
    }
    else
    {
        std::vector<Point> discretized = discretize(vd_edge, points, segments);
        assert(discretized.size() >= 2);
        
        assert(!prev_edge || prev_edge->to);
//         assert(!prev_edge || prev_edge->to == &make_node(*vd_edge.vertex0(), from)); // TODO: investigate whether boost:voronoi can produce multiple verts and violates consistency
        node_t* v0 = (prev_edge)? prev_edge->to : &make_node(*vd_edge.vertex0(), from);
        Point p0 = discretized.front();
        for (size_t p1_idx = 1; p1_idx < discretized.size(); p1_idx++)
        {
            Point p1 = discretized[p1_idx];
            node_t* v1;
            if (p1_idx < discretized.size() - 1)
            {
                graph.nodes.emplace_front(SkeletalTrapezoidationJoint(), p1);
                v1 = &graph.nodes.front();
            }
            else
            {
                v1 = &make_node(*vd_edge.vertex1(), to);
            }

            graph.edges.emplace_front(SkeletalTrapezoidationEdge());
            edge_t* edge = &graph.edges.front();
            edge->from = v0;
            edge->to = v1;
//             edge->twin = nullptr;
            edge->from->some_edge = edge;
            
            if (prev_edge)
            {
                edge->prev = prev_edge;
                prev_edge->next = edge;
            }
            
            prev_edge = edge;
            p0 = p1;
            v0 = v1;
            
            if (p1_idx < discretized.size() - 1)
            { // rib for last segment gets introduced outside this function!
                make_rib(prev_edge, start_source_point, end_source_point);
            }
        }
        assert(prev_edge);
        vd_edge_to_he_edge.emplace(&vd_edge, prev_edge);
    }
}

void SkeletalTrapezoidation::make_rib(edge_t*& prev_edge, Point start_source_point, Point end_source_point)
{
    Point p = LinearAlg2D::getClosestOnLineSegment(prev_edge->to->p, start_source_point, end_source_point);
    coord_t dist = vSize(prev_edge->to->p - p);
    prev_edge->to->data.distance_to_boundary = dist;
    assert(dist >= 0);

    graph.nodes.emplace_front(SkeletalTrapezoidationJoint(), p);
    node_t* node = &graph.nodes.front();
    node->data.distance_to_boundary = 0;
    
    graph.edges.emplace_front(SkeletalTrapezoidationEdge(SkeletalTrapezoidationEdge::EXTRA_VD));
    edge_t* forth_edge = &graph.edges.front();
    graph.edges.emplace_front(SkeletalTrapezoidationEdge(SkeletalTrapezoidationEdge::EXTRA_VD));
    edge_t* back_edge = &graph.edges.front();
    
    prev_edge->next = forth_edge;
    forth_edge->prev = prev_edge;
    forth_edge->from = prev_edge->to;
    forth_edge->to = node;
    forth_edge->twin = back_edge;
    back_edge->twin = forth_edge;
    back_edge->from = node;
    back_edge->to = prev_edge->to;
    node->some_edge = back_edge;
    
    prev_edge = back_edge;
}

std::vector<Point> SkeletalTrapezoidation::discretize(const vd_t::edge_type& vd_edge, const std::vector<Point>& points, const std::vector<Segment>& segments)
{
    const vd_t::cell_type* left_cell = vd_edge.cell();
    const vd_t::cell_type* right_cell = vd_edge.twin()->cell();
    Point start = VoronoiUtils::p(vd_edge.vertex0());
    Point end = VoronoiUtils::p(vd_edge.vertex1());
    
    bool point_left = left_cell->contains_point();
    bool point_right = right_cell->contains_point();
    if ((!point_left && !point_right)
        || vd_edge.is_secondary() // source vert is directly connected to source segment
    )
    {
        return std::vector<Point>({ start, end });
    }
    else if (point_left == !point_right)
    {
        const vd_t::cell_type* point_cell = left_cell;
        const vd_t::cell_type* segment_cell = right_cell;
        if (!point_left)
        {
            std::swap(point_cell, segment_cell);
        }
        Point p = VoronoiUtils::getSourcePoint(*point_cell, points, segments);
        const Segment& s = VoronoiUtils::getSourceSegment(*segment_cell, points, segments);
        return VoronoiUtils::discretizeParabola(p, s, start, end, discretization_step_size, transitioning_angle);
    }
    else
    {
        Point left_point = VoronoiUtils::getSourcePoint(*left_cell, points, segments);
        Point right_point = VoronoiUtils::getSourcePoint(*right_cell, points, segments);
        coord_t d = vSize(right_point - left_point);
        Point middle = (left_point + right_point) / 2;
        Point x_axis_dir = turn90CCW(right_point - left_point);
        coord_t x_axis_length = vSize(x_axis_dir);

        const auto projected_x = [x_axis_dir, x_axis_length, middle](Point from)
        {
            Point vec = from - middle;
            coord_t x = dot(vec, x_axis_dir) / x_axis_length;
            return x;
        };
        coord_t start_x = projected_x(start);
        coord_t end_x = projected_x(end);

        float bound = 0.5 / tan((M_PI - transitioning_angle) * 0.5);
        coord_t marking_start_x = - d * bound;
        coord_t marking_end_x = d * bound;
        Point marking_start = middle + x_axis_dir * marking_start_x / x_axis_length;
        Point marking_end = middle + x_axis_dir * marking_end_x / x_axis_length;
        coord_t dir = 1;
        if (start_x > end_x)
        {
            dir = -1;
            std::swap(marking_start, marking_end);
            std::swap(marking_start_x, marking_end_x);
        }

        Point a = start;
        Point b = end;
        std::vector<Point> ret;
        ret.emplace_back(a);
        
        bool add_marking_start = marking_start_x * dir > start_x * dir;
        bool add_marking_end = marking_end_x * dir > start_x * dir;
        
        Point ab = b - a;
        coord_t ab_size = vSize(ab);
        coord_t step_count = (ab_size + discretization_step_size / 2) / discretization_step_size;
        if (step_count % 2 == 1)
        {
            step_count++; // enforce a discretization point being added in the middle
        }
        for (coord_t step = 1; step < step_count; step++)
        {
            Point here = a + ab * step / step_count;
            coord_t x_here = projected_x(here);
            if (add_marking_start && marking_start_x * dir < x_here * dir)
            {
                ret.emplace_back(marking_start);
                add_marking_start = false;
            }
            if (add_marking_end && marking_end_x * dir < x_here * dir)
            {
                ret.emplace_back(marking_end);
                add_marking_end = false;
            }
            ret.emplace_back(here);
        }
        if (add_marking_end && marking_end_x * dir < end_x * dir)
        {
            ret.emplace_back(marking_end);
        }
        ret.emplace_back(b);
        return ret;
    }
}


bool SkeletalTrapezoidation::computePointCellRange(vd_t::cell_type& cell, Point& start_source_point, Point& end_source_point, vd_t::edge_type*& starting_vd_edge, vd_t::edge_type*& ending_vd_edge, const std::vector<Point>& points, const std::vector<Segment>& segments)
{
    if (cell.incident_edge()->is_infinite())
    {
        return false;
    }
    // check if any point of the cell is inside or outside polygon
    // copy whole cell into graph or not at all
    
    const Point source_point = VoronoiUtils::getSourcePoint(cell, points, segments);
    const PolygonsPointIndex source_point_index = VoronoiUtils::getSourcePointIndex(cell, points, segments);
    Point some_point = VoronoiUtils::p(cell.incident_edge()->vertex0());
    if (some_point == source_point)
    {
        some_point = VoronoiUtils::p(cell.incident_edge()->vertex1());
    }
    if (!LinearAlg2D::isInsideCorner(source_point_index.prev().p(), source_point_index.p(), source_point_index.next().p(), some_point))
    { // cell is outside of polygon
        return false; // don't copy any part of this cell
    }
    bool first = true;
    for (vd_t::edge_type* vd_edge = cell.incident_edge(); vd_edge != cell.incident_edge() || first; vd_edge = vd_edge->next())
    {
        assert(vd_edge->is_finite());
        Point p1 = VoronoiUtils::p(vd_edge->vertex1());
        if (p1 == source_point)
        {
            start_source_point = source_point;
            end_source_point = source_point;
            starting_vd_edge = vd_edge->next();
            ending_vd_edge = vd_edge;
        }
        else
            assert((VoronoiUtils::p(vd_edge->vertex0()) == source_point || !vd_edge->is_secondary()) && "point cells must end in the point! They cannot cross the point with an edge, because collinear edges are not allowed in the input.");
        first = false;
    }
    assert(starting_vd_edge && ending_vd_edge);
    assert(starting_vd_edge != ending_vd_edge);
    return true;
}
void SkeletalTrapezoidation::computeSegmentCellRange(vd_t::cell_type& cell, Point& start_source_point, Point& end_source_point, vd_t::edge_type*& starting_vd_edge, vd_t::edge_type*& ending_vd_edge, const std::vector<Point>& points, const std::vector<Segment>& segments)
{
    const Segment& source_segment = VoronoiUtils::getSourceSegment(cell, points, segments);
    Point from = source_segment.from();
    Point to = source_segment.to();

    // find starting edge
    // find end edge
    bool first = true;
    bool seen_possible_start = false;
    bool after_start = false;
    bool ending_edge_is_set_before_start = false;
    for (vd_t::edge_type* edge = cell.incident_edge(); edge != cell.incident_edge() || first; edge = edge->next())
    {
        if (edge->is_infinite())
        {
            first = false;
            continue;
        }
        bool check_secondary_edge = true;
        Point v0 = VoronoiUtils::p(edge->vertex0());
        Point v1 = VoronoiUtils::p(edge->vertex1());
        assert(!(v0 == to && v1 == from));
        if (v0 == to
            && !after_start // use the last edge which starts in source_segment.to
        )
        {
            starting_vd_edge = edge;
            seen_possible_start = true;
            check_secondary_edge = false;
        }
        else if (seen_possible_start)
            after_start = true;
        if (v1 == from
            && (!ending_vd_edge || ending_edge_is_set_before_start)
        )
        {
            ending_edge_is_set_before_start = !after_start;
            ending_vd_edge = edge;
            check_secondary_edge = false;
        }
        if (false && check_secondary_edge && edge->is_secondary()
               != LinearAlg2D::pointLiesOnTheRightOfLine(v1, from, to)
            &&    LinearAlg2D::pointLiesOnTheRightOfLine(v0, from, to)
            && (LinearAlg2D::getDist2FromLineSegment(v0, from, v1) <= 5 // TODO: magic value
                || LinearAlg2D::getDist2FromLineSegment(v0, to, v1) <= 5)
        )
        { // edge crosses source segment
            // TODO: handle the case where two consecutive line segments are collinear!
            // that's the only case where a voronoi segment doesn't end in a polygon vertex, but goes though it
            if (LinearAlg2D::pointLiesOnTheRightOfLine(VoronoiUtils::p(edge->vertex1()), source_segment.from(), source_segment.to()))
            {
                ending_vd_edge = edge;
            }
            else
            {
                starting_vd_edge = edge;
            }
            first = false;
            continue;
        }
        first = false;
    }
    
    assert(starting_vd_edge && ending_vd_edge);
    assert(starting_vd_edge != ending_vd_edge);
    
    start_source_point = source_segment.to();
    end_source_point = source_segment.from();
}

void SkeletalTrapezoidation::initialize_graph()
{
    std::vector<Point> points; // remains empty

    std::vector<Segment> segments;
    for (size_t poly_idx = 0; poly_idx < polys.size(); poly_idx++)
    {
        ConstPolygonRef poly = polys[poly_idx];
        for (size_t point_idx = 0; point_idx < poly.size(); point_idx++)
        {
            segments.emplace_back(&polys, poly_idx, point_idx);
        }
    }

    vd_t vd;
    construct_voronoi(points.begin(), points.end(),
                                        segments.begin(), segments.end(),
                                        &vd);

#ifdef DEBUG
    VoronoiUtils::debugOutput("output/vd.svg", vd, points, segments);
#endif

    for (const vd_t::edge_type& edge : vd.edges())
    {
        assert(edge.vertex0() == edge.twin()->vertex1());
        assert(edge.vertex1() == edge.twin()->vertex0());
        assert(edge.vertex1() == edge.next()->vertex0());
        assert(edge.vertex0() == edge.prev()->vertex1());
    }

    
    for (vd_t::cell_type cell : vd.cells())
    {
        if (!cell.incident_edge())
        { // there is no spoon
            continue;
        }
        Point start_source_point;
        Point end_source_point;
        vd_t::edge_type* starting_vd_edge = nullptr;
        vd_t::edge_type* ending_vd_edge = nullptr;
        // compute and store result in above variables
        
        if (cell.contains_point())
        {
            bool keep_going = computePointCellRange(cell, start_source_point, end_source_point, starting_vd_edge, ending_vd_edge, points, segments);
            if (!keep_going)
            {
                continue;
            }
        }
        else
        {
            computeSegmentCellRange(cell, start_source_point, end_source_point, starting_vd_edge, ending_vd_edge, points, segments);
        }
        
        if (!starting_vd_edge || !ending_vd_edge)
        {
            assert(false && "each cell should start / end in a polygon vertex");
            continue;
        }
        
        // copy start to end edge to graph
        
        edge_t* prev_edge = nullptr;
        transfer_edge(start_source_point, VoronoiUtils::p(starting_vd_edge->vertex1()), *starting_vd_edge, prev_edge, start_source_point, end_source_point, points, segments);
        node_t* starting_node = vd_node_to_he_node[starting_vd_edge->vertex0()];
        starting_node->data.distance_to_boundary = 0;
        // starting_edge->prev = nullptr;
//         starting_edge->from->data.distance_to_boundary = 0; // TODO

        make_rib(prev_edge, start_source_point, end_source_point);
        for (vd_t::edge_type* vd_edge = starting_vd_edge->next(); vd_edge != ending_vd_edge; vd_edge = vd_edge->next())
        {
            assert(vd_edge->is_finite());
            Point v1 = VoronoiUtils::p(vd_edge->vertex0());
            Point v2 = VoronoiUtils::p(vd_edge->vertex1());
            transfer_edge(v1, v2, *vd_edge, prev_edge, start_source_point, end_source_point, points, segments);

            make_rib(prev_edge, start_source_point, end_source_point);
        }

        transfer_edge(VoronoiUtils::p(ending_vd_edge->vertex0()), end_source_point, *ending_vd_edge, prev_edge, start_source_point, end_source_point, points, segments);
        // ending_edge->next = nullptr;
        prev_edge->to->data.distance_to_boundary = 0;

        debugCheckGraphConsistency(true);
    }

    
    debugCheckGraphCompleteness();
    debugCheckGraphConsistency();
    debugCheckGraphExistance();

    separatePointyQuadEndNodes();

    debugCheckGraphCompleteness();
    debugCheckGraphConsistency();
    debugCheckGraphExistance();

    fixNodeDuplication();

    debugCheckGraphCompleteness();
    debugCheckGraphConsistency();
    debugCheckGraphExistance();
    debugCheckEndpointUniqueness();
    
    collapseSmallEdges();


    debugCheckGraphCompleteness();
    debugCheckGraphConsistency();
    debugCheckGraphExistance();

    { // set [some_edge] the the first possible edge
        // that way we can iterate over all reachable edges from node.some_edge without needing to iterate backward
        for (edge_t& edge : graph.edges)
        {
            if (!edge.prev)
            {
                edge.from->some_edge = &edge;
            }
        }
    }


    debugCheckGraphCompleteness();
    debugCheckGraphConsistency();
    debugCheckGraphStructure();
    debugCheckGraphReachability();
    debugCheckGraphExistance();
    
#ifdef DEBUG
    {
        AABB aabb(polys);
        SVG svg("output/vq.svg", aabb);
        debugOutput(svg, false, false);
    }
    {
        AABB aabb(polys);
        SVG svg("output/radial_dists.svg", aabb);
        debugOutput(svg, false, true);
    }
    {
        AABB aabb(polys);
        SVG svg("output/bead_counts.svg", aabb);
        debugOutput(svg, false, false, true);
    }
    {
        AABB aabb(polys);
        SVG svg("output/locations.svg", aabb);
        debugOutput(svg, false, false, false, true);
    }
    debugCheckGraphCompleteness();
    debugCheckGraphConsistency();

    if (generate_MAT_STL)
    {
        STLwriter stl("output/mat.stl");
        debugOutputSTL(stl);
    }
#endif

    vd_edge_to_he_edge.clear();
    vd_node_to_he_node.clear();
}

void SkeletalTrapezoidation::separatePointyQuadEndNodes()
{
    std::unordered_set<node_t*> visited_nodes;
    for (edge_t& edge : graph.edges)
    {
        if (edge.prev) continue;
        edge_t* quad_start = &edge;
        if (visited_nodes.find(quad_start->from) == visited_nodes.end())
        {
            visited_nodes.emplace(quad_start->from);
        }
        else
        { // needs to be duplicated
            graph.nodes.emplace_back(*quad_start->from);
            node_t* new_node = &graph.nodes.back();
            new_node->some_edge = quad_start;
            quad_start->from = new_node;
            quad_start->twin->to = new_node;
        }
    }

    debugCheckEndpointUniqueness();
}

void SkeletalTrapezoidation::collapseSmallEdges(coord_t snap_dist)
{
    std::unordered_map<edge_t*, std::list<edge_t>::iterator> edge_locator;
    std::unordered_map<node_t*, std::list<node_t>::iterator> node_locator;
    for (auto edge_it = graph.edges.begin(); edge_it != graph.edges.end(); ++edge_it)
        edge_locator.emplace(&*edge_it, edge_it);
    for (auto node_it = graph.nodes.begin(); node_it != graph.nodes.end(); ++node_it )
        node_locator.emplace(&*node_it, node_it );
    auto safelyRemoveEdge = [this, &edge_locator](edge_t* to_be_removed, std::list<edge_t>::iterator& current_edge_it, bool& edge_it_is_updated)
        {
            if (current_edge_it != graph.edges.end()
                && to_be_removed == &*current_edge_it)
            {
                current_edge_it = graph.edges.erase(current_edge_it);
                edge_it_is_updated = true;
            }
            else
            {
                graph.edges.erase(edge_locator[to_be_removed]);
            }
        };

    auto should_collapse = [snap_dist](node_t* a, node_t* b) { return shorterThen(a->p - b->p, snap_dist); };
        
    for (auto edge_it = graph.edges.begin(); edge_it != graph.edges.end();)
    {
        if (edge_it->prev)
        {
            edge_it++;
            continue;
        }
        edge_t* quad_start = &*edge_it;
        edge_t* quad_end = quad_start; while (quad_end->next) quad_end = quad_end->next;
        edge_t* quad_mid = (quad_start->next == quad_end)? nullptr : quad_start->next;

        bool edge_it_is_updated = false;
        bool quad_mid_is_removed = false;
        if (quad_mid && should_collapse(quad_mid->from, quad_mid->to))
        {
            assert(quad_mid->twin);
            int count = 0;
            for (edge_t* edge_from_3 = quad_end; edge_from_3 && edge_from_3 != quad_mid->twin; edge_from_3 = edge_from_3->twin->next)
            {
                edge_from_3->from = quad_mid->from;
                edge_from_3->twin->to = quad_mid->from;
                if (count > 50)
                {
                    std::cerr << edge_from_3->from->p << " - " << edge_from_3->to->p << '\n';
                }
                assert(++count < 100);
                if (count > 1000) break;
            }

            // o-o > collapse top
            // | |
            // | |
            // | |
            // o o
            if (quad_mid->from->some_edge == quad_mid)
            {
                if (quad_mid->twin->next)
                {
                    quad_mid->from->some_edge = quad_mid->twin->next;
                }
                else
                {
                    quad_mid->from->some_edge = quad_mid->prev->twin;
                }
            }
//             if (quad_mid->twin->from->some_edge == quad_mid->twin)
//             {
//                 quad_mid->twin->from->some_edge = quad_mid->next;
//             }
            graph.nodes.erase(node_locator[quad_mid->to]);

            quad_mid->prev->next = quad_mid->next;
            quad_mid->next->prev = quad_mid->prev;
            quad_mid->twin->next->prev = quad_mid->twin->prev;
            quad_mid->twin->prev->next = quad_mid->twin->next;

            safelyRemoveEdge(quad_mid->twin, edge_it, edge_it_is_updated);
            safelyRemoveEdge(quad_mid, edge_it, edge_it_is_updated);
            quad_mid_is_removed = true;
        }

        //  o-o
        //  | | > collapse sides
        //  o o
        if ( should_collapse(quad_start->from, quad_end->to)
            && should_collapse(quad_start->to, quad_end->from))
        { // collapse start and end edges and remove whole cell
            assert(!quad_mid || quad_mid_is_removed);

            quad_start->twin->to = quad_end->to;
            quad_end->to->some_edge = quad_end->twin;
            if (quad_end->from->some_edge == quad_end)
            {
                if (quad_end->twin->next)
                {
                    quad_end->from->some_edge = quad_end->twin->next;
                }
                else
                {
                    quad_end->from->some_edge = quad_end->prev->twin;
                }
            }
            graph.nodes.erase(node_locator[quad_start->from]);

            quad_start->twin->twin = quad_end->twin;
            quad_end->twin->twin = quad_start->twin;
            safelyRemoveEdge(quad_start, edge_it, edge_it_is_updated);
            safelyRemoveEdge(quad_end, edge_it, edge_it_is_updated);
        }
        // if only one side had collapsable length then the cell on the other side of that edge has to collapse
        // if we would collapse that one edge then that would change the quad_start and/or quad_end of neighboring cells
        // this is to do with the constraint that !prev == !twin.next

        if (!edge_it_is_updated)
        {
            edge_it++;
        }
    }
}

void SkeletalTrapezoidation::fixNodeDuplication()
{ // fix duplicate verts
    for (auto node_it = graph.nodes.begin(); node_it != graph.nodes.end();)
    {
        node_t* replacing_node = nullptr;
        for (edge_t* outgoing = node_it->some_edge; outgoing != node_it->some_edge; outgoing = outgoing->twin->next)
        {
            assert(outgoing);
            if (outgoing->from != &*node_it)
            {
                replacing_node = outgoing->from;
            }
            if (outgoing->twin->to != &*node_it)
            {
                replacing_node = outgoing->twin->to;
            }
        }
        if (replacing_node)
        {
            for (edge_t* outgoing = node_it->some_edge; outgoing != node_it->some_edge; outgoing = outgoing->twin->next)
            {
                outgoing->twin->to = replacing_node;
                outgoing->from = replacing_node;
            }
            node_it = graph.nodes.erase(node_it);
        }
        else
        {
            ++node_it;
        }
    }
}

//
// ^^^^^^^^^^^^^^^^^^^^^
//    INITIALIZATION
// =====================
//
// =====================
//       HELPERS
// vvvvvvvvvvvvvvvvvvvvv
//


void SkeletalTrapezoidation::debugCheckGraphCompleteness()
{
#ifdef DEBUG
    for (const node_t& node : graph.nodes)
    {
        if (!node.some_edge)
        {
            assert(false);
        }
    }
    for (const edge_t& edge : graph.edges)
    {
        if (!edge.twin || !edge.from || !edge.to)
        {
            assert(false);
        }
        assert((edge.next == nullptr) == (edge.twin->prev == nullptr));
        assert((edge.prev == nullptr) == (edge.twin->next == nullptr));
        assert(edge.next || edge.to->data.distance_to_boundary == 0);
        assert(edge.prev || edge.from->data.distance_to_boundary == 0);
    }
#endif
}

void SkeletalTrapezoidation::debugCheckEndpointUniqueness()
{
#ifdef DEBUG
    for (edge_t& edge : graph.edges)
    {
        if (edge.prev) continue;
        for (edge_t& e2 : graph.edges)
        {
            assert(e2.from != edge.from || edge == e2);
        }
    }
#endif
}

void SkeletalTrapezoidation::debugCheckGraphExistance()
{
#ifdef DEBUG
    std::unordered_set<edge_t*> all_edges;
    std::unordered_set<node_t*> all_nodes;
    for (node_t& node : graph.nodes)
    {
        all_nodes.emplace(&node);
    }
    for (edge_t& edge : graph.edges)
    {
        all_edges.emplace(&edge);
    }
    
    auto edge_exists = [&all_edges](edge_t* edge)
        {
            assert(edge == nullptr ||
                all_edges.find(edge) != all_edges.end());
        };
    auto node_exists = [&all_nodes](node_t* node)
        {
            assert(node == nullptr ||
                all_nodes.find(node) != all_nodes.end());
        };
    for (node_t& node : graph.nodes)
    {
        edge_exists(node.some_edge);
    }
    for (edge_t& edge : graph.edges)
    {
        edge_t* edge_ = &edge;
        edge_exists(edge.prev);
        edge_exists(edge.next);
        edge_exists(edge.twin);
        node_exists(edge.from);
        node_exists(edge.to);
    }
#endif
}

void SkeletalTrapezoidation::debugCheckGraphStructure()
{
#ifdef DEBUG
    for (edge_t& edge : graph.edges)
    {
        node_t* node = edge.from;
        bool first = true;
        coord_t count = 0;
        bool seen_edge = false;
        for (const edge_t* outgoing = node->some_edge; outgoing && (first || outgoing != node->some_edge); outgoing = outgoing->twin->next)
        {
            assert(++count < 100);
            if (outgoing == &edge) seen_edge = true;
            first = false;
            if (!outgoing->twin) break;
        }
        assert(seen_edge);
    }
#endif
}

void SkeletalTrapezoidation::debugCheckGraphReachability()
{
#ifdef DEBUG
    std::unordered_set<node_t*> reachable_nodes;
    for (edge_t& edge : graph.edges)
    {
        reachable_nodes.emplace(edge.from);
        reachable_nodes.emplace(edge.to);
    }
    for (node_t& node : graph.nodes)
    {
        assert(reachable_nodes.find(&node) != reachable_nodes.end());
    }

    std::unordered_set<edge_t*> reachable_edges;
    for (node_t& node : graph.nodes)
    {
        bool first = true;
        for (edge_t* outgoing = node.some_edge; outgoing && (first || outgoing != node.some_edge); outgoing = outgoing->twin->next)
        {
            reachable_edges.emplace(outgoing);
            first = false;
            if (!outgoing->twin) break;
        }
    }
    for (edge_t& edge : graph.edges)
    {
        edge_t* edge_ = &edge;
        if (reachable_edges.find(&edge) == reachable_edges.end())
        {
            std::cerr << "Cannot find " << edge.from->p << " - " << edge.to->p << " among edges around former!\n";
            bool seen = false;
            bool first = true;
            for (edge_t* outgoing = edge.from->some_edge; outgoing && (first || outgoing != edge.from->some_edge); outgoing = outgoing->twin->next)
            {
                std::cerr << outgoing->to->p << "\n";
                if (outgoing == edge_)
                {
                    seen = true;
                }
                first = false;
                if (!outgoing->twin) break;
            }
            assert(seen);
            assert(std::find(graph.nodes.begin(), graph.nodes.end(), *edge.from) != graph.nodes.end());
        }
        assert(reachable_edges.find(&edge) != reachable_edges.end());
    }
#endif
}

void SkeletalTrapezoidation::debugCheckGraphConsistency(bool ignore_duplication)
{
#ifdef DEBUG
    auto vert_assert = [ignore_duplication](const node_t* first, const node_t* second)
    {
        if (first != second)
        {
            if (first->p == second->p)
            {
                RUN_ONCE(logWarning("Unneccesary duplicatation of SkeletalTrapezoidation nodes!\n"));
                assert(ignore_duplication);
            }
            else
            {
                assert(false && "connected edges don't refer to the same node!");
            }
        }
    };
    
    for (const edge_t& edge : graph.edges)
    {
        const edge_t* edge_p = &edge;
        if (edge_p->twin)
        {
            if (!edge_p->to)
            {
                assert(!edge_p->twin->from);
                vert_assert(edge_p->twin->from, edge_p->to);
            }
            if (!edge_p->from)
            {
                assert(!edge_p->twin->to);
                vert_assert(edge_p->twin->to, edge_p->from);
            }
            assert((edge_p->from == nullptr) == (edge_p->twin->to == nullptr));
            assert((edge_p->to == nullptr) == (edge_p->twin->from == nullptr));
            assert(edge_p->twin->twin == &edge);
            assert(edge_p->twin != edge_p);
        }
        if (edge_p->next)
        {
            vert_assert(edge_p->next->from, edge_p->to);
        }
        if (edge_p->prev)
        {
            vert_assert(edge_p->prev->to, edge_p->from);
        }
    }
    for (const node_t& node : graph.nodes)
    {
        if (node.some_edge)
        {
            const node_t* node_ = &node;
            vert_assert(node.some_edge->from, &node);
            if (node.some_edge->twin)
            {
                for (const edge_t* outgoing = node.some_edge->twin->next; outgoing && outgoing->twin && outgoing != node.some_edge; outgoing = outgoing->twin->next)
                {
                    vert_assert(outgoing->from, &node);
                }
            }
        }
    }

    for (const edge_t& edge : graph.edges)
    {
        int i = 0;
        for (const edge_t* e = &edge; e; e = e->next) assert(++i < 10);
        for (const edge_t* e = &edge; e; e = e->prev) assert(++i < 10);
    }
#endif // DEBUG
}


SVG::ColorObject SkeletalTrapezoidation::getColor(edge_t& edge)
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

void SkeletalTrapezoidation::debugOutput(SVG& svg, bool draw_arrows, bool draw_dists, bool draw_bead_counts, bool draw_locations)
{
    svg.nextLayer();
    svg.writeAreas(polys, SVG::Color::NONE, SVG::Color::BLACK, 3);
    svg.nextLayer();
    for (edge_t& edge : graph.edges)
    {
        Point a = edge.from->p;
        Point b = edge.to->p;
        SVG::ColorObject clr = getColor(edge);
        float stroke_width = 1;
        if (edge.data.markingIsSet() && edge.data.isMarked())
        {
            continue;
        }
        if (draw_arrows)
        {
            constexpr coord_t spacing_dist = 10;
            svg.writeArrow(a + normal(b - a, spacing_dist * 3), b - normal(b - a, spacing_dist), clr, stroke_width, 10, spacing_dist);
        }
        else
        {
            if (edge.to->p < edge.from->p)
            {
                svg.writeLine(a, b, clr, stroke_width);
            }
        }
    }
    svg.nextLayer();
    for (edge_t& edge : graph.edges)
    {
        Point a = edge.from->p;
        Point b = edge.to->p;
        SVG::ColorObject clr = SVG::Color::BLUE;
        float stroke_width = 2;
        if ( ! edge.data.markingIsSet() || ! edge.data.isMarked())
        {
            continue;
        }
        if (draw_arrows)
        {
            constexpr coord_t spacing_dist = 10;
            svg.writeArrow(a + normal(b - a, spacing_dist * 3), b - normal(b - a, spacing_dist), clr, stroke_width, 10, spacing_dist);
        }
        else
        {
            if (edge.to->p < edge.from->p)
            {
                svg.writeLine(a, b, clr, stroke_width);
            }
        }
    }
    svg.nextLayer();
    for (node_t& node : graph.nodes)
    {
        if (draw_arrows)
        {
            svg.writePoint(node.p);
        }
        if (draw_dists && node.data.distance_to_boundary > 0)
        {
            svg.writeText(node.p, std::to_string(node.data.distance_to_boundary));
        }
        if (draw_bead_counts && node.data.bead_count >= 0)
        {
            if (node.data.transition_ratio != 0)
            {
                std::ostringstream ss;
                ss.precision(2);
                ss << (node.data.transition_ratio + node.data.bead_count);
                svg.writeText(node.p, ss.str());
            }
            else
            {
                svg.writeText(node.p, std::to_string(node.data.bead_count));
            }
        }
        if (draw_locations)
        {
            svg.writePoint(node.p, false, 1);
            std::ostringstream ss;
            ss << node.p.X << "," << node.p.Y;
            svg.writeText(node.p, ss.str(), SVG::Color::BLACK, 4);
        }
    }
    svg.nextLayer();
}

void SkeletalTrapezoidation::debugOutputSTL(STLwriter& stl, bool use_bead_count)
{
    auto toPoint3 = [use_bead_count](node_t* node)
        {
            coord_t h = use_bead_count? std::max(float(0), node->data.bead_count + node->data.transition_ratio) * 200 : node->data.distance_to_boundary;
            return Point3(node->p.X, node->p.Y, h);
        };
    for (edge_t& edge : graph.edges)
    {
        if (edge.prev)
        {
            continue;
        }
        
        assert(edge.next);
        if (edge.next->next)
        {
            stl.writeQuad(toPoint3(edge.next->next->to), toPoint3(edge.next->to), toPoint3(edge.from), toPoint3(edge.to));
//             stl.writeTriangle(toPoint3(edge.from), toPoint3(edge.next->to), toPoint3(edge.next->next->to));
//             stl.writeTriangle(toPoint3(edge.from), toPoint3(edge.to), toPoint3(edge.next->to));
        }
        else
        {
            stl.writeTriangle(toPoint3(edge.from), toPoint3(edge.to), toPoint3(edge.next->to));
        }
    }
}

} // namespace arachne
