//Copyright (c) 2019 Ultimaker B.V.


#ifndef UTILS_HALF_EDGE_GRAPH_H
#define UTILS_HALF_EDGE_GRAPH_H


#include <list>
#include <cassert>



#include "HalfEdge.h"
#include "HalfEdgeNode.h"
#include "SVG.h"

namespace arachne
{

template<
typename edge_t,
typename node_t,
class = typename std::enable_if<std::is_base_of<HalfEdge<edge_t, node_t>, edge_t>::value>,
class = typename std::enable_if<std::is_base_of<HalfEdgeNode<edge_t, node_t>, node_t>::value>
>
class HalfEdgeGraph
{
public:
    std::list<edge_t> edges;
    std::list<node_t> nodes;

    void debugOutput(std::string filename);
    void debugOutput(SVG& svg);

    bool bedugCheckDataCompleteness() const;
};



template<typename edge_t, typename node_t, typename edge_t_enabled, typename node_t_enabled>
void HalfEdgeGraph<edge_t, node_t, edge_t_enabled, node_t_enabled>::debugOutput(std::string filename)
{
    AABB aabb;
    for (node_t& node : nodes)
    {
        aabb.include(node.p);
    }
    SVG svg(filename.c_str(), aabb);
    debugOutput(svg);
}

template<typename edge_t, typename node_t, typename edge_t_enabled, typename node_t_enabled>
void HalfEdgeGraph<edge_t, node_t, edge_t_enabled, node_t_enabled>::debugOutput(SVG& svg)
{
//     for (node_t& node : nodes)
//     {
//         svg.writePoint(node.p);
//     }
    for (edge_t& edge : edges)
    {
        svg.writeLine(edge.from->p, edge.to->p, SVG::Color::RED);
    }
}

template<typename edge_t, typename node_t, typename edge_t_enabled, typename node_t_enabled>
bool HalfEdgeGraph<edge_t, node_t, edge_t_enabled, node_t_enabled>::bedugCheckDataCompleteness() const
{
    size_t problems = 0;
    for (const node_t& node : nodes)
    {
        if (!node.some_edge)
        {
            problems++;
            assert(false);
        }
    }
    for (const edge_t& edge : edges)
    {
        if (!edge.twin || !edge.next || !edge.prev || !edge.from || !edge.to)
        {
            problems++;
            assert(false);
        }
    }
    
    assert(problems == 0);
    return problems == 0;
}



} // namespace arachne
#endif // UTILS_HALF_EDGE_GRAPH_H
