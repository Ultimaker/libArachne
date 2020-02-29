//Copyright (c) 2019 Ultimaker B.V.


#ifndef UTILS_HALF_EDGE_NODE_H
#define UTILS_HALF_EDGE_NODE_H

#include <list>

#include "IntPoint.h"

namespace arachne
{

template<
typename edge_t,
typename node_t
>
class HalfEdge;

template<
typename edge_t,
typename node_t
>
class HalfEdgeNode
{
public:
    Point p;
    edge_t* some_edge = nullptr;
    HalfEdgeNode(Point p)
    : p(p)
    {}
    bool operator==(const HalfEdgeNode& other)
    {
        return this == &other;
    }
};




} // namespace arachne
#endif // UTILS_HALF_EDGE_NODE_H
