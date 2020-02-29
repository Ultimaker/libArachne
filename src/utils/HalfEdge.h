//Copyright (c) 2019 Ultimaker B.V.


#ifndef UTILS_HALF_EDGE_H
#define UTILS_HALF_EDGE_H

#include <forward_list>

namespace arachne
{

template<
typename edge_t,
typename node_t
>
class HalfEdgeNode;

template<
typename edge_t,
typename node_t
>
class HalfEdge
{
public:
    edge_t* twin = nullptr;
    edge_t* next = nullptr;
    edge_t* prev = nullptr;
    node_t* from = nullptr;
    node_t* to = nullptr;
    HalfEdge()
    {}
    bool operator==(const HalfEdge& other)
    {
        return this == &other;
    }
};




} // namespace arachne
#endif // UTILS_HALF_EDGE_H
