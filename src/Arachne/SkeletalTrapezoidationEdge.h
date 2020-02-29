//Copyright (c) 2019 Ultimaker B.V.

#ifndef ARACHNE_SKELETAL_TRAPEZOIDATION_EDGE_H
#define ARACHNE_SKELETAL_TRAPEZOIDATION_EDGE_H

#include "utils/HalfEdge.h"
#include "utils/HalfEdgeNode.h"

namespace arachne
{

class SkeletalTrapezoidationJoint;

class SkeletalTrapezoidationEdge : public HalfEdge<SkeletalTrapezoidationEdge, SkeletalTrapezoidationJoint>
{
public:
    enum Type : int_least16_t
    {
        NORMAL = 0,        // from voronoi diagram
        EXTRA_VD = 1,      // introduced to voronoi diagram in order to make the skeletal trapezoidation
        TRANSITION_END = 2 // introduced to voronoi diagram in order to support transitions
    };
    Type type;

    SkeletalTrapezoidationEdge()
    : SkeletalTrapezoidationEdge(NORMAL)
    {}
    SkeletalTrapezoidationEdge(Type type)
    : HalfEdge<SkeletalTrapezoidationEdge, SkeletalTrapezoidationJoint>()
    , type(type)
    , is_marked(-1)
    {}

    bool isMarked() const
    {
        assert(is_marked != -1);
        return is_marked;
    }
    void setMarked(bool b)
    {
        is_marked = b;
    }
    bool markingIsSet() const
    {
        return is_marked >= 0;
    }
private:
    int_least8_t is_marked; //! whether the edge is significant; whether the source segments have a sharp angle; -1 is unknown
};




} // namespace arachne
#endif // ARACHNE_SKELETAL_TRAPEZOIDATION_EDGE_H
