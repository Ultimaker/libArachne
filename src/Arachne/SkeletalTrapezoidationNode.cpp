//Copyright (c) 2020 Ultimaker B.V.

#include "SkeletalTrapezoidationNode.h"
#include "SkeletalTrapezoidationEdge.h"

namespace arachne
{

using edge_t = SkeletalTrapezoidationEdge;
using node_t = SkeletalTrapezoidationNode;

bool SkeletalTrapezoidationNode::isMarked() const
{
    bool first = true;
    for (edge_t* edge = some_edge; first || edge != some_edge; edge = edge->twin->next)
    {
        if (edge->isMarked())
        {
            return true;
        }
        first = false;
        assert(edge->twin); if (!edge->twin) return false;
    }
    return false;
}


bool SkeletalTrapezoidationNode::isLocalMaximum(bool strict) const
{
    if (distance_to_boundary == 0)
    {
        return false;
    }
    bool first = true;
    for (edge_t* edge = some_edge; first || edge != some_edge; edge = edge->twin->next)
    {
        if (edge->canGoUp(strict))
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

} // namespace arachne
