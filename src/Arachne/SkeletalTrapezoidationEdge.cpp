//Copyright (c) 2020 Ultimaker B.V.

#include "SkeletalTrapezoidationEdge.h"
#include "SkeletalTrapezoidationNode.h"

namespace arachne
{

using edge_t = SkeletalTrapezoidationEdge;
using node_t = SkeletalTrapezoidationNode;

    
bool SkeletalTrapezoidationEdge::canGoUp(bool strict) const
{
    if (to->distance_to_boundary > from->distance_to_boundary)
    {
        return true;
    }
    if (to->distance_to_boundary < from->distance_to_boundary
        || strict
    )
    {
        return false;
    }
    // edge is between equidistqant verts; recurse!
    for (const edge_t* outgoing = next; outgoing != twin; outgoing = outgoing->twin->next)
    {
        if (outgoing->canGoUp())
        {
            return true;
        }
        assert(outgoing->twin); if (!outgoing->twin) return false;
        assert(outgoing->twin->next); if (!outgoing->twin->next) return true; // This point is on the boundary?! Should never occur
    }
    return false;
}

std::optional<coord_t> SkeletalTrapezoidationEdge::distToGoUp() const
{
    if (to->distance_to_boundary > from->distance_to_boundary)
    {
        return 0;
    }
    if (to->distance_to_boundary < from->distance_to_boundary)
    {
        return std::optional<coord_t>();
    }
    // edge is between equidistqant verts; recurse!
    std::optional<coord_t> ret;
    for (edge_t* outgoing = next; outgoing != twin; outgoing = outgoing->twin->next)
    {
        std::optional<coord_t> dist_to_up = outgoing->distToGoUp();
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
        ret =  *ret + vSize(to->p - from->p);
    }
    return ret;
}


bool SkeletalTrapezoidationEdge::isUpward() const
{
    if (to->distance_to_boundary > from->distance_to_boundary)
    {
        return true;
    }
    if (to->distance_to_boundary < from->distance_to_boundary)
    {
        return false;
    }
    // equidistant edge case:
    std::optional<coord_t> forward_up_dist = distToGoUp();
    std::optional<coord_t> backward_up_dist = twin->distToGoUp();
    if (forward_up_dist && backward_up_dist)
    {
        return forward_up_dist < backward_up_dist;
    }
    if (forward_up_dist) return true;
    if (backward_up_dist) return false;
    return to->p < from->p; // arbitrary ordering, which returns the opposite for the twin edge
}

} // namespace arachne
