//Copyright (c) 2019 Ultimaker B.V.

#ifndef ARACHNE_SKELETAL_TRAPEZOIDATION_EDGE_H
#define ARACHNE_SKELETAL_TRAPEZOIDATION_EDGE_H

#include <cassert>

#include "utils/HalfEdge.h"
#include "utils/HalfEdgeNode.h"

namespace arachne
{

class SkeletalTrapezoidationNode;

class SkeletalTrapezoidationEdge : public HalfEdge<SkeletalTrapezoidationEdge, SkeletalTrapezoidationNode>
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
    : HalfEdge<SkeletalTrapezoidationEdge, SkeletalTrapezoidationNode>()
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
    
    /*!
     * Check (recursively) whether there is any upward edge from the distance_to_boundary of the from of the \param edge
     * 
     * \param strict Whether equidistant edges can count as a local maximum
     */
    bool canGoUp(bool strict = false) const;

    /*!
     * Calculate the traversed distance until we meet an upward edge.
     * Useful for calling on edges between equidistant points.
     * 
     * If we can go up then the distance includes the length of the \param edge
     */
    std::optional<coord_t> distToGoUp() const;

    /*!
     * Check whether the edge goes from a lower to a higher distance_to_boundary.
     * Effectively deals with equidistant edges by looking beyond this edge.
     */
    bool isUpward() const;
    
private:
    int_least8_t is_marked; //! whether the edge is significant; whether the source segments have a sharp angle; -1 is unknown
};




} // namespace arachne
#endif // ARACHNE_SKELETAL_TRAPEZOIDATION_EDGE_H
