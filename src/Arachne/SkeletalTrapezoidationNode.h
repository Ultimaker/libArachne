//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_SKELETAL_TRAPEZOIDATION_JOINT_H
#define ARACHNE_SKELETAL_TRAPEZOIDATION_JOINT_H

#include "utils/HalfEdge.h"
#include "utils/HalfEdgeNode.h"

#include "utils/IntPoint.h"

namespace arachne
{

class SkeletalTrapezoidationEdge;

class SkeletalTrapezoidationNode : public HalfEdgeNode<SkeletalTrapezoidationEdge, SkeletalTrapezoidationNode>
{
public:
    coord_t distance_to_boundary;
    coord_t bead_count;
    float transition_ratio; //! The distance near the skeleton to leave free because this joint is in the middle of a transition, as a fraction of the inner bead width of the bead at the higher transition.
    SkeletalTrapezoidationNode(Point p)
    : HalfEdgeNode<SkeletalTrapezoidationEdge, SkeletalTrapezoidationNode>(p)
    , distance_to_boundary(-1)
    , bead_count(-1)
    , transition_ratio(0)
    {}

    /*!
     * Whether the node is connected to any marked edge
     */
    bool isMarked() const;
    
    /*!
     * Check whether this node has a locally maximal distance_to_boundary
     * 
     * \param strict Whether equidistant edges can count as a local maximum
     */
    bool isLocalMaximum(bool strict = false) const;

};




} // namespace arachne
#endif // ARACHNE_SKELETAL_TRAPEZOIDATION_JOINT_H
