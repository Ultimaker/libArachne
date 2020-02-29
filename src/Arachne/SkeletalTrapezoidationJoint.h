//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_SKELETAL_TRAPEZOIDATION_JOINT_H
#define ARACHNE_SKELETAL_TRAPEZOIDATION_JOINT_H

#include "utils/HalfEdge.h"
#include "utils/HalfEdgeNode.h"

#include "utils/IntPoint.h"

namespace arachne
{

class SkeletalTrapezoidationJoint : public HalfEdgeNode<SkeletalTrapezoidationEdge, SkeletalTrapezoidationJoint>
{
public:
    coord_t distance_to_boundary;
    coord_t bead_count;
    float transition_ratio; //! The distance near the skeleton to leave free because this joint is in the middle of a transition, as a fraction of the inner bead width of the bead at the higher transition.
    SkeletalTrapezoidationJoint(Point p)
    : HalfEdgeNode<SkeletalTrapezoidationEdge, SkeletalTrapezoidationJoint>(p)
    , distance_to_boundary(-1)
    , bead_count(-1)
    , transition_ratio(0)
    {}
};




} // namespace arachne
#endif // ARACHNE_SKELETAL_TRAPEZOIDATION_JOINT_H
