//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_BEADING_STRATEGIES_INWARD_DISTRIBUTED_BEADING_STRATEGY_H
#define ARACHNE_BEADING_STRATEGIES_INWARD_DISTRIBUTED_BEADING_STRATEGY_H

#include "DistributedBeadingStrategy.h"

namespace arachne
{

/*!
 * Beading strategy which divides the discrepancy between the current thickness and optimal thickness mainly to the inner beads.
 * This causes the outer inset to be constant width almost everywhere.
 */
class InwardDistributedBeadingStrategy : public DistributedBeadingStrategy
{
public:
    /*!
     * \param distribution_radius the radius (in number of beads) over which to distribute the discrepancy between the feature size and the optimal thickness
     */
    InwardDistributedBeadingStrategy(const coord_t optimal_width, coord_t default_transition_length, float transitioning_angle, float distribution_radius)
    : DistributedBeadingStrategy(optimal_width, default_transition_length, transitioning_angle)
    , one_over_distribution_radius_squared(1.0f / distribution_radius * 1.0f / distribution_radius)
    {
    }
    virtual ~InwardDistributedBeadingStrategy() override
    {}
    Beading compute(coord_t thickness, coord_t bead_count) const override;
    virtual std::string toString() const override { return "InwardDistributedBeadingStrategy";}
private:
    float one_over_distribution_radius_squared; // (1 / distribution_radius)^2
};




} // namespace arachne
#endif // ARACHNE_BEADING_STRATEGIES_INWARD_DISTRIBUTED_BEADING_STRATEGY_H
