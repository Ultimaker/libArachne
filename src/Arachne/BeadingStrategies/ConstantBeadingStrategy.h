//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_BEADING_STRATEGIES_CONSTANT_BEADING_STRATEGY_H
#define ARACHNE_BEADING_STRATEGIES_CONSTANT_BEADING_STRATEGY_H

#include "BeadingStrategy.h"

namespace arachne
{

/*!
 * Beading strategy which evenly subdivides the thickness and tries to stay close to the optimal width.
 */
class ConstantBeadingStrategy : public BeadingStrategy
{
public:
    const coord_t bead_count;
    ConstantBeadingStrategy(const coord_t bead_width, const coord_t bead_count, float transitioning_angle)
    : BeadingStrategy(bead_width, /*default_transition_length=(should remain unused)*/-1, transitioning_angle)
    , bead_count(bead_count)
    {
    }
    virtual ~ConstantBeadingStrategy() override
    {}
    Beading compute(coord_t thickness, coord_t bead_count) const override;
    coord_t optimalThickness(coord_t bead_count) const override;
    coord_t transitionThickness(coord_t lower_bead_count) const override;
    coord_t optimalBeadCount(coord_t thickness) const override;
    virtual std::string toString() const override { return "ConstantBeadingStrategy";}
};




} // namespace arachne
#endif // ARACHNE_BEADING_STRATEGIES_CONSTANT_BEADING_STRATEGY_H
