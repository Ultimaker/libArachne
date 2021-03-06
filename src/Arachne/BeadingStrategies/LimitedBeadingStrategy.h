//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_BEADING_STRATEGIES_LIMITED_BEADING_STRATEGY_H
#define ARACHNE_BEADING_STRATEGIES_LIMITED_BEADING_STRATEGY_H

#include "BeadingStrategy.h"
#include "utils/logoutput.h"
#include "utils/macros.h"

namespace arachne
{

/*!
 * Beading strategy which evenly subdivides the thickness and tries to stay close to the optimal width.
 */
class LimitedBeadingStrategy : public BeadingStrategy
{
public:
    const coord_t max_bead_count;
    const BeadingStrategy* parent;
    LimitedBeadingStrategy(const coord_t max_bead_count, BeadingStrategy* parent)
    : BeadingStrategy(parent->optimal_width, /*default_transition_length=*/-1, parent->transitioning_angle)
    , max_bead_count(max_bead_count)
    , parent(parent)
    {
        if (max_bead_count % 2 == 1)
        {
            RUN_ONCE(logWarning("LimitedBeadingStrategy with odd bead count is odd indeed!\n"));
        }
    }
    virtual ~LimitedBeadingStrategy() override
    {
        delete parent;
    }
    Beading compute(coord_t thickness, coord_t bead_count) const override;
    coord_t optimalThickness(coord_t bead_count) const override;
    coord_t transitionThickness(coord_t lower_bead_count) const override;
    coord_t optimalBeadCount(coord_t thickness) const override;
    virtual std::string toString() const override { return std::string("LimitedBeadingStrategy+") + parent->toString();}
    coord_t getTransitioningLength(coord_t lower_bead_count) const override
    {
        return parent->getTransitioningLength(lower_bead_count);
    }
    float getTransitionAnchorPos(coord_t lower_bead_count) const override
    {
        return parent->getTransitionAnchorPos(lower_bead_count);
    }
};




} // namespace arachne
#endif // LIMITED_DISTRIBUTED_BEADING_STRATEGY_H
