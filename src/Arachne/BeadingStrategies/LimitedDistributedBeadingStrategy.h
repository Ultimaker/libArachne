//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_BEADING_STRATEGIES_LIMITED_DISTRIBUTED_BEADING_STRATEGY_H
#define ARACHNE_BEADING_STRATEGIES_LIMITED_DISTRIBUTED_BEADING_STRATEGY_H

#include "DistributedBeadingStrategy.h"
#include "utils/logoutput.h"
#include "utils/macros.h"

namespace arachne
{

/*!
 * Beading strategy which evenly subdivides the thickness and tries to stay close to the optimal width.
 */
class LimitedDistributedBeadingStrategy : public DistributedBeadingStrategy
{
public:
    const coord_t max_bead_count;
    LimitedDistributedBeadingStrategy(const coord_t optimal_width, coord_t default_transition_length, const coord_t max_bead_count, float transitioning_angle)
    : DistributedBeadingStrategy(optimal_width, default_transition_length, transitioning_angle)
    , max_bead_count(max_bead_count)
    {
        if (max_bead_count % 2 == 1)
        {
            RUN_ONCE(logWarning("LimitedDistributedBeadingStrategy with odd bead count is odd indeed!\n"));
        }
    }
    virtual ~LimitedDistributedBeadingStrategy() override
    {}
    Beading compute(coord_t thickness, coord_t bead_count) const override;
    coord_t optimalThickness(coord_t bead_count) const override;
    coord_t transitionThickness(coord_t lower_bead_count) const override;
    coord_t optimalBeadCount(coord_t thickness) const override;
    virtual std::string toString() const override { return "LimitedDistributedBeadingStrategy";}
};




} // namespace arachne
#endif // ARACHNE_BEADING_STRATEGIES_LIMITED_DISTRIBUTED_BEADING_STRATEGY_H
