//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_BEADING_STRATEGIES_WIDENING_BEADING_STRATEGY_H
#define ARACHNE_BEADING_STRATEGIES_WIDENING_BEADING_STRATEGY_H

#include "BeadingStrategy.h"

namespace arachne
{

/*!
 * Beading strategy which evenly subdivides the thickness and tries to stay close to the optimal width.
 */
class WideningBeadingStrategy : public BeadingStrategy
{
public:
    BeadingStrategy* parent;
    const coord_t min_input_width;
    const coord_t min_output_width;
    /*!
     * Takes responsibility for deleting \param parent
     */
    WideningBeadingStrategy(BeadingStrategy* parent, const coord_t min_input_width, const coord_t min_output_width)
    : BeadingStrategy(parent->optimal_width, /*default_transition_length=*/-1, parent->transitioning_angle)
    , parent(parent)
    , min_input_width(min_input_width)
    , min_output_width(min_output_width)
    {
    }
    virtual ~WideningBeadingStrategy() override
    {
        if (parent) delete parent;
    }
    virtual Beading compute(coord_t thickness, coord_t bead_count) const override;
    virtual coord_t optimalThickness(coord_t bead_count) const override;
    virtual coord_t transitionThickness(coord_t lower_bead_count) const override;
    virtual coord_t optimalBeadCount(coord_t thickness) const override;
    virtual coord_t getTransitioningLength(coord_t lower_bead_count) const override;
    virtual float getTransitionAnchorPos(coord_t lower_bead_count) const override;
    virtual std::vector<coord_t> getNonlinearThicknesses(coord_t lower_bead_count) const override;
    virtual std::string toString() const override { return std::string("Widening+") + parent->toString();}
};

} // namespace arachne
#endif // ARACHNE_BEADING_STRATEGIES_WIDENING_BEADING_STRATEGY_H
