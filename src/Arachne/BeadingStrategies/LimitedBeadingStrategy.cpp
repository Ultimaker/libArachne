//Copyright (c) 2019 Ultimaker B.V.
#include <cassert>

#include "LimitedBeadingStrategy.h"

namespace arachne
{

LimitedBeadingStrategy::Beading LimitedBeadingStrategy::compute(coord_t thickness, coord_t bead_count) const
{
    if (bead_count <= max_bead_count)
    {
        return parent->compute(thickness, bead_count);
    }
    assert(bead_count == max_bead_count + 1);

    coord_t optimal_thickness = parent->optimalThickness(max_bead_count);
    Beading ret = parent->compute(optimal_thickness, max_bead_count);
    ret.left_over += thickness - ret.total_thickness;
    ret.total_thickness = thickness;
    
    // enforece symmetry
    if (bead_count % 2 == 1)
    {
        ret.toolpath_locations[bead_count / 2] = thickness / 2;
        ret.bead_widths[bead_count / 2] = thickness - optimal_thickness;
    }
    for (coord_t bead_idx = 0; bead_idx < (bead_count + 1) / 2; bead_idx++)
    {
        ret.toolpath_locations[bead_count - 1 - bead_idx] = thickness - ret.toolpath_locations[bead_idx];
    }
    return ret;
}

coord_t LimitedBeadingStrategy::optimalThickness(coord_t bead_count) const
{
    if (bead_count <= max_bead_count)
    {
        return parent->optimalThickness(bead_count);
    }
    return 10000000; // 10 meter
}

coord_t LimitedBeadingStrategy::transitionThickness(coord_t lower_bead_count) const
{
    if (lower_bead_count < max_bead_count)
    {
        return parent->transitionThickness(lower_bead_count);
    }
    if (lower_bead_count == max_bead_count)
    {
        return parent->optimalThickness(lower_bead_count + 1) - 10;
    }
    return 9000000; // 9 meter
}

coord_t LimitedBeadingStrategy::optimalBeadCount(coord_t thickness) const
{
    coord_t parent_bead_count = parent->optimalBeadCount(thickness);
    if (parent_bead_count <= max_bead_count)
    {
        return parent->optimalBeadCount(thickness);
    }
    else if (parent_bead_count == max_bead_count + 1)
    {
        if (thickness < parent->optimalThickness(max_bead_count + 1) - 10)
            return max_bead_count;
        else 
            return max_bead_count + 1;
    }
    else return max_bead_count + 1;
}

} // namespace arachne
