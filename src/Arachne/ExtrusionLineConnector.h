//Copyright (c) 2019 Ultimaker B.V.


#ifndef ARACHNE_EXTRUSION_LINE_CONNECTOR_H
#define ARACHNE_EXTRUSION_LINE_CONNECTOR_H

#include "utils/polygon.h"
#include "utils/ExtrusionSegment.h"
#include "utils/ExtrusionJunction.h"
#include "utils/ExtrusionLine.h"

namespace arachne
{

/*!
 * Connecting ExtrusionLines together into chains / polygons.
 * 
 * If 3 (or more) ExtrusionLines end in the same location,
 * two of the mwill be connected together and the third one will get a bit of its end choppped off a.k.a. 'reduced'.
 * 
 * TODO: there's a bug which connects thw wrong end of a short segment because our algorithm is gready
 * When adding a segment we should explicitly check whether it would be better to connect the segment the other way around.
 */
class ExtrusionLineConnector
{
public:
    /*!
     * \param[in,out] polygons_per_index the resulting polygons by connecting end poitns of polylines together into loops
     * \param[in,out] polylines_per_index The input polylines and the output connected polylines which are not combined into polygons
     * 
     * \param reduce_overlapping_segments Whether to chop off part of the end of an extrusion line if it ends in the same location as other extrusion lines.
     * \param connect_odd_lines_to_polygons Whether to preferably connect an odd count segment to a polygon, rather than closing the polygonal toolpath
     */
    static void optimize(std::vector<std::list<ExtrusionLine>>& polygons_per_index, std::vector<std::list<ExtrusionLine>>& polylines_per_index, bool reduce_overlapping_segments = true, bool connect_odd_lines_to_polygons = true);
private:
    /*!
     * \param polylines_per_index the input polylines which are not yet connected into polygons+polylines. At 3-way intersections no 2 polylins should already be connected.
     * \param intersection_overlap ratio of reduction compared to the toolpath width from where we reduce. At 25% a 0.4mm wide line will be reduced by 0.1mm.
     */
    ExtrusionLineConnector(std::vector<std::list<ExtrusionLine>>& polylines_per_index, float intersection_overlap = 0.25);

    /*!
     * Either end of an \ref ExtrusionLine
     */
    struct ExtrusionLineEndRef
    {
        coord_t inset_idx; //!< The inset index of the original extrusion line
        std::list<ExtrusionLine>::iterator polyline; //!< reference to an extrusionline either in \ref polygons_per_indexor in \ref polylines_per_index
        bool front; //!< Whether this end references the start or the end of the \ref polyline list
        ExtrusionLineEndRef(coord_t inset_idx, std::list<ExtrusionLine>::iterator polyline, bool front)
        : inset_idx(inset_idx)
        , polyline(polyline)
        , front(front)
        {}
        /*!
         * Get the location of the end point
         */
        Point p() const
        {
            return front? polyline->junctions.front().p : polyline->junctions.back().p;
        }
        bool operator==(const ExtrusionLineEndRef& other)
        {
            return inset_idx == other.inset_idx
            && polyline == other.polyline
            && front == other.front;
        }
    };
    
    float intersection_overlap; //!< Amount of reduction as a ratio relative to the extrusion width of the extrusionline which is being reduced
    static constexpr coord_t snap_dist = 10; //!< Distance at which two end points are considered to be at the same location

    std::vector<std::list<ExtrusionLine>>& polylines_per_index; //!< input and output extrusion polylines which aren't connected into polygons

    /*!
     * Connecting polylines together, but allowing for rounding erros in the end points
     * 
     * Reduce unconnected polylines away from the intersection locations as well
     * 
     * \param[out] polygons_per_index output polygons when ends of the same polyline are connected together
     * \param reduce_overlapping_segments whether to chop off ends of extrusion lines if they are overlapping with other extrusion lines
     * \param connect_odd_lines_to_polygons Whether to prefer connecting single line odd segments to normal inset segments over connecting normal inset segments into polygons
     */
    void fuzzyConnect(std::vector<std::list<ExtrusionLine>>& polygons_per_index, bool reduce_overlapping_segments, bool connect_odd_lines_to_polygons);

    /*!
     * Reduce the end of an extrusion line, because it was overlapping with other extrusion line (ends).
     * 
     * Recursive function.
     * 
     * \param polyline the polyline on which to perform the reduction
     * \param polyline_start_it iterator from the end point 'inward'. Either a \ref std::list<ExtrusionJunction>::iterator or a \ref std::list<ExtrusionJunction>::reverse_iterator 
     * \param traveled_dist The distance we have travelled up until this call. (cause this function calls itself recursively)
     * \param reduction_length The total length that this \p polyline should be reduced
     * \param reduction_source Any end point of another polyline with which this \p polyline overlaps
     */
    template<typename directional_iterator>
    void reduceIntersectionOverlap( ExtrusionLine& polyline, directional_iterator polyline_start_it, coord_t traveled_dist, coord_t reduction_length, ExtrusionLineEndRef& reduction_source);

    /*!
     * Get the iterator to the position before which we should insert if we want to insert just 'outward' of \p it 
     * 
     * \param it iterator to element one location more 'inward' than the returned iterator. Either a \ref std::list<ExtrusionJunction>::iterator or a \ref std::list<ExtrusionJunction>::reverse_iterator 
     */
    template<typename directional_iterator>
    static std::list<ExtrusionJunction>::iterator getInsertPosIt(directional_iterator it);


    /*!
     * Get the forward iterator of a forward or backward iterator
     * 
     * \param it Either a \ref std::list<ExtrusionJunction>::iterator or a \ref std::list<ExtrusionJunction>::reverse_iterator 
     */
    template<typename directional_iterator>
    static std::list<ExtrusionJunction>::iterator getSelfPosIt(directional_iterator it);

    /*!
     * Whether the iterator is at the end of its iteration
     */
    template<typename directional_iterator>
    bool isEnd(directional_iterator it, ExtrusionLine& polyline);
    
    /*!
     * Check internal consistency
     */
    void debugCheck();
};




} // namespace arachne
#endif // ARACHNE_EXTRUSION_LINE_CONNECTOR_H
