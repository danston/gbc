// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Barycentric coordinates along a line.

    This class depends on:
    1. VertexExpressionsR2.hpp
    2. VertexR2.hpp

*/

#ifndef GBC_SEGMENTCOORDINATESR2_HPP
#define GBC_SEGMENTCOORDINATESR2_HPP

// STL includes.
#include <utility>
#include <cassert>

// Local includes.
#include "VertexR2.hpp"

namespace gbc {

    // Barycentric coordinates along a line/segment in R2.
    class SegmentCoordinatesR2 {

    public:
        // Constructor.
        SegmentCoordinatesR2(const VertexR2 &v1, const VertexR2 &v2) : _v1(v1), _v2(v2) { }

        // Function that computes coordinates (basically linear interpolation). If the point p is not on the line that
        // passes through the segment, then the coordinates are computed with respect
        // to the projection of p on this line.
        void compute(const VertexR2 &p, std::pair<double, double> &b) const {

            // Project point on the segment and compute the first coordinate.
            const double opposite_scalar_product = (p - _v2).scalarProduct(_v1 - _v2);
            const double b1 = opposite_scalar_product / (_v1 - _v2).squaredLength();

            assert(fabs((_v1 - _v2).squaredLength()) > 0.0);

            // Compute the second coordinate, using the partition of unity property.
            b = std::make_pair(b1, 1.0 - b1);
        }

    private:
        // Vertices of the segment.
        const VertexR2 &_v1;
        const VertexR2 &_v2;
    };

} // namespace gbc

#endif // GBC_SEGMENTCOORDINATESR2_HPP
