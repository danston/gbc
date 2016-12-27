// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Barycentric coordinates for any triangle in the plane.

    This class depends on:
    1. VertexExpressionsR2.hpp
    2. VertexR2.hpp

*/

#ifndef GBC_TRIANGLECOORDINATESR2_HPP
#define GBC_TRIANGLECOORDINATESR2_HPP

// STL includes.
#include <vector>
#include <cassert>

// Local includes.
#include "VertexR2.hpp"

namespace gbc {

    // Triangle coordinates in R2.
    class TriangleCoordinatesR2 {

    public:
        // Constructor.
        TriangleCoordinatesR2(const VertexR2 &v1, const VertexR2 &v2, const VertexR2 &v3) : _v1(v1), _v2(v2), _v3(v3) { }

        // Function that computes coordinates. This implementation is based on the following article:
        // https://en.wikipedia.org/wiki/Barycentric_coordinate_system
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();
            b.resize(3);

            // Compute some related sub-areas.
            const double A2 = 0.5 * (_v2 - p).crossProduct(_v3 - p);
            const double A3 = 0.5 * (_v3 - p).crossProduct(_v1 - p);

            // Compute the total area of the triangle.
            const double A = 0.5 * (_v1 - _v2).crossProduct(_v1 - _v3);
            
            // Invert this area.
            assert(fabs(A) > 0.0);
            const double invA = 1.0 / A;

            // Compute the coordinates.
            b[0] = A2 * invA;
            b[1] = A3 * invA;

            // Find the last coordinate using the partition of unity property.
            b[2] = 1.0 - b[0] - b[1];
        }

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p) const {
            compute(p, p.b());
        }

        // Compute coordinates at all points p in the vector.
        void compute(std::vector<VertexR2> &p) const {
            
            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) compute(p[i], p[i].b());
        }

    private:
        // Vertices of the triangle.
        const VertexR2 &_v1;
        const VertexR2 &_v2;
        const VertexR2 &_v3;
    };

} // namespace gbc

#endif // GBC_TRIANGLECOORDINATESR2_HPP
