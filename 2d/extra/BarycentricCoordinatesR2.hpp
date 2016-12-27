// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    This is the base class for all generalized barycentric coordinates.

    This class depends on:
    1. SegmentCoordinatesR2.hpp
    2. VertexExpressionsR2.hpp
    3. VertexR2.hpp

*/

#ifndef GBC_BARYCENTRICCOORDINATESR2_HPP
#define GBC_BARYCENTRICCOORDINATESR2_HPP

// STL includes.
#include <vector>
#include <utility>

// Local includes.
#include "VertexR2.hpp"
#include "SegmentCoordinatesR2.hpp"

namespace gbc {

    // Generalized barycentric coordinates in R2.
    class BarycentricCoordinatesR2 {

    public:
        // Constructor.
        BarycentricCoordinatesR2(const std::vector<VertexR2> &v, const double tol) 
        : _v(v), _tol(tol) {
            
            assert(!_v.empty());
        }

        // This is a virtual function that is used to compute all coordinate classes at the same time.
        virtual void bc(std::vector<VertexR2> &p) = 0; 

        // Return name of the current coordinate function.
        virtual std::string name() const = 0;

    protected:
        // Vertices of the polygon.
        const std::vector<VertexR2> &_v;

        // Tolerance.
        const double _tol;

        // Compute boundary coordinates (linearly interpolate along the polygon's edges).
        bool computeBoundaryCoordinates(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();

            std::vector<VertexR2> s(n);
            std::vector<double> r(n);

            // Vertices.
            b.resize(n, 0.0);

            for (size_t i = 0; i < n; ++i) {
                
                s[i] = _v[i] - p;
                r[i] = s[i].length();

                if (fabs(r[i]) < _tol) {
                    b[i] = 1.0;
                    return true;
                }
            }

            std::vector<double> A(n);
            std::vector<double> D(n);

            // Edges.
            for (size_t i = 0; i < n; ++i) {
                
                size_t ip = (i + 1) % n;

                A[i] = 0.5 * s[i].crossProduct(s[ip]);
                D[i] = s[i].scalarProduct(s[ip]);

                if (fabs(A[i]) < _tol && D[i] < 0.0) {

                    SegmentCoordinatesR2 segcoords(_v[i], _v[ip]);
                    std::pair<double, double> sc;

                    segcoords.compute(p, sc);
                    b[i]  = sc.first;
                    b[ip] = sc.second;

                    return true;
                }
            }
            return false;
        }

        // Compute boundary coordinates at p using the internal storage from the VertexR2 class.
        inline void computeBoundaryCoordinates(VertexR2 &p) const {
            computeBoundaryCoordinates(p, p.b());
        }
    };

} // namespace gbc

#endif // GBC_BARYCENTRICCOORDINATESR2_HPP
