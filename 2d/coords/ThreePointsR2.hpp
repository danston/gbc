// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2107.

// README:
/*

    Three-point coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef GBC_THREEPOINTSR2_HPP
#define GBC_THREEPOINTSR2_HPP

// STL includes.
#include <vector>
#include <cassert>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Three-point coordinates in R2.
    class ThreePointsR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        ThreePointsR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol), _p(1.0) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "ThreePointsR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // M. S. Floater, K. Hormann, and G. Kos. A general construction of barycentric coordinates over convex polygons.
        // Advances in Computational Mathematics, 24(1-4):311-331, 2006.
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            std::vector<VertexR2> s(n);
            std::vector<double> r(n);

            for (size_t i = 0; i < n; ++i) {
                s[i] = _v[i] - p;
                r[i] = s[i].length();
            }

            std::vector<double> A(n);
            std::vector<double> B(n);

            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const size_t ip = (i + 1) % n;

                A[i] = 0.5 * s[i].crossProduct(s[ip]);
                B[i] = 0.5 * s[im].crossProduct(s[ip]);
            }

            std::vector<double> w(n);
            double W = 0.0;

            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const size_t ip = (i + 1) % n;

                assert(fabs(A[im] * A[i]) > 0.0);

                w[i] = (std::pow(r[ip], _p) * A[im] - std::pow(r[i], _p) * B[i] + std::pow(r[im], _p) * A[i]) / (A[im] * A[i]);
                W += w[i];
            }

            assert(fabs(W) > 0.0);
            const double invW = 1.0 / W;

            for (size_t i = 0; i < n; ++i) b[i] = w[i] * invW;
        }

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p) const {
            compute(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) const {

            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) compute(p[i], bb[i]);
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p) const {

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) compute(p[i], p[i].b());
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

        // Set power used in three-point coordinates.
        inline void setPower(const double p) {
            _p = p;
        }

        // Get currently used value of the power.
        inline double power() const {
            return _p;
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Power. Default power is 1.0 that is mean value coordinates.
        double _p;
    };

} // namespace gbc

#endif // GBC_THREEPOINTSR2_HPP
