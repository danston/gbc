// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Wachspress coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef GBC_WACHSPRESSR2_HPP
#define GBC_WACHSPRESSR2_HPP

// STL includes.
#include <vector>
#include <cassert>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Wachspress coordinates in R2.
    class WachspressR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        WachspressR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "WachspressR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // M. S. Floater, K. Hormann and G. Kos. A general construction of barycentric coordinates over convex polygons.
        // Advances in Computational Mathematics 24(1-4):311-331, 2006.
        // The formula is modified as in:
        // http://doc.cgal.org/latest/Barycentric_coordinates_2/index.html#Chapter_2D_Generalized_Barycentric_Coordinates (see Section 4.3).
        // Note that this is a more robust but slower O(n^2) version of Wachspress coordinates. For the fast O(n) version see below.
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            std::vector<VertexR2> s(n);
            std::vector<VertexR2> e(n);

            for (size_t i = 0; i < n; ++i) {
                const size_t ip = (i + 1) % n;

                s[i] = _v[i] - p;
                e[i] = _v[ip] - _v[i];
            }

            std::vector<double> A(n, 1.0);
            std::vector<double> C(n);

            for (size_t i = 0; i < n; ++i) {
                const size_t im = (i + n - 1) % n;

                for (size_t j = 0; j < n; ++j) {
                    if (j != im && j != i) {

                        const size_t jp = (j + 1) % n;
                        A[i] *= 0.5 * s[j].crossProduct(s[jp]);
                    }
                }

                C[i] = 0.5 * e[i].crossProduct(-e[im]);
            }

            std::vector<double> w(n);
            double W = 0.0;

            for (size_t i = 0; i < n; ++i) {

                w[i] = C[i] * A[i];
                W += w[i];
            }

            assert(fabs(W) > 0.0);
            const double invW = 1.0 / W;

            for (size_t i = 0; i < n; ++i) b[i] = w[i] * invW;
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // M. Meyer, H. Lee, A. Barr, and M. Desbrun. Generalized barycentric coordinates on irregular polygons.
        // Journal of Graphics Tools, 7(1):13-22, 2002.
        void computeFast(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            std::vector<double> w(n);
            double W = 0.0;

            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const size_t ip = (i + 1) % n;

                w[i] = (cotangent(p - _v[i], _v[im] - _v[i]) + cotangent(p - _v[i], _v[ip] - _v[i])) /
                       (_v[i] - p).squaredLength();

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

        inline void computeFast(VertexR2 &p) const {
            computeFast(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) const {

            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) compute(p[i], bb[i]);
        }

        void computeFast(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) const {

            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) computeFast(p[i], bb[i]);
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p) const {

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) compute(p[i], p[i].b());
        }

        void computeFast(std::vector<VertexR2> &p) const {

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) computeFast(p[i], p[i].b());
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Compute cotangent.
        double cotangent(const VertexR2 &a, const VertexR2 &b) const {
            return a.scalarProduct(b) / fabs(a.crossProduct(b));
        }
    };

} // namespace gbc

#endif // GBC_WACHSPRESSR2_HPP
