// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2107.

// README:
/*

    Metric coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef GBC_METRICR2_HPP
#define GBC_METRICR2_HPP

// STL includes.
#include <vector>
#include <cassert>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Metric coordinates in R2.
    class MetricR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        MetricR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "MetricR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // K. Hormann and M. S. Floater. Mean value coordinates for arbitrary planar polygons.
        // ACM Transactions on Graphics, 25(4):1424-1441, 2006 (see Section 3).
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            std::vector<VertexR2> s(n);
            std::vector<VertexR2> e(n);

            std::vector<double> r(n);

            for (size_t i = 0; i < n; ++i) {
                
                const size_t ip = (i + 1) % n;

                s[i] = _v[i]  - p;
                e[i] = _v[ip] - _v[i];

                r[i] = s[i].length();
            }

            std::vector<double> A(n);
            std::vector<double> B(n);
            std::vector<double> C(n);
            std::vector<double> q(n);

            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const size_t ip = (i + 1) % n;

                A[i] = 0.5 * s[i].crossProduct(s[ip]);
                B[i] = 0.5 * s[im].crossProduct(s[ip]);
                C[i] = 0.5 * (-e[im].crossProduct(e[i]));

                q[i] = r[i] + r[ip] - e[i].length();
            }

            std::vector<double> w(n);
            double W = 0.0;

            for (size_t i = 0; i < n; ++i) {

                const size_t imm = (i + n - 2) % n;
                const size_t im  = (i + n - 1) % n;
                const size_t ip  = (i + 1) % n;

                assert(fabs(q[imm] * q[im] * C[im]) > 0.0);

                assert(fabs(q[im] * q[i] * C[i]) > 0.0);

                assert(fabs(q[i] * q[ip] * C[ip]) > 0.0);

                w[i] = A[imm] / (q[imm] * q[im] * C[im]) - B[i] / (q[im] * q[i] * C[i]) + A[ip] / (q[i] * q[ip] * C[ip]);
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

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;
    };

} // namespace gbc

#endif // GBC_METRICR2_HPP
