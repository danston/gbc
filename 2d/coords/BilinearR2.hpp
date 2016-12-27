// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Bilinear coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef GBC_BILINEARR2_HPP
#define GBC_BILINEARR2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Bilinear coordinates in R2.
    class BilinearR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        BilinearR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) {
            
            assert(v.size() == 4);
        }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "BilinearR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // M. S. Floater. Generalized barycentric coordinates and applications.
        // Acta Numerica, 24:161-214, 2015 (see Section 3).
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();

            assert(n == 4);
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            std::vector<VertexR2> s(n);

            std::vector<double> A(n);
            std::vector<double> B(n);
            std::vector<double> E(n);

            for (size_t i = 0; i < n; ++i) 
                s[i] = _v[i] - p;

            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const size_t ip = (i + 1) % n;

                A[i] = 0.5 * s[i].crossProduct(s[ip]);
                B[i] = 0.5 * s[im].crossProduct(s[ip]);
            }

            const double D = B[0] * B[0] + B[1] * B[1] + 2.0 * (A[0] * A[2] + A[1] * A[3]);

            for (size_t i = 0; i < n; ++i) {

                const size_t ip = (i + 1) % n;

                E[i] = 2.0 * A[i] - B[i] - B[ip] + std::sqrt(D);
            }

            for (size_t i = 0; i < n; ++i) {

                const size_t ip  = (i + 1) % n;
                const size_t ipp = (i + 2) % n;

                assert(fabs(E[ip] * E[ipp]) > 0.0);

                b[i] = (4.0 * A[ip] * A[ipp]) / (E[ip] * E[ipp]);
            }
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

        // Evaluate bilinear functions x * y, (1 - x) * y, (1 - x) * (1 - y), and x * (1 - y) at p on a unit square.
        void evaluateOnUnitSquare(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            assert(n == 4);

            b.resize(n, 0.0);

            assert(_v[0] == VertexR2(0.0, 0.0));
            assert(_v[1] == VertexR2(1.0, 0.0));
            assert(_v[2] == VertexR2(1.0, 1.0));
            assert(_v[3] == VertexR2(0.0, 1.0));

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            b[0] = p.x() * p.y();
            b[1] = (1.0 - p.x()) * p.y();
            b[2] = (1.0 - p.x()) * (1.0 - p.y());
            b[3] = p.x() * (1.0 - p.y());
        }

        // Evaluate bilinear functions bb on a unit square at all points p in the vector.
        void evaluateOnUnitSquare(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) const {

            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) evaluateOnUnitSquare(p[i], bb[i]);
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;
    };

} // namespace gbc

#endif // GBC_BILINEARR2_HPP
