// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Mean value coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef GBC_MEANVALUER2_HPP
#define GBC_MEANVALUER2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Mean value coordinates in R2.
    class MeanValueR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        MeanValueR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "MeanValueR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // M. S. Floater. Generalized barycentric coordinates and applications.
        // Acta Numerica, 24:161-214, 2015 (see Section 5.2).
        // The formula from the paper works only for convex polygons, however here it works for any simple polygon
        // due to the modification proposed in:
        // http://doc.cgal.org/latest/Barycentric_coordinates_2/index.html#Chapter_2D_Generalized_Barycentric_Coordinates (see Section 4.5).
        // Note that this is a more robust but slower O(n^2) version of mean value coordinates. For the fast O(n) version see below.
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

                w[i] = r[im] * r[ip] - s[im].scalarProduct(s[ip]);

                for (size_t j = 0; j < n; ++j) {
                    if (j != im && j != i) {

                        const size_t jp = (j + 1) % n;
                        w[i] *= r[j] * r[jp] + s[j].scalarProduct(s[jp]);
                    }
                }

                w[i] = fabs(w[i]);
                assert(w[i] >= 0.0);

                w[i] = mvsign(A[im], A[i], B[i]) * sqrt(w[i]);

                W += w[i];
            }

            assert(fabs(W) > 0.0);
            const double invW = 1.0 / W;

            for (size_t i = 0; i < n; ++i) b[i] = w[i] * invW;
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following pseudocode:
        // http://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf (see Slide 19).
        void computeFast(const VertexR2 &p, std::vector<double> &b) const {
            
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
            std::vector<double> D(n);
            std::vector<double> t(n);

            for (size_t i = 0; i < n; ++i) {
                
                const size_t ip = (i + 1) % n;

                A[i] = 0.5 * s[i].crossProduct(s[ip]);
                D[i] = s[i].scalarProduct(s[ip]);

                assert(fabs(r[i] * r[ip] + D[i]) > 0.0);

                t[i] = A[i] / (r[i] * r[ip] + D[i]);
            }

            std::vector<double> w(n);
            double W = 0.0;

            for (size_t i = 0; i < n; ++i) {
                
                const size_t im = (i + n - 1) % n;

                assert(fabs(r[i]) > 0.0);

                w[i] = (t[im] + t[i]) / r[i];
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

        // Returns the sign of the mean value weight.
        inline int mvsign(const double Am, const double A, const double B) const {

            if (Am > 0.0 && A > 0.0 && B <= 0.0) return  1;
            if (Am < 0.0 && A < 0.0 && B >= 0.0) return -1;
            
            if (B > 0.0) return  1;
            if (B < 0.0) return -1;

            return 0;
        }
    };

} // namespace gbc

#endif // GBC_MEANVALUER2_HPP
