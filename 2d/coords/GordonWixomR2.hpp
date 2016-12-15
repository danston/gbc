// Author: Dmitry Anisimov, danston@ymail.com.
// Copyright Dmitry Anisimov (c) 2016.

// README:
/*

    This class depends on four other classes that can be found in the extra folder:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

    This class also depends on the external integration library: gauss_legendre.hpp!
    http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/

    Please not that this is not the closed-form implementation of Gordon-Wixom coordinates but rather
    a direct numerical evaluation of the underlying integral, so the code may be slow and the final quality
    of the coordinates depends on the number of samples used in the integration. 
    The default number of samples is 1024.

*/

#ifndef GBC_GORDONWIXOMR2_HPP
#define GBC_GORDONWIXOMR2_HPP

// STL includes.
#include <vector>
#include <cassert>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

// Libs.
#include "../extra/libs/gaussinteg/gauss_legendre.hpp"

namespace gbc {

    // Gordon-Wixom coordinates in R2.
    class GordonWixomR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        GordonWixomR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol), _numSamples(1024) { }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // A. Belyaev. On transfinite barycentric coordinates.
        // In Proceedings of the Fourth Eurographics Symposium on Geometry Processing, pages 89-99, 2006.
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            computeAll(p, b);
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

        // Set number of samples used in the numerical integration.
        inline void setIntegralPrecision(const size_t numSamples) {
            _numSamples = numSamples;
        }

        // Get number of samples currently used in the numerical integration.
        inline size_t numSamples() const {
            return _numSamples;
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Integral precision.
        size_t _numSamples;

        // Compute coordinates with respect too all vertices of the polygon.
        void computeAll(const VertexR2 &p, std::vector<double> &b) const {

            const size_t n = _v.size();

            for (size_t i = 0; i < n; ++i) {

                std::vector<double> f(n, 0.0);
                f[i] = 1.0;

                b[i] = computeGWC(p, f);
            }
        }

        // Compute the Gordon-Wixom interpolant.
        inline double computeGWC(const VertexR2 &p, const std::vector<double> &f) const {

            const double a = 0.0;
            const double b = 2.0 * M_PI;

            const double value = estimateIntegral(p, f, a, b);

            return value / (2.0 * M_PI);
        }

        // Estimate the Gordon-Wixom integral.
        double estimateIntegral(const VertexR2 &p, const std::vector<double> &f,
                                const double a, const double b) const {

            double *xx = NULL;
            double *ww = NULL;
            double A, B, Axx, s;

            const size_t m = (_numSamples + 1) >> 1;

            bool allocated = false;
            for (size_t i = 0; i < GLAWSIZE; ++i) {
                if (glaw[i].n == (int) _numSamples) {

                    xx = glaw[i].x;
                    ww = glaw[i].w;

                    break;
                }
            }

            if (xx == NULL) {
                allocated = true;

                xx = (double *) malloc(m * sizeof(double));
                ww = (double *) malloc(m * sizeof(double));

                gauss_legendre_tbl(int(_numSamples), xx, ww, _tol);
            }

            A = 0.5 * (b - a);
            B = 0.5 * (b + a);

            if (_numSamples & 1) { /* _numSamples - odd */

                s = ww[0] * estimateFunction(p, f, B);
                for (size_t i = 1; i < m; ++i) {

                    Axx = A * xx[i];
                    s += ww[i] * (estimateFunction(p, f, B + Axx) + estimateFunction(p, f, B - Axx));
                }

            } else { /* _numSamples - even */

                s = 0.0;
                for (size_t i = 0; i < m; ++i) {

                    Axx = A * xx[i];
                    s += ww[i] * (estimateFunction(p, f, B + Axx) + estimateFunction(p, f, B - Axx));
                }
            }

            if (allocated) {
                free(xx);
                free(ww);
            }

            return A * s;
        }

        // Estimate the function under the integral.
        double estimateFunction(const VertexR2 &p, const std::vector<double> &f, const double angle) const {

            const double r = 1000.0;

            VertexR2 infP( r * cos(angle),  r * sin(angle));
            VertexR2 infM(-r * cos(angle), -r * sin(angle));

            const size_t n = _v.size();

            assert(f.size() == n);

            VertexR2 res1, res2;

            std::vector<std::pair<VertexR2, int> > y1;
            std::vector<std::pair<VertexR2, int> > y2;

            std::vector<double> ff1, ff2;
            std::vector<double> ro1, ro2;

            for (size_t i = 0; i < n; ++i) {
                const size_t ip = (i + 1) % n;

                if (intersectLines(p, infP, _v[i], _v[ip], res1)) {

                    if (isOnSegment(p, infP, res1) && isOnSegment(_v[i], _v[ip], res1)) {

                        const size_t turn = checkTurn(_v[i], _v[ip], p);
                        int sign = 0;

                        switch(turn) {
                            case 0:
                                sign = -1;
                                break;
                            case 1:
                                sign = 1;
                                break;
                            case 2:
                                sign = 0;
                                break;
                            default:
                                break;
                        }

                        y1.push_back(std::make_pair(res1, sign));

                        const double ri  = (_v[i]  - res1).length();
                        const double rip = (_v[ip] - res1).length();

                        assert(fabs(ri + rip) > 0.0);

                        ff1.push_back((rip * f[i] + ri * f[ip]) / (ri + rip));
                        ro1.push_back((res1 - p).length());
                    }
                }

                if (intersectLines(p, infM, _v[i], _v[ip], res2)) {

                    if (isOnSegment(p, infM, res2) && isOnSegment(_v[i], _v[ip], res2)) {

                        const size_t turn = checkTurn(_v[i], _v[ip], p);
                        int sign = 0;

                        switch(turn) {
                            case 0:
                                sign = -1;
                                break;
                            case 1:
                                sign = 1;
                                break;
                            case 2:
                                sign = 0;
                                break;
                            default:
                                break;
                        }

                        y2.push_back(std::make_pair(res2, sign));

                        const double ri  = (_v[i]  - res2).length();
                        const double rip = (_v[ip] - res2).length();

                        assert(fabs(ri + rip) > 0.0);

                        ff2.push_back((rip * f[i] + ri * f[ip]) / (ri + rip));
                        ro2.push_back((res2 - p).length());
                    }
                }
            }

            const size_t m1 = ff1.size();
            const size_t m2 = ff2.size();

            assert(m1 > 0);
            assert(m2 > 0);

            double sumF1 = 0.0, sumF2 = 0.0;
            for (size_t i = 0; i < m1; ++i) {

                assert(fabs(ro1[i]) > 0.0);
                const double val = double(y1[i].second) / ro1[i];

                sumF1 += ff1[i] * val;
                sumF2 += val;
            }

            for (size_t i = 0; i < m2; ++i) {

                assert(fabs(ro2[i]) > 0.0);
                const double val = double(y2[i].second) / ro2[i];

                sumF1 += ff2[i] * val;
                sumF2 += val;
            }

            assert(fabs(sumF2) > 0.0);
            const double sumF = sumF1 / sumF2;

            return sumF;
        }

        // Intersect two lines.
        inline bool intersectLines(const VertexR2 &a1, const VertexR2 &a2,
                                   const VertexR2 &b1, const VertexR2 &b2,
                                   VertexR2 &query) const {

            bool state = true;
            const double lambda = intersectionCoefficient(a1, a2, b1, b2, state);
            if (!state) return false;

            query = a1 + (a2 - a1) * lambda;
            return true;
        }

        // Intersection coefficient for two lines.
        inline double intersectionCoefficient(const VertexR2 &a1, const VertexR2 &a2,
                                              const VertexR2 &b1, const VertexR2 &b2, bool &state) const {
            
            const VertexR2 u = a2 - a1;
            const VertexR2 v = b2 - b1;
            const VertexR2 w = b2 - a1;

            const double denom = u[0] * v[1] - u[1] * v[0];

            if (fabs(denom) < _tol) state = false;
            else state = true;

            return -(-w[0] * v[1] + w[1] * v[0]) / denom;
        }

        // Check if the point is on the segment.
        inline bool isOnSegment(const VertexR2 &v1, const VertexR2 &v2, const VertexR2 &query) const {

            const double length = (v2 - v1).length();
            const double actualDist = pointLineDistance(v1, v2, query);

            const double dist1 = (v1 - query).length();
            const double dist2 = (v2 - query).length();

            return actualDist < _tol && dist1 <= length && dist2 <= length;
        }

        // The distance from a point to a line.
        inline double pointLineDistance(const VertexR2 &v1, const VertexR2 &v2, const VertexR2 &query) const {

            VertexR2 dir = v2 - v1;

            assert(fabs(dir.length()) > 0.0);
            dir *= 1.0 / dir.length();

            return (v1 - query - (v1 - query).scalarProduct(dir) * dir).length();
        }

        // Check turns. The turn is checked with respect to the last vertex v3.
        inline size_t checkTurn(const VertexR2 &v1, const VertexR2 &v2, const VertexR2 &v3) const {

            const VertexR2 a = v3 - v2;
            const VertexR2 b = v1 - v2;

            const double area = a.crossProduct(b);

            size_t result = 2; // collinear

            if (area < 0.0) result = 0; // right turn
            if (area > 0.0) result = 1; // left turn

            return result;
        }
    };

} // namespace gbc

#endif // GBC_GORDONWIXOMR2_HPP
