// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2107.

// README:
/*

    Positive mean value coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef GBC_POSITIVEMEANVALUER2_HPP
#define GBC_POSITIVEMEANVALUER2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Compute a part of the polygon, which is visible from a point p inside this polygon.
    class VisibilityPolygonR2 {

    public:
        // Constructor.
        VisibilityPolygonR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : _v(v), _tol(tol) { }

        // Compute the visibility region vr from a point p with O(n^3) performance. Naive approach that may bug!
        void compute(const VertexR2 &p, std::vector<VertexR2> &vr) {

            vr.clear();

            const size_t n = _v.size();

            VertexR2 res;
            std::vector<std::pair<VertexR2, int> > extra;

            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const size_t ip = (i + 1) % n;

                std::vector<std::pair<VertexR2, int> > tmp;

                bool found1 = false;
                for (size_t j = 0; j < n; ++j) {
                    if (j != im && j != i) {

                        const size_t jp = (j + 1) % n;
                        const bool state = intersectLines(p, _v[i], _v[j], _v[jp], res);
                        assert(state);

                        if (state && isOnSegment(_v[j], _v[jp], res)) {

                            if (isOnSegment(p, _v[i], res)) {
                                
                                found1 = true;
                                continue;

                            } else if (isInFront(res, p, _v[i])) {

                                const int sign = orientation(p, _v[j], _v[jp]);

                                assert(sign != 0);
                                tmp.push_back(std::make_pair(res, sign * j));
                            }
                        }
                    }
                }

                for (size_t j = 0; j < n; ++j) {
                    if (j != i && j != ip) {

                        const size_t jp = (j + 1) % n;
                        const bool state = intersectLines(p, _v[ip], _v[j], _v[jp], res);
                        assert(state);

                        if (state && isOnSegment(_v[j], _v[jp], res)) {
                            if (isOnSegment(p, _v[ip], res)) {

                                break;
                            }
                        }
                    }
                }

                if (!found1) {

                    vr.push_back(_v[i]);

                    if (!tmp.empty()) {
                        const size_t sizeTmp = tmp.size();

                        int closest = -1;
                        double r = std::numeric_limits<double>::max();

                        for (size_t j = 0; j < sizeTmp; ++j) {
                            const double tmpR = (_v[i] - tmp[j].first).squaredLength();

                            if (tmpR < r) {
                                closest = (int) j;
                                r = tmpR;
                            }
                        }

                        assert(closest != -1);
                        extra.push_back(tmp[closest]);
                        tmp.clear();
                    }
                }

                if (!extra.empty()) {

                    std::vector<std::pair<VertexR2, int> > tmpSet;
                    eraseFun(i, extra, tmpSet);

                    if (!tmpSet.empty()) cleanFun(i, tmpSet, vr);
                }
            }

            if (!extra.empty()) {

                for (size_t i = 0; i < extra.size(); ++i) {
                    for (size_t j = 0; j < n; ++j) {

                        if (extra[i].second == (int) j) {

                            bool found = false;
                            const size_t vs = vr.size();

                            for (size_t k = 0; k < vs; ++k) {
                                if (_v[j] == vr[k]) {

                                    vr.insert(vr.begin() + (k + 1) % vs, extra[i].first);

                                    found = true;
                                    break;

                                } else {

                                    if (_v[(j + 1) % n] == vr[k]) {

                                        vr.insert(vr.begin() + k % vs, extra[i].first);

                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if (found) break;
                        }
                    }
                }
                extra.clear();
            }
        }

    private:
        // Vertices of the polygon.
        const std::vector<VertexR2> &_v;

        // Tolerance.
        const double _tol;

        // Erase some not needed vertices.
        void eraseFun(const size_t ind,
                      std::vector<std::pair<VertexR2, int> > &extra,
                      std::vector<std::pair<VertexR2, int> > &tmpSet) {

            const size_t numExtra = extra.size();

            for (size_t j = 0; j < numExtra; ++j) {
                if (extra[j].second == (int) ind) {

                    tmpSet.push_back(extra[j]);
                    extra.erase(extra.begin() + j);

                    eraseFun(ind, extra, tmpSet);

                    return;
                }
            }
        }

        // Clean some temporary sets of vertices.
        void cleanFun(const size_t ind,
                      std::vector<std::pair<VertexR2, int> > &tmpSet,
                      std::vector<VertexR2> &vr) {

            const size_t ts = tmpSet.size();

            double minR = std::numeric_limits<double>::max();
            int foundInd = -1;

            for (size_t j = 0; j < ts; ++j) {
                if ((tmpSet[j].first - _v[ind]).squaredLength() < minR) {
                    foundInd = (int) j;
                }
            }

            assert(foundInd != -1);
            vr.push_back(tmpSet[foundInd].first);
            tmpSet.erase(tmpSet.begin() + foundInd);

            if (!tmpSet.empty()) cleanFun(ind, tmpSet, vr);
            return;
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

        // Check orientation of the query point with respect to the given segment.
        inline int orientation(const VertexR2 &query, const VertexR2 &s1, const VertexR2 &s2) const {

            const size_t turn = checkTurn(s1, s2, query);

            switch (turn) {
                case 0:
                    return -1;
                case 1:
                    return 1;
                case 2:
                    return 0;
                default:
                    return 0;
            }
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

        // Check if a point p is after s2 on the line through the segment [s1, s2] in direction of s2
        // or before s1 on this line.
        inline bool isInFront(const VertexR2 &p, const VertexR2 &s1, const VertexR2 &s2) const {

            std::pair<double, double> b;
            SegmentCoordinatesR2 sc(s1, s2);
            sc.compute(p, b);

            return (b.first <= 0.0 && b.second >= 0.0);
        }
    };

    // Positive mean value coordinates in R2.
    class PositiveMeanValueR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructors.
        PositiveMeanValueR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "PositiveMeanValueR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // Y. Lipman, J. Kopf, D. Cohen-Or, and D. Levin. GPU-assisted positive mean value coordinates for mesh deformations.
        // In Proceedings of SGP 2007, Eurographics Symposium Proceedings, pages 117-123, 2007. A Belyaev and M. Garland editors.
        // Note that this is actually a naive (NOT GPU-based) implementation of the coordinates, which also might be slow and 
        // sometimes produce buggy results.
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.

            // Find a region vr of the polygon that is visible from p.
            std::vector<VertexR2> vr;

            VisibilityPolygonR2 vp(_v);
            vp.compute(p, vr);

            // Compute coordinates with respect to this region.
            computeAll(vr, p, b);
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

        // Compute coordinates with respect too all vertices of the polygon.
        void computeAll(const std::vector<VertexR2> &vr, const VertexR2 &p, std::vector<double> &b) const {

            const size_t n = _v.size();
            const size_t m = vr.size();

            for (size_t i = 0; i < n; ++i) {

                std::vector<double> global(n, 0.0);
                global[i] = 1.0;

                std::vector<double> f(m, 0.0);

                for (size_t j = 0; j < n; ++j) {
                    const size_t jp = (j + 1) % n;

                    for (size_t k = 0; k < m; ++k) {

                        if ((_v[j] - vr[k]).squaredLength() <= _tol) {

                            f[k] = global[j];
                            continue;
                        }

                        if (isOnSegment(_v[j], _v[jp], vr[k]) && (_v[jp] - vr[k]).squaredLength() > _tol) {

                            const double rj  = (_v[j]  - vr[k]).length();
                            const double rjp = (_v[jp] - vr[k]).length();

                            assert(fabs(rj + rjp) > 0.0);

                            f[k] = (rjp * global[j] + rj * global[jp]) / (rj + rjp);
                        }
                    }
                }

                b[i] = computePMVC(vr, p, f);
            }
        }

        // Compute the positive mean value interpolant. This implementation is based on the following paper:
        // M. S. Floater. Generalized barycentric coordinates and applications.
        // Acta Numerica, 24:161-214, 2015 (see Section 5.2).
        // The formula from the paper works only for convex polygons, however here it works for any simple polygon
        // due to the modification proposed in:
        // http://doc.cgal.org/latest/Barycentric_coordinates_2/index.html#Chapter_2D_Generalized_Barycentric_Coordinates (see Section 4.5).
        // Note that this is a more robust but slower O(n^2) version of mean value coordinates.
        double computePMVC(const std::vector<VertexR2> &vr, const VertexR2 &p, const std::vector<double> &f) const {

            const size_t n = vr.size();
            assert(n == f.size());

            std::vector<VertexR2> s(n);
            std::vector<double> r(n);

            for (size_t i = 0; i < n; ++i) {

                s[i] = vr[i] - p;
                r[i] = s[i].length();

                if (fabs(r[i]) < _tol) return f[i];
            }

            std::vector<double> A(n);
            std::vector<double> B(n);
            std::vector<double> D(n);

            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const size_t ip = (i + 1) % n;

                A[i] = 0.5 * s[i].crossProduct(s[ip]);
                D[i] = s[i].scalarProduct(s[ip]);

                if (fabs(A[i]) < _tol && D[i] < 0.0) {
                    assert(fabs(r[i] + r[ip]) > 0.0);

                    return (r[ip] * f[i] + r[i] * f[ip]) / (r[i] + r[ip]);
                }

                B[i] = 0.5 * s[im].crossProduct(s[ip]);
            }

            double F = 0.0, W = 0.0;
            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const size_t ip = (i + 1) % n;

                double w = r[im] * r[ip] - s[im].scalarProduct(s[ip]);

                for (size_t j = 0; j < n; ++j) {
                    if (j != im && j != i) {

                        const size_t jp = (j + 1) % n;
                        w *= r[j] * r[jp] + s[j].scalarProduct(s[jp]);
                    }
                }

                w = fabs(w);
                assert(w >= 0.0);

                w = mvsign(A[im], A[i], B[i]) * sqrt(w);

                F += w * f[i];
                W += w;
            }

            return F / W;
        }

        // Returns the sign of the mean value weight.
        inline int mvsign(const double Am, const double A, const double B) const {

            if (Am > 0.0 && A > 0.0 && B <= 0.0) return  1;
            if (Am < 0.0 && A < 0.0 && B >= 0.0) return -1;

            if (B > 0.0) return  1;
            if (B < 0.0) return -1;

            return 0;
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
    };

} // namespace gbc

#endif // GBC_POSITIVEMEANVALUER2_HPP
