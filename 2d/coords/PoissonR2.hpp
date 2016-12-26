// Authors: Xian-Ying Li and Shi-Min Hu, http://cg.cs.tsinghua.edu.cn/people/~xianying/
// Updated by: Dmitry Anisimov danston@ymail.com.
// Copyright Xian-Ying Li, Shi-Min Hu, and Dmitry Anisimov (c) 2016-2017.

// README:
/*

    Poisson coordinates and related data structures.

    This code was originally written by Xian-Ying Li and Shi-Min Hu. I modified this code to fit my framework
    with different generalized barycentric coordinates.

*/

#ifndef GBC_POISSONR2_HPP
#define GBC_POISSONR2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Circle structure.
    class Circle {

    public:
        // Constructor.
        Circle(const VertexR2 &c = VertexR2(), const double r = 0.0) : _c(c), _r(r) { };

        // Return center of the circle.
        inline VertexR2 &c() {
            return _c;
        }

        // Return center of the circle (const).
        inline VertexR2 c() const {
            return _c;
        }

        // Return radius of the circle.
        inline double &r() {
            return _r;
        }

        // Return radius of the circle (const).
        inline double r() const {
            return _r;
        }

        // Return x coordinate of the circle's center.
        inline double &cx() {
            return _c.x();
        }

        // Return x coordinate of the circle's center (const).
        inline double cx() const {
            return _c.x();
        }

        // Return y coordinate of the circle's center.
        inline double &cy() {
            return _c.y();
        }

        // Return y coordinate of the circle's center (const).
        inline double cy() const {
            return _c.y();
        }

    private:
        // Circle centre.
        VertexR2 _c;

        // Circle radius.
        double _r;
    };

    // May not always work. If it does not work, slightly modify the polygon!
    // Poisson coordinates in R2.
    class PoissonR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        PoissonR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "PoissonR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // X.-Y. Li and S.-M. Hu. Poisson coordinates.
        // IEEE Transactions on Visualization and Computer Graphics, 19(2):344-352, 2013.
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.

            // Define a circle.
            Circle circle;
            setMinCircle(circle);

            // Compute coordinates.
            interiorCoordinates(p, b, circle);
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

        // Set minimum circle.
        void setMinCircle(Circle &circle) const {

            VertexR2 center = _v[0];
            circle.r() = 0.0;

            for (size_t i = 1; i < _v.size(); ++i) {
                if ((center - _v[i]).length() > circle.r()) {

                    center = _v[i];
                    circle.r() = 0.0;

                    for (size_t j = 1; j <= i - 1; ++j) {
                        if ((center - _v[j]).length() > circle.r()) {

                            center.x() = 0.5 * (_v[i].x() + _v[j].x());
                            center.y() = 0.5 * (_v[i].y() + _v[j].y());

                            circle.r() = 0.5 * (_v[i] - _v[j]).length();

                            for (size_t k = 1; k <= j - 1; ++k) {
                                if ((center - _v[k]).length() > circle.r()) {

                                    intersect(_v[j].x() - _v[i].x(), _v[j].y() - _v[i].y(),
                                              0.5 * (_v[j].x() * _v[j].x() + _v[j].y() * _v[j].y() - _v[i].x() * _v[i].x() -
                                                     _v[i].y() * _v[i].y()),
                                              _v[k].x() - _v[i].x(), _v[k].y() - _v[i].y(),
                                              0.5 * (_v[k].x() * _v[k].x() + _v[k].y() * _v[k].y() - _v[i].x() * _v[i].x() -
                                                     _v[i].y() * _v[i].y()),
                                              center);

                                    circle.r() = (center - _v[k]).length();
                                }
                            }
                        }
                    }
                }
            }
            circle.c() = center;
        }

        // Compute coordinates in the polygon's interior.
        void interiorCoordinates(const VertexR2 &p, std::vector<double> &b, Circle &circle) const {

            const size_t n = _v.size();

            if (fabs(circle.r()) < _tol) {

                circle.c() = circle.c() + p;
                circle.r() = 1.0;

            } else {

                circle.c() = p + (circle.c() - p) * (1.0 / circle.r());
                circle.r() = 1.0;

            }

            std::vector<VertexR2> zeta(n);
            std::vector<VertexR2> xi(n);

            std::vector<double> w(n);

            for (size_t i = 0; i < n; ++i) zeta[i] = _v[i] - p;

            if (fabs(circle.cx() - p.x()) < _tol && fabs(circle.cy() - p.y()) < _tol) {

                for (size_t i = 0; i < n; ++i) xi[i] = zeta[i] * (1.0 / zeta[i].length());

                for (size_t i = 0; i < n; ++i) {
                    
                    const size_t j = (i + 1) % n;

                    const double areaIJ = zeta[i].crossProduct(zeta[j]);
                    const double squaredLengthIJ = (zeta[i] - zeta[j]).squaredLength();

                    if (fabs(areaIJ) < _tol * squaredLengthIJ) {
                        if (crossOrigin(zeta[i], zeta[j])) {

                            boundaryCoordinates(p, b, i, j);
                            return;

                        } else continue;
                    }

                    // Mean Value coordinates.
                    VertexR2 upsilonIJ = rotL(xi[j] - xi[i]);

                    w[i] += zeta[j].crossProduct(upsilonIJ) / areaIJ;
                    w[j] -= zeta[i].crossProduct(upsilonIJ) / areaIJ;
                }

            } else {

                VertexR2 tmpP = p - circle.c();

                const double C = tmpP.squaredLength() - 1.0;

                VertexR2 tau_kappa = tmpP * (1.0 / (C + 1.0));
                VertexR2 tau = tau_kappa + circle.c() - p;

                for (size_t i = 0; i < n; ++i) {

                    const double A = zeta[i].squaredLength();
                    const double B = tmpP.scalarProduct(zeta[i]);

                    assert((B * B - A * C) >= 0.0);

                    xi[i] = zeta[i] * ((-B + sqrt(B * B - A * C)) / A) - tau;
                }

                for (size_t i = 0; i < n; ++i) {

                    const size_t j = (i + 1) % n;

                    const double areaIJ = zeta[i].crossProduct(zeta[j]);
                    const double squaredLengthIJ = (zeta[i] - zeta[j]).squaredLength();

                    if (fabs(areaIJ) < _tol * squaredLengthIJ) {
                        if (crossOrigin(zeta[i], zeta[j])) {

                            boundaryCoordinates(p, b, i, j);
                            return;

                        } else continue;
                    }

                    // Poisson coordinates.
                    VertexR2 logIJ = logg(div(xi[i], xi[j]));
                    VertexR2 upsilonIJ = rotL(mul(tau_kappa, logIJ));

                    w[i] += zeta[j].crossProduct(upsilonIJ) / areaIJ;
                    w[j] -= zeta[i].crossProduct(upsilonIJ) / areaIJ;
                }
            }

            double W = 0.0;
            for (size_t i = 0; i < n; ++i) W += w[i];

            assert(fabs(W) > 0.0);
            const double invW = 1.0 / W;

            for (size_t i = 0; i < n; ++i) b[i] = w[i] * invW;
        }

        // Compute boundary coordinates.
        void boundaryCoordinates(const VertexR2 &p, std::vector<double> &b, const int i, const int j) const {

            const double lengthI = (p - _v[i]).length();
            const double lengthJ = (p - _v[j]).length();

            const double denom = lengthI + lengthJ;

            b[i] = lengthJ / denom;
            b[j] = lengthI / denom;
        }

        // Cross the origin.
        bool crossOrigin(const VertexR2 &a, const VertexR2 &b) const {

            const double areaAB = fabs(b.crossProduct(a));
            const double squaredLengthAB = (a - b).squaredLength();
            const double maxInner = (1.0 + _tol) * squaredLengthAB;

            return areaAB < _tol * squaredLengthAB &&
                   (a - b).scalarProduct(a) < maxInner &&
                   (b - a).scalarProduct(b) < maxInner;
        }

        // Intersection function.
        void intersect(const double a, const double b, const double c, const double d, const double e, const double f,
                       VertexR2 &v) const {
            
            v.x() = (c * e - f * b) / (a * e - b * d);
            v.y() = (c * d - f * a) / (b * d - e * a);
        }

        // Rotation.
        inline VertexR2 rotL(const VertexR2 &a) const {
            return VertexR2(-a.y(), a.x());
        }

        // Multiplication.
        inline VertexR2 mul(const VertexR2 &a, const VertexR2 &b) const {
            return VertexR2(a.x() * b.x() - a.y() * b.y(), a.x() * b.y() + a.y() * b.x());
        }

        // Division.
        inline VertexR2 div(const VertexR2 &a, const VertexR2 &b) const {
            return VertexR2(a.x() * b.x() + a.y() * b.y(), a.y() * b.x() - a.x() * b.y()) * (1.0 / b.squaredLength());
        }

        // Logarithm.
        VertexR2 logg(const VertexR2 &a) const {

            const double R = 0.5 * std::log(a.squaredLength());
            double I = 0.0;

            if (fabs(a.x()) < _tol && a.y() > 0.0) {
                I = 0.5 * M_PI;
            } else if (fabs(a.x()) < _tol && a.y() < 0.0) {
                I = -0.5 * M_PI;
            } else if (a.x() > 0.0) {
                I = atan(a.y() / a.x());
            } else if (a.x() < 0.0 && (a.y() > 0.0 || fabs(a.y()) < _tol)) {
                I = atan(a.y() / a.x()) + M_PI;
            } else if (a.x() < 0.0 && a.y() < 0.0) {
                I = atan(a.y() / a.x()) - M_PI;
            }

            return VertexR2(R, I);
        }
    };

} // namespace gbc

#endif // GBC_POISSONR2_HPP
