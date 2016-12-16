// Authors: Xian-Ying Li, Tao Ju, and Shi-Min Hu, http://cg.cs.tsinghua.edu.cn/people/~xianying/
// Updated by: Dmitry Anisimov, danston@ymail.com.
// Copyright Xian-Ying Li, Tao Ju, Shi-Min Hu, Dmitry Anisimov (c) 2016.

// README:
/*

    This code was originally written by Xian-Ying Li, Tao Ju, and Shi-Min Hu. I modified this code to fit my framework
    with different generalized barycentric coordinates.

    This class depends on four other classes that can be found in the extra folder:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/

#ifndef GBC_CUBICMEANVALUER2_HPP
#define GBC_CUBICMEANVALUER2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // Cubic mean value coordinates in R2.
    class CubicMeanValueR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        CubicMeanValueR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "CubicMeanValueR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // X.-Y. Li, T. Ju, and S.-M. Hu. Cubic mean value coordinates.
        // ACM Transactions on Graphics, 32(4):126:1-10, 2013.
        // Since cubic mean value have three different sets of coefficients, I save all of them in the internal
        // vectors _a, _b, and _c. However, if you need only the barycentric coordinates associated with the
        // vertices of the polygon and you want to increase the performance, you can outcomment the parts 
        // of the function below, which involve these coefficients.
        void compute(const VertexR2 &p, std::vector<double> &b) {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Clear coefficients.
            _a.clear();
            _b.clear();
            _c.clear();

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) {
             
                // Save all coefficients.
                _a = b;
                _b.resize(2 * n, 0.0);
                _c.resize(2 * n, 0.0);

                return;
            }

            // Interior.
            std::vector<double> gn(2 * n);
            std::vector<double> gt(2 * n);

            // Compute coordinates.
            interiorCoordinates(p, b, gn, gt);

            // Save all coefficients.
            _a = b;
            _b = gt;
            _c = gn;
        }

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p) {
            compute(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) {

            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) compute(p[i], bb[i]);
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p) {

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) compute(p[i], p[i].b());
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

        // Return all coefficients associated with cubic mean value coordinates.
        void coefficients(std::vector<double> &a, 
                          std::vector<double> &b, 
                          std::vector<double> &c) const {
            
            a = _a;
            b = _b;
            c = _c;
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Coefficients.
        std::vector<double> _a;
        std::vector<double> _b;
        std::vector<double> _c;

        // Compute cubic mean value coordinates in the polygon's interior.
        void interiorCoordinates(const VertexR2 &p,
                                 std::vector<double> &b,
                                 std::vector<double> &gn,
                                 std::vector<double> &gt) const {

            const size_t n = _v.size();

            std::vector<VertexR2> y(n);
            std::vector<VertexR2> z(n);

            std::vector<double> bC[3];

            std::vector<double> gnC[3];
            std::vector<double> gtC[3];

            for (size_t i = 0; i < n; ++i) {

                y[i] = _v[i] - p;
                const double yiLength = y[i].length();

                assert(fabs(yiLength) > 0.0);
                if (yiLength > _tol) z[i] = y[i] * (1.0 / yiLength);
            }

            double A[3] = {0.0, 0.0, 0.0};
            double B[3] = {0.0, 0.0, 0.0};
            double C[3] = {0.0, 0.0, 0.0};

            for (size_t i = 0; i < 3; ++i) {

                bC[i].resize(n);

                gnC[i].resize(2 * n);
                gtC[i].resize(2 * n);

                for (size_t j = 0; j < n; ++j) bC[i][j] = 0.0;

                for (size_t j = 0; j < 2 * n; ++j) {
                    
                    gnC[i][j] = 0.0;
                    gtC[i][j] = 0.0;
                }
            }

            for (size_t i = 0; i < n; ++i) {

                const size_t j = (i + 1) % n;

                const double areaIJ = y[i].crossProduct(y[j]);
                const double squaredLengthIJ = (y[i] - y[j]).squaredLength();

                if (std::fabs(areaIJ) < _tol * squaredLengthIJ) {
                    if (crossOrigin(y[i], y[j])) {

                        boundaryCoordinates(p, b, gn, gt, i, j);
                        return;

                    } else continue;
                }

                assert(squaredLengthIJ >= 0.0);
                const double lengthIJ = sqrt(squaredLengthIJ);

                assert(fabs(lengthIJ) > 0.0);
                const double invLengthIJ = 1.0 / lengthIJ;

                assert(fabs(areaIJ) > 0.0);

                VertexR2 alphaI    = y[j] * (1.0 / areaIJ);
                VertexR2 kappaI[2] = {VertexR2(0.5 * alphaI.y(),  0.5 * alphaI.x()),
                                      VertexR2(0.5 * alphaI.y(), -0.5 * alphaI.x())};

                VertexR2 alphaJ    = -y[i] * (1.0 / areaIJ);
                VertexR2 kappaJ[2] = {VertexR2(0.5 * alphaJ.y(),  0.5 * alphaJ.x()),
                                      VertexR2(0.5 * alphaJ.y(), -0.5 * alphaJ.x())};

                VertexR2 kappa[2] = {kappaI[0] + kappaJ[0], kappaI[1] + kappaJ[1]};
                VertexR2 intgIJ0  = conj(z[i]) - conj(z[j]);
                VertexR2 intgIJ1  = z[j] - z[i];
                VertexR2 intgIJ2  = (mul(z[j], mul(z[j], z[j])) - mul(z[i], mul(z[i], z[i]))) * (1.0 / 3.0);
                VertexR2 vecIJ    = (y[j] - y[i]) * invLengthIJ;

                VertexR2 kappaSquare  = mul(kappa[0], kappa[0]);
                VertexR2 kappaIntgIJ1 = mul(kappa[0], intgIJ1);

                VertexR2 kappaConjIntgIJ0 = mul(kappa[1], intgIJ0);
                VertexR2 kappaIJ = mul(kappaI[0], kappaJ[0]);
                VertexR2 kappaI2 = mul(kappaI[0], kappaI[0]);
                VertexR2 kappaI3 = mul(kappaI[0], kappaI2);
                VertexR2 kappaIIntgIJ2 = mul(kappaI[0], intgIJ2);

                const double kappaIAbsSquare = mul(kappaI[0], kappaI[1]).x();

                VertexR2 kappaJ2 = mul(kappaJ[0], kappaJ[0]);
                VertexR2 kappaJ3 = mul(kappaJ[0], kappaJ2);
                VertexR2 kappaJIntgIJ2 = mul(kappaJ[0], intgIJ2);

                const double kappaJAbsSquare = mul(kappaJ[0], kappaJ[1]).x();

                VertexR2 tmpI = mul(kappa[1], kappaI[0]);
                VertexR2 tmpJ = mul(kappa[1], kappaJ[0]);

                double valueI =
                        2.0 * (mul(kappaSquare, kappaIIntgIJ2) + mul(tmpI, kappaIntgIJ1) + 2.0 * tmpI.x() * kappaIntgIJ1).y();
                double valueJ =
                        2.0 * (mul(kappaSquare, kappaJIntgIJ2) + mul(tmpJ, kappaIntgIJ1) + 2.0 * tmpJ.x() * kappaIntgIJ1).y();

                bC[0][i] += 2.0 * valueI;
                bC[0][j] += 2.0 * valueJ;

                A[0] += 2.0 * (valueI + valueJ);

                VertexR2 intgI = rotR(
                        mul(kappa[0], kappaIIntgIJ2) + 2.0 * tmpI.x() * intgIJ1 + mul(kappaI[1], kappaConjIntgIJ0));
                VertexR2 intgJ = rotR(
                        mul(kappa[0], kappaJIntgIJ2) + 2.0 * tmpJ.x() * intgIJ1 + mul(kappaJ[1], kappaConjIntgIJ0));

                bC[1][i] += 3.0 * intgI.x();
                bC[1][j] += 3.0 * intgJ.x();
                bC[2][i] += 3.0 * intgI.y();
                bC[2][j] += 3.0 * intgJ.y();

                B[0] += (intgI.x() + intgJ.x());
                C[0] += (intgI.y() + intgJ.y());

                A[1] += 3.0 * (intgI.x() + intgJ.x());
                A[2] += 3.0 * (intgI.y() + intgJ.y());

                double gIJ = -vecIJ.x() * (intgI.x() + intgJ.x()) - vecIJ.y() * (intgI.y() + intgJ.y());

                gnC[0][2 * i]     += (vecIJ.y() * intgI.x() - vecIJ.x() * intgI.y());
                gnC[0][2 * i + 1] += (vecIJ.y() * intgJ.x() - vecIJ.x() * intgJ.y());

                bC[0][i] -= gIJ * invLengthIJ;
                bC[0][j] += gIJ * invLengthIJ;

                intgI = rotR(kappaIIntgIJ2 + mul(kappaI[1], intgIJ1));
                intgJ = rotR(kappaJIntgIJ2 + mul(kappaJ[1], intgIJ1));

                valueI = 2.0 * mul(kappaI[0], intgIJ1).y();
                valueJ = 2.0 * mul(kappaJ[0], intgIJ1).y();

                double coscosI = 0.5 * (valueI + intgI.x());
                double sinsinI = 0.5 * (valueI - intgI.x());
                double sincosI = 0.5 * intgI.y();

                double coscosJ = 0.5 * (valueJ + intgJ.x());
                double sinsinJ = 0.5 * (valueJ - intgJ.x());
                double sincosJ = 0.5 * intgJ.y();

                B[1] += 2.0 * (coscosI + coscosJ);
                B[2] += 2.0 * (sincosI + sincosJ);
                C[1] += 2.0 * (sincosI + sincosJ);
                C[2] += 2.0 * (sinsinI + sinsinJ);

                gIJ = -vecIJ.x() * (coscosI + coscosJ) - vecIJ.y() * (sincosI + sincosJ);

                gnC[1][2 * i]     += (vecIJ.y() * coscosI - vecIJ.x() * sincosI);
                gnC[1][2 * i + 1] += (vecIJ.y() * coscosJ - vecIJ.x() * sincosJ);

                bC[1][i] -= gIJ * invLengthIJ;
                bC[1][j] += gIJ * invLengthIJ;

                gIJ = -vecIJ.x() * (sincosI + sincosJ) - vecIJ.y() * (sinsinI + sinsinJ);

                gnC[2][2 * i]     += (vecIJ.y() * sincosI - vecIJ.x() * sinsinI);
                gnC[2][2 * i + 1] += (vecIJ.y() * sincosJ - vecIJ.x() * sinsinJ);

                bC[2][i] -= gIJ * invLengthIJ;
                bC[2][j] += gIJ * invLengthIJ;

                tmpI = mul(kappaI[1], kappaJ[0]);
                tmpJ = mul(kappaJ[1], kappaI[0]);

                valueI = 2.0 * lengthIJ *
                         (mul(kappaI2, kappaJIntgIJ2) + mul(mul(tmpI, kappaI[0]) + 2.0 * tmpI.x() * kappaI[0], intgIJ1)).y();
                valueJ = 2.0 * lengthIJ *
                         (mul(kappaJ2, kappaIIntgIJ2) + mul(mul(tmpJ, kappaJ[0]) + 2.0 * tmpJ.x() * kappaJ[0], intgIJ1)).y();

                gtC[0][2 * i]     += 2.0 * valueI;
                gtC[0][2 * i + 1] += 2.0 * valueJ;

                VertexR2 tmpIntgII = mul(kappaI2, intgIJ2) + 2.0 * kappaIAbsSquare * intgIJ1 + mul(conj(kappaI2), intgIJ0);
                VertexR2 tmpIntgJJ = mul(kappaJ2, intgIJ2) + 2.0 * kappaJAbsSquare * intgIJ1 + mul(conj(kappaJ2), intgIJ0);
                VertexR2 tmpIntgIJ = mul(kappaIJ, intgIJ2) + (tmpI + tmpJ).x() * intgIJ1 + mul(conj(kappaIJ), intgIJ0);

                intgI = rotR(tmpIntgII - 2.0 * tmpIntgIJ);
                intgJ = rotR(tmpIntgJJ - 2.0 * tmpIntgIJ);

                gtC[0][2 * i]     -= (vecIJ.x() * intgI.x() + vecIJ.y() * intgI.y());
                gtC[0][2 * i + 1] += (vecIJ.x() * intgJ.x() + vecIJ.y() * intgJ.y());

                const double kappaAbs = kappa[0].length();

                assert(fabs(kappaAbs) > 0.0);
                const double invKappaAbs = 1.0 / kappaAbs;

                VertexR2 invKappa = kappa[1] * invKappaAbs * invKappaAbs;
                VertexR2 kappaR   = mul(kappa[1], invKappa);

                VertexR2 intgZ1 = (atangent(mul(z[j], kappa[0]) * invKappaAbs) -
                                   atangent(mul(z[i], kappa[0]) * invKappaAbs)) * invKappaAbs;
                VertexR2 intgZ0 = mul(intgIJ1, invKappa) - mul(intgZ1, kappaR);

                VertexR2 sI = mul(kappaI3, invKappa);
                VertexR2 tI = mul(kappaI2, 3.0 * kappaI[1] - mul(kappaI[0], kappaR));

                VertexR2 sJ = mul(kappaJ3, invKappa);
                VertexR2 tJ = mul(kappaJ2, 3.0 * kappaJ[1] - mul(kappaJ[0], kappaR));

                intgI = lengthIJ * rotR(tmpIntgII -
                                        (mul(sI, intgIJ2) + mul(conj(sI), intgIJ0) + mul(tI, intgZ0) + mul(conj(tI), intgZ1)));
                intgJ = lengthIJ * rotR(tmpIntgJJ -
                                        (mul(sJ, intgIJ2) + mul(conj(sJ), intgIJ0) + mul(tJ, intgZ0) + mul(conj(tJ), intgZ1)));

                gtC[1][2 * i]     += 3.0 * intgI.x();
                gtC[1][2 * i + 1] += 3.0 * intgJ.x();
                gtC[2][2 * i]     += 3.0 * intgI.y();
                gtC[2][2 * i + 1] += 3.0 * intgJ.y();

                sI = 3.0 * mul(kappaI2, invKappa) - 2.0 * kappaI[0];
                double uI = 6.0 * kappaIAbsSquare - 6.0 * mul(kappaI2, kappaR).x();

                sJ = 3.0 * mul(kappaJ2, invKappa) - 2.0 * kappaJ[0];
                double uJ = 6.0 * kappaJAbsSquare - 6.0 * mul(kappaJ2, kappaR).x();

                intgI = rotR(mul(sI, intgIJ2) + mul(conj(sI), intgIJ1) + uI * intgZ0);
                intgJ = rotR(mul(sJ, intgIJ2) + mul(conj(sJ), intgIJ1) + uJ * intgZ0);

                valueI = 2.0 * mul(sI, intgIJ1).y() + 2.0 * uI * intgZ1.y();
                valueJ = 2.0 * mul(sJ, intgIJ1).y() + 2.0 * uJ * intgZ1.y();

                coscosI = 0.5 * (valueI + intgI.x());
                sinsinI = 0.5 * (valueI - intgI.x());
                sincosI = 0.5 * intgI.y();
                coscosJ = 0.5 * (valueJ + intgJ.x());
                sinsinJ = 0.5 * (valueJ - intgJ.x());
                sincosJ = 0.5 * intgJ.y();

                gtC[1][2 * i]     -= (vecIJ.x() * coscosI + vecIJ.y() * sincosI);
                gtC[1][2 * i + 1] += (vecIJ.x() * coscosJ + vecIJ.y() * sincosJ);
                gtC[2][2 * i]     -= (vecIJ.x() * sincosI + vecIJ.y() * sinsinI);
                gtC[2][2 * i + 1] += (vecIJ.x() * sincosJ + vecIJ.y() * sinsinJ);

                for (size_t k = 0; k < 3; ++k) {

                    bC[k][i] += (gtC[k][2 * i] - gtC[k][2 * i + 1]) * invLengthIJ;
                    bC[k][j] -= (gtC[k][2 * i] - gtC[k][2 * i + 1]) * invLengthIJ;
                }
            }

            double lambda[3] = {B[1] * C[2] - B[2] * C[1], B[2] * C[0] - B[0] * C[2], B[0] * C[1] - B[1] * C[0]};
            const double sum = A[0] * lambda[0] + A[1] * lambda[1] + A[2] * lambda[2];
            
            assert(fabs(sum) > 0.0);
            const double invSum = 1.0 / sum;

            lambda[0] *= invSum;
            lambda[1] *= invSum;
            lambda[2] *= invSum;

            for (size_t i = 0; i < n; ++i) 
                b[i] = lambda[0] * bC[0][i] + lambda[1] * bC[1][i] + lambda[2] * bC[2][i];

            for (size_t j = 0; j < 2 * n; ++j) {

                gn[j] = lambda[0] * gnC[0][j] + lambda[1] * gnC[1][j] + lambda[2] * gnC[2][j];
                gt[j] = lambda[0] * gtC[0][j] + lambda[1] * gtC[1][j] + lambda[2] * gtC[2][j];
            }
        }

        // Compute boundary coordinates.
        void boundaryCoordinates(const VertexR2 &p,
                                 std::vector<double> &b,
                                 std::vector<double> &gn,
                                 std::vector<double> &gt,
                                 const int i,
                                 const int j) const {

            const double lengthI = (p - _v[i]).length();
            const double lengthJ = (p - _v[j]).length();

            const double lengthIJ = (_v[i] - _v[j]).length();

            const double denom = lengthI + lengthJ;

            const double alphaI = lengthJ / denom;
            const double alphaJ = lengthI / denom;

            const double cubicI = alphaI * alphaI * alphaJ;
            const double cubicJ = alphaJ * alphaJ * alphaI;

            b[i] = alphaI + cubicI - cubicJ;
            b[j] = alphaJ + cubicJ - cubicI;

            gn[2 * i]     = 0.0;
            gn[2 * i + 1] = 0.0;

            gt[2 * i]     = lengthIJ * cubicI;
            gt[2 * i + 1] = lengthIJ * cubicJ;
        }

        // Cross the origin.
        bool crossOrigin(const VertexR2 &a, const VertexR2 &b) const {

            const double areaAB = std::fabs(b.crossProduct(a));
            const double squaredLengthAB = (a - b).squaredLength();
            const double maxInner = (1.0 + _tol) * squaredLengthAB;

            return areaAB < _tol * squaredLengthAB &&
                   (a - b).scalarProduct(a) < maxInner &&
                   (b - a).scalarProduct(b) < maxInner;
        }

        // Multiplication.
        inline VertexR2 mul(const VertexR2 &a, const VertexR2 &b) const {
            return VertexR2(a.x() * b.x() - a.y() * b.y(), a.x() * b.y() + a.y() * b.x());
        }

        // Division.
        inline VertexR2 div(const VertexR2 &a, const VertexR2 &b) const {
            return VertexR2(a.x() * b.x() + a.y() * b.y(), a.y() * b.x() - a.x() * b.y()) * (1.0 / b.squaredLength());
        }

        // Left rotation.
        inline VertexR2 rotL(const VertexR2 &a) const {
            return VertexR2(-a.y(), a.x());
        }

        // Right rotation.
        inline VertexR2 rotR(const VertexR2 &a) const {
            return VertexR2(a.y(), -a.x());
        }

        // Conjugation.
        inline VertexR2 conj(const VertexR2 &a) const {
            return VertexR2(a.x(), -a.y());
        }

        // Atangent.
        inline VertexR2 atangent(const VertexR2 &a) const {
            return 0.5 * rotR(logg(div(VertexR2(1.0, 0.0) + rotL(a), VertexR2(1.0, 0.0) - rotL(a))));
        }

        // Logarithm.
        VertexR2 logg(const VertexR2 &a) const {

            const double R = 0.5 * log(a.squaredLength());
            double I = 0.0;

            if (std::fabs(a.x()) < _tol && a.y() > 0.0) {
                I = 0.5 * M_PI;
            } else if (std::fabs(a.x()) < _tol && a.y() < 0.0) {
                I = 3.0 * 0.5 * M_PI;
            } else if (a.x() > 0.0 && a.y() >= 0.0) {
                I = atan(a.y() / a.x());
            } else if (a.x() > 0.0 && a.y() < 0.0) {
                I = atan(a.y() / a.x()) + 2.0 * M_PI;
            } else if (a.x() < 0.0) {
                I = atan(a.y() / a.x()) + M_PI;
            }

            return VertexR2(R, I);
        }
    };

} // namespace gbc

#endif // GBC_CUBICMEANVALUER2_HPP
