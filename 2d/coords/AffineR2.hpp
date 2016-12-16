// Author: Dmitry Anisimov, danston@ymail.com.
// Copyright Dmitry Anisimov (c) 2016.

// README:
/*

    This class depends on four other classes that can be found in the extra folder:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

    This class also depends on the external library: Eigen http://eigen.tuxfamily.org!

    Note that this code works with scattered points rather than with polygons. Though you can also treat a polygon
    as a scattered set.

*/

#ifndef GBC_AFFINER2_HPP
#define GBC_AFFINER2_HPP

// STL includes.
#include <vector>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

// Libs.
#include "../extra/libs/eigen/Eigen/Core"
#include "../extra/libs/eigen/Eigen/Dense"

namespace gbc {

    // Affine generalized barycentric coordinates for scattered points in R2.
    class AffineR2 : public BarycentricCoordinatesR2 {

    private:
        // Eigen typedefs.
        typedef Eigen::VectorXd VectorXd;
        typedef Eigen::MatrixXd MatrixXd;

    public:
        // Constructor.
        AffineR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "AffineR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // S. Waldron. Affine generalized barycentric coordinates.
        // Jaen Journal on Approximation, 3(2):209-226, 2011.
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Coordinates.

            // Compute the barycenter.
            VertexR2 c;
            barycentre(c);

            // Set the matrices.
            VertexR2 tmp;
            MatrixXd V(2, n);

            for (size_t i = 0; i < n; ++i) {
                tmp = _v[i] - c;

                V(0, i) = tmp.x();
                V(1, i) = tmp.y();
            }

            const MatrixXd Vs  = V.adjoint();
            const MatrixXd mat = V * Vs;
            const MatrixXd inv = mat.inverse();

            VectorXd vert(2);
            VertexR2 diff;

            for (size_t i = 0; i < n; ++i) {
                diff = p - c;

                vert(0) = V(0, i);
                vert(1) = V(1, i);

                VectorXd res = inv * vert;

                b[i] = diff.x() * res(0) + diff.y() * res(1) + 1.0 / n;
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

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Compute the barycentre of a set of points.
        inline void barycentre(VertexR2 &c) const {

            const size_t n = _v.size();

            for (size_t i = 0; i < n; ++i) c += _v[i];
            c *= 1.0 / n;
        }
    };

} // namespace gbc

#endif // GBC_AFFINER2_HPP
