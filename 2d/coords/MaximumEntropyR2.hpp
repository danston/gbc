// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2107.

// README:
/*

    Maximum entropy coordinates and the corresponding Newton solver.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

    This class also depends on the external library: Eigen http://eigen.tuxfamily.org!

*/

#ifndef GBC_MAXIMUMENTROPYR2_HPP
#define GBC_MAXIMUMENTROPYR2_HPP

// STL includes.
#include <cmath>
#include <vector>
#include <cassert>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

// Libs.
#include "../extra/libs/eigen/Eigen/Core"
#include "../extra/libs/eigen/Eigen/Dense"

namespace gbc {

    // Partition function for maximum entropy coordinates in R2.
    class PartitionFunction {

    private:
        // Eigen typedefs.
        typedef Eigen::VectorXd VectorXd;
        typedef Eigen::MatrixXd MatrixXd;

    public:
        // Overload the operator ().
        inline double operator()(const MatrixXd &vtilde, const VectorXd &m, const VectorXd &lambda, const int index) const {
            
            assert(index >= 0);
            const double dotProduct = lambda(0) * vtilde(index, 0) + lambda(1) * vtilde(index, 1);

            return m(index) * exp(-dotProduct);
        }
    };

    // Solver type used in the Newton solver below.
    // CL - closed form, EG - Eigen.
    enum SolverType {CF, EG};

    // Newton solver for maximum entropy coordinates in R2.
    class NewtonMEC {

    private:
        // Eigen typedefs.
        typedef Eigen::VectorXd VectorXd;
        typedef Eigen::MatrixXd MatrixXd;

    public:
        // Constructor.
        NewtonMEC(const MatrixXd &vtilde,
                  const VectorXd &m,
                  const int n,
                  const int maxNumIter = 1000,
                  const double tol = 1.0e-12)
                : _vtilde(vtilde), _m(m), _n(n), _maxNumIter(maxNumIter), _tol(tol) { 

                    assert(n >= 0);
                    assert(maxNumIter >= 0);
                }

        // Main function.
        inline void solve(VectorXd &lambda, const SolverType solverType = EG) const {
            optimizeParameters(lambda, solverType);
        }

    private:
        // Fields.
        const MatrixXd &_vtilde;
        const VectorXd &_m;
        const int _n;
        const int _maxNumIter;
        const double _tol;

        // Main optimization function.
        void optimizeParameters(VectorXd &lambda, const SolverType solverType) const {

            const double alpha = 1.0;
            for (int k = 0; k < _maxNumIter; k++) {

                VectorXd g(2);
                computeGradient(lambda, g);

                const double gNorm = g.norm();
                if (gNorm < _tol) break;

                MatrixXd H(2, 2);
                computeHessian(lambda, H);

                VectorXd deltaLambda(2);
                solveLinearSystem(g, H, deltaLambda, solverType);

                lambda = lambda + alpha * deltaLambda;
            }
        }

        // Compute first derivative.
        void computeGradient(const VectorXd &lambda, VectorXd &g) const {

            double dZ1 = 0.0, dZ2 = 0.0;
            for (int i = 0; i < _n; ++i) {
                
                PartitionFunction Zi;
                const double Zival = Zi(_vtilde, _m, lambda, i);

                dZ1 += Zival * _vtilde(i, 0);
                dZ2 += Zival * _vtilde(i, 1);
            }

            g(0) = -dZ1;
            g(1) = -dZ2;
        }

        // Compute second derivative.
        void computeHessian(const VectorXd &lambda, MatrixXd &H) const {

            double dZ11 = 0.0, dZ12 = 0.0, dZ22 = 0.0;
            for (int i = 0; i < _n; ++i) {

                PartitionFunction Zi;
                const double Zival = Zi(_vtilde, _m, lambda, i);

                dZ11 += Zival * _vtilde(i, 0) * _vtilde(i, 0);
                dZ12 += Zival * _vtilde(i, 0) * _vtilde(i, 1);
                dZ22 += Zival * _vtilde(i, 1) * _vtilde(i, 1);
            }

            H(0, 0) = dZ11;
            H(0, 1) = H(1, 0) = dZ12;
            H(1, 1) = dZ22;
        }

        // Function that solves linear system.
        void solveLinearSystem(const VectorXd &g,
                               const MatrixXd &H,
                               VectorXd &deltaLambda,
                               const SolverType solverType) const {

            switch (solverType) {
                case CF: {

                    const double denom0 = H(0, 0);
                    const double denom1 = H(1, 1) * H(0, 0) - H(1, 0) * H(0, 1);

                    assert(fabs(denom0) > 0.0 && fabs(denom1) > 0.0);

                    deltaLambda(1) = (H(1, 0) * g(0) - g(1) * H(0, 0)) / denom1;
                    deltaLambda(0) = (-g(0) - H(0, 1) * deltaLambda(1)) / denom0;

                    break;
                }

                case EG:
                    deltaLambda = H.fullPivLu().solve(-g);
                    break;

                default:
                    break;
            };
        }
    };

    // Maximum entropy coordinates in R2.
    class MaximumEntropyR2 : public BarycentricCoordinatesR2 {

    private:
        // Eigen typedefs.
        typedef Eigen::VectorXd VectorXd;
        typedef Eigen::MatrixXd MatrixXd;

    public:
        // Constructor.
        MaximumEntropyR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "MaximumEntropyR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // K. Hormann and N. Sukumar. Maximum entropy coordinates for arbitrary polytopes.
        // Computer Graphics Forum, 27(5):1513-1520, 2008.
        void compute(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            std::vector<VertexR2> s(n);
            Eigen::MatrixXd vtilde(n, 2);
            
            for (size_t i = 0; i < n; ++i) {
                s[i] = _v[i] - p;

                vtilde(i, 0) = s[i].x();
                vtilde(i, 1) = s[i].y();
            }

            VectorXd m(n);
            computePriorFunctions(p, m);

            VectorXd lambda = VectorXd::Zero(2);
            solveOptimizationProblem(vtilde, m, lambda);

            std::vector<double> z(n);
            double Z = 0.0;

            for (size_t i = 0; i < n; ++i) {
                
                PartitionFunction Zi;
                z[i] = Zi(vtilde, m, lambda, (int) i);

                Z += z[i];
            }

            assert(fabs(Z) > 0.0);
            for (size_t i = 0; i < n; ++i) {
                
                b[i] = z[i] / Z;
                assert(!isnan(b[i]));
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

        // Compute all necessary prior functions.
        void computePriorFunctions(const VertexR2 &p, VectorXd &m) const {

            const size_t n = _v.size();

            std::vector<double> r(n);
            std::vector<double> e(n);

            for (size_t i = 0; i < n; ++i) {
                
                const size_t ip = (i + 1) % n;

                r[i] = (_v[i]  -     p).length();
                e[i] = (_v[ip] - _v[i]).length();
            }

            std::vector<double> ro(n);

            for (size_t i = 0; i < n; ++i) {

                const size_t ip = (i + 1) % n;
                ro[i] = r[i] + r[ip] - e[i];
            }

            double PItilde = 0.0;
            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const double denom = ro[im] * ro[i];

                assert(fabs(denom) > 0.0);

                m(i) = 1.0 / denom;
                PItilde += m(i);
            }

            assert(fabs(PItilde) > 0.0);

            for (size_t i = 0; i < n; ++i) m(i) /= PItilde;
        }

        // Function that solves optimization problem.
        void inline solveOptimizationProblem(const MatrixXd &vtilde, const VectorXd &m, VectorXd &lambda) const {

            const size_t n = _v.size();
            NewtonMEC newtonSolver(vtilde, m, (int) n);
            newtonSolver.solve(lambda);
        }
    };

} // namespace gbc

#endif // GBC_MAXIMUMENTROPYR2_HPP
