// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2107.

// README:
/*

    Harmonic coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp
    5. MeshR2.hpp
    6. Face.hpp

    This code also depends on the external triangulation library:
    Triangle, https://www.cs.cmu.edu/~quake/triangle.html
    I use the wrapper TriangulatorR2.hpp to include the library in the code.

    and it depends on the external linear algebra library: Eigen http://eigen.tuxfamily.org!

*/

#ifndef GBC_HARMONICR2_HPP
#define GBC_HARMONICR2_HPP

// STL includes.
#include <vector>
#include <cassert>

// Local includes.
#include "../extra/Face.hpp"
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"
#include "../extra/MeshR2.hpp"

// Libs.
#include "../extra/libs/eigen/Eigen/Core"
#include "../extra/libs/eigen/Eigen/Dense"
#include "../extra/libs/eigen/Eigen/Sparse"

#include "../extra/TriangulatorR2.hpp"

namespace gbc {

    // Harmonic coordinates in R2.
    class HarmonicR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        HarmonicR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) 
        : super(v, tol), _isMeshCreated(false), _areCoordinatesComputed(false), _isMeshGiven(false) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "HarmonicR2";
        }

        // Function that computes coordinates bb at all points p.
        // Here the set of points p must exclude the polygon's vertices _v.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) {

            // Create internal triangle mesh.
            createMesh(p);

            // Compute coordinates.
            computeCoordinates(bb);
        }

        // Compute coordinates bb at the vertices of the internal mesh with the given edgeLength
        // of the average triangle in this mesh. The number of returned coordinates is the number of the polygon's vertices 
        // + the number of the points p generated by the TriangulatorR2.hpp. All the coordinates are ordered as 
        // _v[0].b(), ..., _v[n-1].b(), p[0].b(), ..., p[N-1].b(), where n = _v.size() and N = p.size().
        void compute(const double edgeLength, std::vector<std::vector<double> > &bb) {

            // Create internal triangle mesh.
            createMesh(edgeLength);

            // Compute coordinates.
            computeCoordinates(bb);
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        // Note that here the function returns a slightly different set of points p from the given one.
        void compute(std::vector<VertexR2> &p) {

            // Compute coordinates.
            std::vector<std::vector<double> > bb;
            compute(p, bb);

            // Wrap vertices with coordinates.
            p = _mesh.vertices();
        }

        // Implementation of the virtual function to compute all coordinates.
        // Here the set of points p must exclude the polygon's vertices _v 
        // apart from the case when the given mesh is used.
        inline void bc(std::vector<VertexR2> &p) {
            
            // Create internal triangle mesh.
            if (!_isMeshGiven) createMesh(p);

            // Compute coordinates.
            std::vector<std::vector<double> > bb;
            computeCoordinates(bb);

            // Wrap vertices with coordinates.
            p = _mesh.vertices();
        }

        // Evaluate coordinates b at any point p.
        void evaluate(const VertexR2 &p, std::vector<double> &b) const {

            assert(_isMeshCreated);
            assert(_areCoordinatesComputed);

            std::vector<double> lambda;
            const int fi = _mesh.findFace(p, lambda);

            assert(lambda.size() == 3);

            const std::vector<VertexR2> &tp = _mesh.vertices();
            const std::vector<Face> &tf = _mesh.faces();

            const int i0 = tf[fi].v[0];
            const int i1 = tf[fi].v[1];
            const int i2 = tf[fi].v[2];

            const VertexR2 &tp0 = tp[i0];
            const VertexR2 &tp1 = tp[i1];
            const VertexR2 &tp2 = tp[i2];

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            for (size_t i = 0; i < n; ++i)
                b[i] = lambda[0] * tp0.b()[i] + lambda[1] * tp1.b()[i] + lambda[2] * tp2.b()[i];
        }

        // Evalaute coordinates at p using the internal storage from the VertexR2 class.
        inline void evaluate(VertexR2 &p) const {
            evaluate(p, p.b());
        }

        // Evaluate coordinates bb at all points p in the vector.
        void evaluate(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) const {

            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) evaluate(p[i], bb[i]);
        }

        // Evaluate coordinates at all points p in the vector using 
        // the internal storage from the VertexR2 class.
        void evaluate(std::vector<VertexR2> &p) const {

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) evaluate(p[i], p[i].b());
        }

        // Set mesh if any.
        void setMesh(const std::vector<VertexR2> &tp, const std::vector<Face> &tf) {

            _mesh.clear();
            _mesh.initialize(tp, tf);

            _isMeshGiven = true;
            _isMeshCreated = false;
        }

        // Clear.
        inline void clear() {
            
            _mesh.clear();

            _isMeshGiven = false;
            _isMeshCreated = false;
            _areCoordinatesComputed = false;
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Eigen typedefs.
        typedef Eigen::MatrixXd MatrixXd;

        // Triangle mesh.
        MeshR2 _mesh;

        // Flags.
        bool _isMeshCreated;
        bool _areCoordinatesComputed;
        bool _isMeshGiven;

        // Given a set of points p, create the internal triangle mesh.
        void createMesh(const std::vector<VertexR2> &p) {

            assert(!p.empty());

            std::vector<VertexR2> tp;
            std::vector<Face> tf;

            TriangulatorR2 tri(_v);

            tri.setPlanarGraph(true);

            tri.allowBoundaryRefinement(false);
            tri.allowEdgeRefinement(false);

            tri.setPoints(p);
            tri.generate(tp, tf);

            _mesh.clear();
            _mesh.initialize(tp, tf);

            _isMeshGiven = false;
            _isMeshCreated = true;
        }

        // Given the edgeLength of the average triangle, create the internal triangle mesh.
        void createMesh(const double edgeLength) {

            assert(edgeLength > 0.0);

            std::vector<VertexR2> tp;
            std::vector<Face> tf;

            TriangulatorR2 tri(_v, edgeLength, true);
            tri.setPlanarGraph(true);

            tri.generate(tp, tf);

            _mesh.clear();
            _mesh.initialize(tp, tf);

            _isMeshGiven = false;
            _isMeshCreated = true;
        }

        // Compute coordinates. This implementation is based on the following paper:
        // P. Joshi, M. Meyer, T. DeRose, B. Green, and T. Sanocki. Harmonic coordinates for character articulation.
        // ACM Transactions on Graphics, 26(3):71:1-9, 2007.
        void computeCoordinates(std::vector<std::vector<double> > &bb) {

            assert(_isMeshCreated || _isMeshGiven);

            // Boundary and interior.
            const size_t n = _v.size();
            const size_t N = _mesh.numVertices();

            assert(N != 0);
            bb.resize(N);

            std::vector<VertexR2> &p = _mesh.vertices();
            std::vector<size_t> indices(N, 0);

            size_t numB = 0, numI = 0;
            for (size_t i = 0; i < N; ++i) {

                if (p[i].type != INTERIOR) indices[i] = numB++;
                else indices[i] = numI++;
            }

            MatrixXd boundary(numB, n);

            Eigen::SparseMatrix<double> A(numI, numI);

            MatrixXd x = MatrixXd::Zero(numI, n);
            MatrixXd b = MatrixXd::Zero(numI, n);

            for (size_t i = 0; i < N; ++i) {

                p[i].b().clear();
                p[i].b().resize(n, 0.0);

                bb[i].clear();
                bb[i].resize(n, 0.0);

                if (p[i].type != INTERIOR) {
                    computeBoundaryCoordinates(p[i]);
                    for (size_t j = 0; j < n; ++j) boundary(indices[i], j) = p[i].b()[j];
                }
            }

            typedef Eigen::Triplet<double> T;
            std::vector<T> tripletList;
            tripletList.reserve(numI * 7);

            for (size_t i = 0; i < N; ++i) {
                if (p[i].type == INTERIOR) {

                    std::vector<int> neighbours;
                    _mesh.getRing(i, neighbours);

                    const size_t nn = neighbours.size();

                    std::vector<double> alphaCot(nn), betaCot(nn);

                    for (size_t j = 0; j < nn; ++j) {

                        const size_t jp = (j + 1) % nn;

                        VertexR2 s1 = p[i] - p[neighbours[j]];
                        VertexR2 s2 = p[neighbours[jp]] - p[neighbours[j]];

                        alphaCot[j] = cotangent(s2, s1);

                        s1 = p[neighbours[j]] - p[neighbours[jp]];
                        s2 = p[i] - p[neighbours[jp]];

                        betaCot[j] = cotangent(s2, s1);
                    }

                    double W = 0.0;
                    for (size_t j = 0; j < nn; ++j) {

                        const size_t jp  = (j + 1) % nn;
                        const size_t idx = neighbours[jp];

                        const double w = -(alphaCot[j] + betaCot[jp]);
                        W -= w;

                        if (p[idx].type != INTERIOR) {
                            for (size_t k = 0; k < n; ++k)
                                b(indices[i], k) -= boundary(indices[idx], k) * w;
                        } else {
                            tripletList.push_back(T(indices[i], indices[idx], w));
                        }
                    }
                    tripletList.push_back(T(indices[i], indices[i], W));
                }
            }

            A.setFromTriplets(tripletList.begin(), tripletList.end());
            A.makeCompressed();
            solveLinearSystem(A, b, x);

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    
                    if (p[j].type != INTERIOR) {
                    
                        p[j].b()[i] = boundary(indices[j], i);
                        bb[j][i] = p[j].b()[i];

                    } else {
                        
                        p[j].b()[i] = x(indices[j], i);
                        bb[j][i] = p[j].b()[i];
                    }
                }
            }

            _areCoordinatesComputed = true;
        }

        // Function that solves the linear system.
        void inline solveLinearSystem(const Eigen::SparseMatrix<double> &A, const MatrixXd &b, MatrixXd &x) const {

            // LDLT solver.
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldlt;

            ldlt.compute(A);
            x = ldlt.solve(b);
        }

        // Compute cotangent.
        double cotangent(const VertexR2 &a, const VertexR2 &b) const {
            return a.scalarProduct(b) / fabs(a.crossProduct(b));
        }
    };

} // namespace gbc

#endif // GBC_HARMONICR2_HPP
