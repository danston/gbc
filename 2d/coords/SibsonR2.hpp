// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Sibson coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp
    5. MeshR2.hpp

    This code also depends on the external triangulation library:
    Triangle, https://www.cs.cmu.edu/~quake/triangle.html
    I use the wrapper TriangulatorR2.hpp to include the library in the code.

    Note that this code works with scattered points rather than with polygons. Though you can also treat a polygon
    as a scattered set. In particular, be careful with the boundaries of the convex hull of the given scattered set.
    Along certain lines the code may bug.

*/

#ifndef GBC_SIBSONR2_HPP
#define GBC_SIBSONR2_HPP

// STL includes.
#include <set>
#include <string>
#include <vector>
#include <cassert>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"
#include "../extra/MeshR2.hpp"

// Libs.
#include "../extra/TriangulatorR2.hpp"

namespace gbc {

    // Sibson coordinates for scattered points in R2.
    class SibsonR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        SibsonR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "SibsonR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following paper:
        // N. Sukumar. Voronoi cell finite difference method for the diffusion operator on arbitrary unstructured grids.
        // International Journal for Numerical Methods in Engineering, 57(1):1-34, 2003.
        // However, it is naive and not optimized implementation. So it may be slow for some configurations.
        void compute(const VertexR2 &p, std::vector<double> &b) {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary. If you treat the given set of scattered points as a polygon,
            // better use this function. Otherwise, you can outcomment it.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.

            // Triangulate.
            triangulate();

            // Search for natural neighbours and create the corresponding circumcircles.
            std::vector<VertexR2> nn;
            std::vector<VertexR2> cc;
            std::vector<size_t> t;

            naturalNeighbourSearch(p, nn, cc, t);

            // Compute Sibson's areas.
            const size_t nns = nn.size();
            std::vector<double> A(nns, 0.0);

            computeAreas(p, nn, cc, t, A);

            cc.clear();
            t.clear();

            // Compute coordinates.
            std::vector<double> w(n, 0.0);
            double W = 0.0;

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < nns; ++j) {

                    if (_v[i] == nn[j]) {

                        w[i] = A[j];
                        break;
                    }
                }
                W += w[i];
            }

            nn.clear();

            assert(fabs(W) > 0.0);
            const double invW = 1.0 / W;
            for (size_t i = 0; i < n; ++i) b[i] = w[i] * invW;
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

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Triangle mesh.
        MeshR2 _mesh;

        // Triangulate a given set of points.
        void triangulate() {

            assert(_v.size() > 2);

            TriangulatorR2 tri(_v);

            tri.setPlanarGraph(false);

            tri.allowBoundaryRefinement(false);
            tri.allowEdgeRefinement(false);
            tri.allowSteinerPoints(false);

            std::vector<VertexR2> tp;
            std::vector<Face> tf;

            tri.generate(tp, tf);

            _mesh.clear();
            _mesh.initialize(tp, tf);
            _mesh.createFaces();

            tp.clear();
            tf.clear();
        }

        // Search for natural neighbours in the mesh and compute the corresponding circumcircles.
        // Here we also save the corresponding triangles.
        void naturalNeighbourSearch(const VertexR2 &p,
                                    std::vector<VertexR2> &nn,
                                    std::vector<VertexR2> &cc,
                                    std::vector<size_t> &t) const {

            std::set<size_t> neigh;
            std::set<size_t>::iterator it;
            std::pair<std::set<size_t>::iterator, bool> ret;

            const size_t numF = _mesh.numFaces();

            const std::vector<VertexR2> &v = _mesh.vertices();
            const std::vector<Face> &f = _mesh.faces();

            for (size_t i = 0; i < numF; ++i) {

                const size_t i0 = (size_t) f[i].v[0];
                const size_t i1 = (size_t) f[i].v[1];
                const size_t i2 = (size_t) f[i].v[2];

                const VertexR2 &v0 = v[i0];
                const VertexR2 &v1 = v[i1];
                const VertexR2 &v2 = v[i2];

                VertexR2 c;
                circumCentre(v0, v1, v2, c);

                if (inCircumCircle(p, v2, c)) {

                    neigh.insert(i0);
                    neigh.insert(i1);
                    neigh.insert(i2);

                    t.push_back(i);
                    cc.push_back(c);
                }
            }

            size_t count = 0;

            nn.resize(neigh.size());
            for (it = neigh.begin(); it != neigh.end(); ++it) nn[count++] = v[*it];

            assert(!nn.empty());
        }

        // Compute Sibson's areas. May not work if p is on the edge of the triangle in _mesh!
        // In case it does not work, just slightly move the vertices of the polygon
        // or increase the resolution of the mesh.
        void computeAreas(const VertexR2 &p,
                          const std::vector<VertexR2> &nn,
                          const std::vector<VertexR2> &cc,
                          const std::vector<size_t> &t,
                          std::vector<double> &A) const {

            const size_t tn = t.size();
            assert(tn == cc.size());

            const std::vector<VertexR2> &v = _mesh.vertices();
            const std::vector<Face> &f = _mesh.faces();

            const size_t nns = nn.size();
            assert(nns == A.size());

            for (size_t i = 0; i < tn; ++i) {

                std::vector<VertexR2> lcc(3);
                for (size_t j = 0; j < 3; ++j) {

                    const size_t jp = (j + 1) % 3;
                    const size_t jm = (j + 2) % 3;

                    const VertexR2 &v0 = v[f[t[i]].v[jp]];
                    const VertexR2 &v1 = v[f[t[i]].v[jm]];

                    circumCentre(v0, v1, p, lcc[j]);
                }

                for (size_t j = 0; j < 3; ++j) {

                    const size_t jp = (j + 1) % 3;
                    const size_t jm = (j + 2) % 3;

                    const VertexR2 &tmpV = v[f[t[i]].v[j]];

                    for (size_t k = 0; k < nns; ++k) {
                        if (tmpV == nn[k]) {

                            A[k] += triangleArea(lcc[jp], lcc[jm], cc[i]);
                            break;
                        }
                    }
                }
            }
        }

        // Test if point is inside a circumcircle of a triangle.
        inline bool inCircumCircle(const VertexR2 &p, const VertexR2 &v, const VertexR2 &c) const {
            
            return (p - c).squaredLength() < (v - c).squaredLength();
        }

        // Return the circumcentre of a triangle.
        double circumCentre(const VertexR2 &v0, const VertexR2 &v1, const VertexR2 &v2, VertexR2 &c) const {

            const double D = 2.0 * (v0 - v2).crossProduct(v1 - v2);

            const double v0sl = v0.squaredLength();
            const double v1sl = v1.squaredLength();
            const double v2sl = v2.squaredLength();

            const double a1 = (v0sl - v2sl) * (v1.y() - v2.y());
            const double b1 = (v1sl - v2sl) * (v0.y() - v2.y());

            const double a2 = (v1sl - v2sl) * (v0.x() - v2.x());
            const double b2 = (v0sl - v2sl) * (v1.x() - v2.x());

            assert(fabs(D) > 0.0);

            const double x = (a1 - b1) / D;
            const double y = (a2 - b2) / D;

            c = VertexR2(x, y);

            return D;
        }

        // Return area of the triangle.
        inline double triangleArea(const VertexR2 &v0, const VertexR2 &v1, const VertexR2 &v2) const {

            return 0.5 * ((v1 - v0).crossProduct(v2 - v0));
        }
    };

} // namespace gbc

#endif // GBC_SIBSONR2_HPP
