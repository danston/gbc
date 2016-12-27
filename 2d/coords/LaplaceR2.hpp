// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Laplace coordinates.

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

#ifndef GBC_LAPLACER2_HPP
#define GBC_LAPLACER2_HPP

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

    // Laplace coordinates for scattered points in R2.
    class LaplaceR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        LaplaceR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol), _mesh() { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "LaplaceR2";
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
            std::vector<VertexR2> tp(_v);
            triangulate(tp);

            // Search for natural neighbours.
            naturalNeighbourSearch(p, tp);

            const size_t k = tp.size();

            tp.push_back(p);
            triangulate(tp);

            // Compute distances.
            std::vector< std::pair<double, VertexR2> > s(k);
            std::vector<double> h(k);

            computeDistances(s, h);

            // Compute coordinates.
            std::vector<double> w(n);
            double W = 0.0;

            for (size_t i = 0; i < n; ++i) {

                bool found = false;
                for (size_t j = 0; j < k; ++j) {
                    if (_v[i] == s[j].second) {

                        assert(fabs(h[j]) > 0.0);

                        w[i] = s[j].first / h[j];
                        found = true;

                        break;
                    }
                }

                if (!found) w[i] = 0.0;
                W += w[i];
            }

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
        void triangulate(std::vector<VertexR2> &tp) {

            assert(tp.size() > 2);

            TriangulatorR2 tri(tp);

            tri.setPlanarGraph(false);

            tri.allowBoundaryRefinement(false);
            tri.allowEdgeRefinement(false);
            tri.allowSteinerPoints(false);

            std::vector<Face> tf;
            tri.generate(tp, tf);

            _mesh.clear();
            _mesh.initialize(tp, tf);
            _mesh.createFaces();

            tp.clear();
        }

        // Search for natural neighbours in the mesh.
        void naturalNeighbourSearch(const VertexR2 &p, std::vector<VertexR2> &tp) const {

            std::set<size_t> nn;
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

                if (inCircumCircle(p, v0, v1, v2)) {
                    nn.insert(i0);
                    nn.insert(i1);
                    nn.insert(i2);
                }
            }

            size_t count = 0;

            tp.resize(nn.size());
            for (it = nn.begin(); it != nn.end(); ++it) tp[count++] = v[*it];
        }

        // Compute lengths of the Voronoi edges and distances to natural neighbours.
        void computeDistances(std::vector< std::pair<double, VertexR2> > &s, std::vector<double> &h) const {

            const size_t k = _mesh.numVertices();

            assert(k > 0);

            const std::vector<VertexR2> &v = _mesh.vertices();

            const VertexR2 &p = v[k - 1];

            std::vector<int> neighs;
            _mesh.getRing(k - 1, neighs);

            const size_t n = neighs.size();

            assert(n == s.size() && n == h.size());

            std::vector<VertexR2> cc(n);

            for (size_t i = 0; i < n; ++i) {
                const size_t ip = (i + 1) % n;

                const VertexR2 &vv = v[neighs[i]];
                const VertexR2 &vp = v[neighs[ip]];

                circumCentre(vv, vp, p, cc[i]);
            }

            for (size_t i = 0; i < n; ++i) {

                const size_t im = (i + n - 1) % n;
                const VertexR2 &vv = v[neighs[i]];

                s[i].first = (cc[i] - cc[im]).length();
                s[i].second = vv;

                h[i] = (vv - p).length();
            }
        }

        // Test if point is inside a circumcircle of a triangle.
        inline bool inCircumCircle(const VertexR2 &p, 
                                   const VertexR2 &v0, const VertexR2 &v1, const VertexR2 &v2) const {
            
            VertexR2 c;
            circumCentre(v0, v1, v2, c);
            return (p - c).squaredLength() < (v2 - c).squaredLength();
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
    };

} // namespace gbc

#endif // GBC_LAPLACER2_HPP
