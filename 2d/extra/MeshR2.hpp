// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Implementation of the triangle/quad mesh in R2.

    This class depends on:
    1. VertexExpressionsR2.hpp
    2. VertexR2.hpp
    3. TriangleCoordinatesR2.hpp
    4. Face.hpp
    5. Halfedge.hpp

*/

#ifndef GBC_MESH_HPP
#define GBC_MESH_HPP

// STL includes.
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <map>

// Local includes.
#include "VertexR2.hpp"
#include "TriangleCoordinatesR2.hpp"
#include "Face.hpp"
#include "Halfedge.hpp"

namespace gbc {

    // Mesh (triangle or quad based) class in R2.
    class MeshR2 {

    public:
        // Constructor.
        MeshR2() : _fn(0), _rval(0), _v(), _he(), _f(), _tol(1.0e-10) { }

        // Initialize mesh.
        void initialize(const std::vector<VertexR2> &v, const std::vector<Face> &f) {

            // Set defaults.
            const int numV = (int) v.size();
            const int numF = (int) f.size();

            assert(numV != 0 && numF != 0);

            assert(f[0].v[0] != -1);
            assert(f[0].v[1] != -1);
            assert(f[0].v[2] != -1);

            _fn = 0;
            for (int i = 0; i < 4; ++i)
                if (f[0].v[i] != -1) ++_fn;

            assert(_fn == 3 || _fn == 4);

            if (_fn == 3) _rval = 6;
            else _rval = 4;

            const int numHE = _fn * numF;

            _v.resize( (size_t)  numV);
            _f.resize( (size_t)  numF);
            _he.resize((size_t) numHE);

            // Copy vertices.
            for (int i = 0; i < numV; ++i) {

                // Assertions.
                assert(v[i].val   == 0);
                assert(v[i].out   == -1);
                assert(v[i].alpha == -1.0);
                assert(v[i].type  == INTERIOR);

                _v[i] = v[i];
            }

            // Copy faces.
            for (int i = 0; i < numF; ++i) {

                // Assertions.
                assert(f[i].f[0] == -1);
                assert(f[i].f[1] == -1);
                assert(f[i].f[2] == -1);
                assert(f[i].f[3] == -1);

                _f[i] = f[i];
            }

            // Copy vertex indices of each triangle in the mesh.
            std::vector<int> vIndices;
            vIndices.resize((size_t) numHE);

            int ind = 0;
            for (int i = 0; i < numF; ++i)
                for (int j = 0; j < _fn; ++j)
                    vIndices[ind++] = f[i].v[j];

            int base, indV = 0, indE = 0;
            for (int i = 0; i < numF; i++) {
                base = indE;
                for (int j = 0; j < _fn; j++) {
                    _he[indE].prev = base + (j + _fn - 1) % _fn;
                    _he[indE].next = base + (j + 1) % _fn;
                    _he[indE++].dest = vIndices[indV++];
                }
            }

            // Build the halfedge connectivity.
            buildHEConnectivity();
        }

        // Get a one ring neighbourhood around the vertex.
        void getRing(const size_t vInd, std::vector<int> &neighbours) const {

            int n = 0, prev, curr, next;

            // Find neighbours of the vertex.
            curr = next = _v[vInd].out;

            assert(_v[vInd].val > 0);
            neighbours.resize((size_t) _v[vInd].val);
            do {
                neighbours[n] = _he[curr].dest;

                prev = _he[curr].prev;
                curr = _he[prev].neigh;
                n++;
            } while ((curr >= 0) && (curr != next));

            // Add one more neighbour for boundary vertices.
            if (curr < 0) {
                curr = _he[prev].prev;
                neighbours[n] = _he[curr].dest;
            }
        }

        // Create faces in the mesh.
        void createFaces() {

            const int numHE = (int) numHalfedges();
            const int numF = numHE / _fn;

            _f.resize((size_t) numF);

            assert(numHE != 0 && numF != 0);

            const int first = _he[0].next;

            // Create faces.
            if (first == 1) {
                int indHE = 0;
                for (int i = 0; i < numF; ++i)
                    for (int j = 0; j < _fn; ++j)
                        _f[i].v[j] = _he[_he[indHE++].prev].dest;
            } else {
                int indHE = 0;
                for (int i = 0; i < numF; ++i) {
                    int ind = _he[indHE].prev;
                    for (int j = 0; j < _fn; ++j) {
                        _f[i].v[j] = _he[ind].dest;
                        ind = _he[ind].next;
                    }
                    indHE++;
                }
            }

            // Create face neighbours.
            int indHE = 0;
            for (int i = 0; i < numF; ++i) {
                for (int j = 0; j < _fn; ++j) {
                    const int neigh = _he[indHE++].neigh;
                    _f[i].f[j] = findFaceFromHE(neigh);
                }
            }
        }

        // Find face that contains the query point.
        inline int findFace(const VertexR2 &query) const {

            std::vector<double> lambda;
            return findFace(query, lambda);
        }

        // Find face that contains the query point and return
        // the corresponding barycentric coordinates.
        int findFace(const VertexR2 &query, std::vector<double> &lambda) const {

            const size_t numF = numFaces();
            for (size_t i = 0; i < numF; ++i) {

                const int i0 = _f[i].v[0];
                const int i1 = _f[i].v[1];
                const int i2 = _f[i].v[2];

                const VertexR2 &v0 = _v[i0];
                const VertexR2 &v1 = _v[i1];
                const VertexR2 &v2 = _v[i2];

                TriangleCoordinatesR2 tc(v0, v1, v2);
                tc.compute(query, lambda);

                if (lambda[0] >= 0.0 && 
                    lambda[1] >= 0.0 && 
                    lambda[2] >= 0.0) return (int) i;
            }

            return -1;
        }

        // Get halfedges that are shared between the face with the index faceInd
        // and its neighbouring faces.
        // The function returns the number of obtained neighbours.
        size_t getNeighbours(const size_t faceInd, std::vector<size_t> &neighs) const {

            assert(numFaces() != 0);

            int heInd = (int) faceInd * 3;
            if (heInd >= 0 && _he[heInd].neigh != -1) neighs.push_back((size_t) heInd);

            heInd = _he[heInd].next;
            if (heInd >= 0 && _he[heInd].neigh != -1) neighs.push_back((size_t) heInd);

            heInd = _he[heInd].next;
            if (heInd >= 0 && _he[heInd].neigh != -1) neighs.push_back((size_t) heInd);

            return neighs.size();
        }

        // Internal data.

        // Return vertices of the mesh.
        inline std::vector<VertexR2> &vertices() {
            return _v;
        }

        // Return const vertices of the mesh.
        inline const std::vector<VertexR2> &vertices() const {
            return _v;
        }

        // Return halfedges of the mesh.
        inline std::vector<Halfedge> &halfedges() {
            return _he;
        }

        // Return const halfedges of the mesh.
        inline const std::vector<Halfedge> &halfedges() const {
            return _he;
        }

        // Return faces of the mesh.
        inline std::vector<Face> &faces() {
            return _f;
        }

        // Return const faces of the mesh.
        inline const std::vector<Face> &faces() const {
            return _f;
        }

        // Return number of vertices in the mesh.
        inline size_t numVertices() const {
            return _v.size();
        }

        // Return number of halfedges in the mesh.
        inline size_t numHalfedges() const {
            return _he.size();
        }

        // Return number of faces in the mesh.
        inline size_t numFaces() const {
            return _f.size();
        }

        // Clear mesh.
        inline void clear() {
            _v.clear();
            _he.clear();
            _f.clear();
        }

        // Check if mesh is empty.
        inline bool isEmpty() const {
            return _v.empty();
        }

        // Tolerance.

        // Tolerance used internally in the class.
        inline double tolerance() const {
            return _tol;
        }

        // Set new user-defined tolerance.
        inline void setTolerance(const double newTol) {
            _tol = newTol;
        }

    private:
        // Number of face vertices.
        int _fn;

        // Regular vertex valency.
        int _rval;

        // Mesh.
        std::vector<VertexR2> _v;  // stores vertices  of the mesh
        std::vector<Halfedge> _he; // stores halfedges of the mesh
        std::vector<Face> _f;      // stores faces     of the mesh

        // Tolerance.
        double _tol;

        // Functions.

        // Build the halfedge connectivity.
        void buildHEConnectivity() {

            const int numV = (int) numVertices();
            const int numHE = (int) numHalfedges();

            typedef std::map<int, int> halfedgeList;
            std::vector<halfedgeList> halfedgeTable((size_t) numV);
            halfedgeList *destList;
            int source, dest;

            // Build the halfedge connectivity.
            for (int i = 0; i < numHE; ++i) {

                // Source and destination of the current halfedge.
                source = _he[_he[i].prev].dest;
                dest = _he[i].dest;

                // Is halfedge from destination to source already in the edge table?
                destList = &(halfedgeTable[dest]);
                std::map<int, int>::iterator it = destList->find(source);

                if (it != destList->end()) {

                    _he[i].neigh = it->second;
                    _he[it->second].neigh = i;
                    destList->erase(it);

                } else {

                    // Put a halfedge in the edge table.
                    halfedgeTable[source].insert(std::make_pair(dest, i));
                }
            }

            // Determine valency and some outgoing halfedge for each vertex.
            // Mark and count the boundary vertices.
            VertexR2 *destV;
            for (int i = 0; i < numHE; ++i) {

                // Destination vertex of the current halfedge.
                destV = &(_v[_he[i].dest]);

                // Increase valency of destination vertex.
                destV->val++;

                // Take next of the current halfedge as the outgoing halfedge.
                destV->out = _he[i].next;

                if (_he[i].neigh < 0) {

                    // NOTE: boundary vertices have one more neighbour than the outgoing edges.
                    destV->type = FLAT; // FLAT vertices, that is collinear
                    destV->val++;
                }
            }

            // Get the "rightmost" outgoing halfedge for boundary vertices.
            for (int i = 0; i < numV; ++i) {
                if (_v[i].type != INTERIOR) {
                    while (_he[_v[i].out].neigh >= 0) {

                        // Move the outgoing halfedge "one to the right".
                        _v[i].out = _he[_he[_v[i].out].neigh].next;
                    }
                }
            }

            // Update type of each boundary vertex in the mesh.
            VertexR2 *prevV;
            VertexR2 *nextV;
            for (int i = 0; i < numHE; ++i) {
                if (_he[i].neigh < 0) {

                    destV = &(_v[_he[i].dest]);

                    const int prev = destV->out;
                    prevV = &(_v[_he[prev].dest]);

                    const int next = _he[i].prev;
                    nextV = &(_v[_he[next].dest]);

                    defineBoundaryVertexType(prevV, destV, nextV);
                }
            }
        }

        // Define type of a vertex on the mesh boundary.
        void defineBoundaryVertexType(VertexR2 *prevV, VertexR2 *destV, VertexR2 *nextV) const {

            const double x1 = prevV->x();
            const double y1 = prevV->y();
            const double x2 = destV->x();
            const double y2 = destV->y();
            const double x3 = nextV->x();
            const double y3 = nextV->y();

            // Compute the following determinant:
            //
            //       | x1 y1 1 |
            // det = | x2 y2 1 |
            //       | x3 y3 1 |
            //
            // det = 0 if three points are collinear;
            // det > 0 if they create a concave corner;
            // det < 0 if they create a convex corner;

            const double det = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);

            // Compute the external signed angle for each corner.
            const double dotProd = (x1 - x2) * (x3 - x2) + (y1 - y2) * (y3 - y2);
            const double extAngle = atan2(det, dotProd); // positive in case of concave and negative in case of convex corner

            const double eps = tolerance();

            if (det > eps) { // CONCAVE corner

                destV->type = CONCAVE;
                destV->alpha = extAngle; // external angle

            } else { // CONVEX corner

                destV->type = CONVEX;
                destV->alpha = 2.0 * M_PI + extAngle; // internal angle
            }
        }

        // Given a halfedge, find a face to which it belongs.
        int findFaceFromHE(const int indHE) const {
            
            const size_t numF = numFaces();

            int tmpHE = 0;
            for (size_t i = 0; i < numF; ++i)
                for (size_t j = 0; j < (size_t) _fn; ++j) {
                    if (indHE == tmpHE) return (int) i;
                    tmpHE++;
                }

            return -1;
        }
    };

} // namespace gbc

#endif // GBC_MESHR2_HPP
