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
    6. MeanValueR2.hpp

*/

#ifndef GBC_MESHR2_HPP
#define GBC_MESHR2_HPP

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
#include "../coords/MeanValueR2.hpp"

namespace gbc {

    // Mesh (triangle or quad-based) class in R2.
    class MeshR2 {

    public:
        // Constructor.
        MeshR2() : _fn(0), _rval(0), _v(), _he(), _f(), _tol(1.0e-10) { }

        // Load mesh from a file.
        void load(const std::string &path) {

            clear();

            if (path.empty()) {
                std::cerr << "ERROR: The mesh file is not provided!\n" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::ifstream readFile(path.c_str(), std::ios_base::in);

            if (!readFile) {
                std::cerr << "ERROR: Error reading an .obj file!\n" << std::endl;
                exit(EXIT_FAILURE);
            }

            std::string flag;

            double x, y, z;
            size_t i1, i2, i3, i4;

            VertexR2 tmpV;
            Face     tmpF;

            // Read an .obj file.
            // Here we always ignore the third "z" coordinate.
            while (!readFile.eof()) {
                readFile >> flag;

                // Vertex flag.
                if (flag == "v") {
                    readFile >> x >> y >> z;

                    tmpV.x() = x;
                    tmpV.y() = y;

                    _v.push_back(tmpV);
                }

                // Face flag.
                if (flag == "f") {
                    readFile >> i1 >> i2 >> i3;

                    assert(i1 > 0 && i2 >0 && i3 > 0);

                    tmpF.v[0] = i1 - 1;
                    tmpF.v[1] = i2 - 1;
                    tmpF.v[2] = i3 - 1;

                    if (readFile.peek() == 32) {
                        readFile >> i4;

                        assert(i4 > 0);
                        tmpF.v[3] = i4 - 1;
                    }

                    _f.push_back(tmpF);
                }
            }

            assert(_v.size() != 0);
            assert(_f.size() != 0);

            readFile.close();

            setMeshFlags(_f);
            initializeHalfedges(_f);
            buildHEConnectivity();
        }

        // Load barycentric coordinates from a file.
        void loadBarycentricCoordinates(const std::string &path) {

            std::ifstream readFile(path.c_str(), std::ios_base::in);

            if (!readFile) {
                std::cerr << "ERROR: Error loading file!\n" << std::endl;
                exit(EXIT_FAILURE);
            }

            size_t numV, numC;

            readFile >> numV;
            readFile >> numC;

            assert(_v.size() == numV);

            for (size_t i = 0; i < numV; ++i) {
                _v[i].b().resize(numC);

                for (size_t j = 0; j < numC; ++j) 
                    readFile >> _v[i].b()[j];
            }

            readFile.close();
        }

        // Initialize mesh.
        void initialize(const std::vector<VertexR2> &v, const std::vector<Face> &f) {

            clear();

            // Set defaults.
            setMeshFlags(f);

            const size_t numV = v.size();
            const size_t numF = f.size();

            assert(numV != 0 && numF != 0);

            _v.resize(numV);
            _f.resize(numF);

            // Copy vertices.
            for (size_t i = 0; i < numV; ++i) {

                // Assertions.
                assert(v[i].val   == 0);
                assert(v[i].out   == -1);
                assert(v[i].alpha == -1.0);
                assert(v[i].type  == INTERIOR);

                _v[i] = v[i];
            }

            // Copy faces.
            for (size_t i = 0; i < numF; ++i) {

                // Assertions.
                assert(f[i].f[0] == -1);
                assert(f[i].f[1] == -1);
                assert(f[i].f[2] == -1);
                assert(f[i].f[3] == -1);

                _f[i] = f[i];
            }

            // Initialize halfedges.
            initializeHalfedges(f);

            // Build the halfedge connectivity.
            buildHEConnectivity();
        }

        // Initialize barycentric coordinates.
        void initializeBarycentricCoordinates(const std::vector< std::vector<double> > &bb)
        {
            const size_t numV = numVertices();
            assert(bb.size() == numV);
            const size_t numC = bb[0].size();
            
            for (size_t i = 0; i < numV; ++i) {

                _v[i].b().resize(numC);
                for (size_t j = 0; j < numC; ++j) _v[i].b()[j] = bb[i][j];
            }
        }

        // Get a one ring neighbourhood around the vertex.
        void getRing(const size_t vInd, std::vector<int> &neighs) const {

            size_t nSize = 0;
            int prev, curr, next;

            // Find neighbours of the vertex.
            curr = next = _v[vInd].out;

            assert(_v[vInd].val > 0);
            neighs.resize((size_t) _v[vInd].val);
            do {
                neighs[nSize] = _he[curr].dest;

                prev = _he[curr].prev;
                curr = _he[prev].neigh;
                nSize++;
            } while ((curr >= 0) && (curr != next));

            // Add one more neighbour for boundary vertices.
            if (curr < 0) {
                curr = _he[prev].prev;
                neighs[nSize] = _he[curr].dest;
            }
        }

        // Create faces in the mesh.
        void createFaces(const bool makeFaceNeighbours = true) {

            const size_t numHE = numHalfedges();
            const size_t numF = numHE / _fn;

            _f.resize(numF);

            assert(numHE != 0 && numF != 0);

            const int first = _he[0].next;

            // Create faces.
            if (first == 1) {
                size_t indHE = 0;
                for (size_t i = 0; i < numF; ++i)
                    for (size_t j = 0; j < _fn; ++j)
                        _f[i].v[j] = _he[_he[indHE++].prev].dest;
            } else {
                size_t indHE = 0;
                for (size_t i = 0; i < numF; ++i) {
                    int ind = _he[indHE].prev;
                    for (size_t j = 0; j < _fn; ++j) {
                        _f[i].v[j] = _he[ind].dest;
                        ind = _he[ind].next;
                    }
                    indHE++;
                }
            }

            // Create face neighbours.
            // Slow code, turn it off if it is not necessary!
            if (makeFaceNeighbours) {

                size_t indHE = 0;
                for (size_t i = 0; i < numF; ++i) {
                    for (size_t j = 0; j < _fn; ++j) {
                        const int neigh = _he[indHE++].neigh;
                        _f[i].f[j] = findFaceFromHE(neigh);
                    }
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

            if (isTriangleMesh()) {

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
            } else if (isQuadMesh()) {

                const size_t numF = numFaces();
                for (size_t i = 0; i < numF; ++i) {
                    
                    const int i0 = _f[i].v[0];
                    const int i1 = _f[i].v[1];
                    const int i2 = _f[i].v[2];
                    const int i3 = _f[i].v[3];

                    std::vector<VertexR2> quad(4);

                    quad[0] = _v[i0];
                    quad[1] = _v[i1];
                    quad[2] = _v[i2];
                    quad[3] = _v[i3];

                    MeanValueR2 mvc(quad);
                    mvc.compute(query, lambda);

                    if (lambda[0] >= 0.0 && 
                        lambda[1] >= 0.0 && 
                        lambda[2] >= 0.0 &&
                        lambda[3] >= 0.0) return (int) i;
                }
            }

            return -1;
        }

        // Get halfedges that are shared between the face with the index faceInd
        // and its neighbouring faces.
        // The function returns the number of obtained neighbours.
        size_t getFaceNeighbours(const size_t faceInd, std::vector<size_t> &neighs) const {

            if (isQuadMesh()) {

                // This functionality is not yet implemented!

                assert(false);
                return 0;
            }

            assert(numFaces() != 0);

            int heInd = (int) faceInd * 3;
            if (heInd >= 0 && _he[heInd].neigh != -1) neighs.push_back((size_t) heInd);

            heInd = _he[heInd].next;
            if (heInd >= 0 && _he[heInd].neigh != -1) neighs.push_back((size_t) heInd);

            heInd = _he[heInd].next;
            if (heInd >= 0 && _he[heInd].neigh != -1) neighs.push_back((size_t) heInd);

            return neighs.size();
        }

        // Subdivision.

        // Initialize mesh with a polygon.
        // In principle, any set of vertices can be used.
        virtual void setFromPolygon(const std::vector<VertexR2>&) {
            // Override me for the corresponding subdivision scheme.
        }

        // Vertex adjustment near the polygon's corners.
        virtual void preprocess(const bool) {
            // Override me for the corresponding subdivision scheme.
        }

        // Subdivide mesh.
        inline void subdivide(const size_t timesToSubdivide, const bool midpoint = false) {
            for (size_t i = 0; i < timesToSubdivide; ++i) 
                subdivideMesh(midpoint);
        }

        // Do one subdivision step.
        virtual void subdivideMesh(const bool) {
            // Override me for the corresponding subdivision scheme.
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

        // Return number of edges in the mesh.
        inline size_t numEdges() const {
            return numVertices() + (numHalfedges() / _fn) - 1;
        }

        // Clear mesh.
        inline void clear() {
            _fn = _rval = 0;

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

        // Type of the mesh.

        // Triangle-based mesh.
        inline bool isTriangleMesh() const {
            return _fn == 3;
        }

        // Quad-based mesh.
        inline bool isQuadMesh() const {
            return _fn == 4;
        }

    protected:
        // Number of face vertices.
        size_t _fn;

        // Regular vertex valency.
        size_t _rval;

        // Mesh.
        std::vector<VertexR2> _v;  // stores vertices  of the mesh
        std::vector<Halfedge> _he; // stores halfedges of the mesh
        std::vector<Face> _f;      // stores faces     of the mesh

        // Tolerance.
        double _tol;

        // Functions.

        // Set mesh flags.
        void setMeshFlags(const std::vector<Face> &f) {

            assert(f[0].v[0] != -1);
            assert(f[0].v[1] != -1);
            assert(f[0].v[2] != -1);

            _fn = 0;
            for (size_t i = 0; i < 4; ++i)
                if (f[0].v[i] != -1) ++_fn;

            assert(_fn == 3 || _fn == 4);

            if (_fn == 3) _rval = 6;
            else _rval = 4;
        }

        // Initialize halfedges.
        void initializeHalfedges(const std::vector<Face> &f) {

            // Copy vertex indices of each triangle in the mesh.
            assert(_fn > 2);

            const size_t numF = f.size();
            const size_t numHE = _fn * numF;

            _he.resize(numHE);

            std::vector<int> vIndices;
            vIndices.resize(numHE);

            size_t ind = 0;
            for (size_t i = 0; i < numF; ++i)
                for (size_t j = 0; j < _fn; ++j)
                    vIndices[ind++] = f[i].v[j];

            size_t base, indV = 0, indE = 0;
            for (size_t i = 0; i < numF; i++) {
                base = indE;

                for (size_t j = 0; j < _fn; j++) {
                    _he[indE].prev = int(base + (j + _fn - 1) % _fn);
                    _he[indE].next = int(base + (j + 1) % _fn);
                    _he[indE++].dest = vIndices[indV++];
                }
            }          
        }

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

        // Set barycentric coordinates with respect to the Lagrange property.
        void setInitialBarycentricCoordinates() {
            
            const size_t numV = numVertices();

            for (size_t i = 0; i < numV; ++i) {
                _v[i].b().resize(numV, 0.0);
                _v[i].b()[i] = 1.0;
            }
        }

    private:
        // Functions.

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
