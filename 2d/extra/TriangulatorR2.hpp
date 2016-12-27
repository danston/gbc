// Copyright Teseo Schneider and Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Wrapper for the Triangle library to compute Delaunay triangulation of the given polygon.

    This class depends on:
    1. VertexExpressionsR2.hpp
    2. VertexR2.hpp
    3. Face.hpp

    This code also depends on the external triangulation library:
    Triangle, https://www.cs.cmu.edu/~quake/triangle.html

*/

#ifndef GBC_TRIANGULATOR2_HPP
#define GBC_TRIANGULATOR2_HPP

#ifdef SINGLE
#define REAL float
#else
#define REAL double
#endif

static const double SQRT3_OVER_4 = 0.43301270189222;

// STL includes.
#include <vector>
#include <sstream>

// Local includes.
#include "VertexR2.hpp"
#include "Face.hpp"

// Libs.
#include "libs/triangle/triangle.hpp"

namespace gbc {

    // Triangulator class in R2.
    class TriangulatorR2 {

    private:
        // Triangle typedefs.
        typedef struct tl::triangulateio Triangulateio;

    public:
        // Set default flags.
        inline void allowBoundaryRefinement(const bool value) {
            _allowBoundaryRefinement = value;
        }

        inline void allowEdgeRefinement(const bool value) {
            _allowEdgeRefinement = value;
        }

        inline void constrainAngle(const bool value) {
            _constrainAngle = value;
        }

        inline void setPlanarGraph(const bool value) {
            _planarGraph = value;
        }

        inline void recordOrder(const bool value) {
            _recordOrder = value;
        }

        inline void allowSteinerPoints(const bool value) {
            _allowSteinerPoints = value;
        }

        inline void setBoundaries(const std::vector<VertexR2> &v) {
            _v = v;
        }

        inline void setPoints(const std::vector<VertexR2> &p) {
            _p = p;
        }

        inline void setEdgeTag(const bool value) {
            _edgeTag = value;
        }

        inline void setMaxEdgeLength(const double value) {
            if (value > 0.0) _maxArea = value * value * SQRT3_OVER_4;
            else _maxArea = -1.0;
        }

        // Constructors.
        TriangulatorR2() {
            resetFlags();
        }

        TriangulatorR2(const std::vector<VertexR2> &v,
                       const double maxEdgeLength = -1.0,
                       const bool refineBoundary = false) {

            resetFlags();
            setBoundaries(v);

            setMaxEdgeLength(maxEdgeLength);
            allowBoundaryRefinement(refineBoundary);
        }

        // Generate triangulation.
        void generate(std::vector<VertexR2> &v, std::vector<Face> &f) {

            Triangulateio out;
            Triangulateio in;

            assignValues(in);
            triangulate(in, out);

            freeIn(in);

            setTriangulation(out, v, f);
            freeOut(out);

            if (!_edgeTag) return;
        }

    private:
        // Boundary vertices.
        std::vector<VertexR2> _v;

        // Extra points.
        std::vector<VertexR2> _p;

        // Flags.
        bool _allowBoundaryRefinement, _allowEdgeRefinement;
        bool _constrainAngle, _planarGraph;

        double _maxArea;
        int _allowSteinerPoints;
        bool _edgeTag, _recordOrder, _verbose;

        // Reset all flags.
        void resetFlags() {

            _allowBoundaryRefinement = true;
            _allowEdgeRefinement = true;

            _constrainAngle = true;
            _planarGraph = false;

            _maxArea = -1.0;

            _allowSteinerPoints = -1;

            _edgeTag = false;
            _recordOrder = false;
            _verbose = false;
        }

        // Check if point is on the boundary.
        bool inline isOnBoundary(const VertexR2 &p) const {

            const size_t n = _v.size();
            for (size_t i = 0; i < n; ++i)
                if (isOnSegment(i, p)) return true;
            return false;
        }

        // Check if point is on the boundary segment with index i.
        bool inline isOnSegment(const size_t i, const VertexR2 &p) const {

            const size_t n  = _v.size();
            const size_t ip = (i + 1) % n;

            const VertexR2 &v1 = _v[i];
            const VertexR2 &v2 = _v[ip];

            VertexR2 dir;
            dir = v2 - v1;
            const double length = dir.length();

            assert(fabs(length) > 0.0);
            dir *= 1.0 / length;

            const double actualDist = (v1 - p - (v1 - p).scalarProduct(dir) * dir).length();

            const double dist1 = (v1 - p).length();
            const double dist2 = (v2 - p).length();

            return actualDist < 1.0e-10 && dist1 <= length && dist2 <= length;
        }

        // Assign values to the Triangle variables.
        void assignValues(Triangulateio &in) const {

            const int n = (int) _v.size();
            const int N = n + (int) _p.size();

            in.numberofpoints = N;
            in.numberofpointattributes = _recordOrder ? 1 : 0;

            in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
            in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));

            if (_recordOrder)
                in.pointattributelist = (double *) malloc(in.numberofpointattributes * in.numberofpoints * sizeof(double));

            in.numberofsegments = n;
            in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
            in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));

            if (_recordOrder) {
                for (int i = 0; i < in.numberofpoints; ++i)
                    in.pointattributelist[i] = i;
            }

            for (int i = 0; i < n; ++i) {
                in.pointlist[2 * i] = _v[i].x();
                in.pointlist[2 * i + 1] = _v[i].y();

                in.pointmarkerlist[i] = 1;
                in.segmentmarkerlist[i] = 1;

                in.segmentlist[2 * i] = i;

                if (i == n - 1) in.segmentlist[2 * i + 1] = 0;
                else in.segmentlist[2 * i + 1] = i + 1;
            }

            for (int i = n; i < N; ++i) {

                in.pointlist[2 * i] = _p[i - n].x();
                in.pointlist[2 * i + 1] = _p[i - n].y();
                in.pointmarkerlist[i] = isOnBoundary(_p[i - n]) ? 1 : 0;
            }

            in.numberofholes = 0;
            in.numberofregions = 0;
        }

        // Assign flags and triangulate.
        void triangulate(Triangulateio &mid, Triangulateio &out) const {

            out.pointlist = (double *) NULL;
            out.pointmarkerlist = (int *) NULL;

            if (_recordOrder) out.pointattributelist = (double *) NULL;

            out.trianglelist = (int *) NULL;
            out.triangleattributelist = (double *) NULL;

            out.segmentlist = (int *) NULL;
            out.segmentmarkerlist = (int *) NULL;

            out.edgelist = (int *) NULL;
            out.edgemarkerlist = (int *) NULL;

            std::stringstream buf;
            buf.precision(100);
            buf.setf(std::ios::fixed, std::ios::floatfield);

            if (!_verbose) buf << "Q";
            else buf << "V";

            buf << "ez";

            if (!_allowEdgeRefinement) buf << "YY";
            else if (!_allowBoundaryRefinement) buf << "Y";

            if (_constrainAngle) buf << "q";

            if (_planarGraph) buf << "p";

            if (_allowSteinerPoints >= 0) buf << "S" << _allowSteinerPoints;

            if (_maxArea > 0.0) buf << "a" << _maxArea;

            char *str = new char[buf.str().size() + 1];
            strcpy(str, buf.str().c_str());

            tl::triangulate(str, &mid, &out, (Triangulateio *) NULL);

            delete[] str;
        }

        // Set/save vertices and faces obtained from Triangle.
        void setTriangulation(Triangulateio &tri, std::vector<VertexR2> &v, std::vector<Face> &f) const {

            v.clear();
            f.clear();

            VertexR2 newV;
            for (int i = 0; i < tri.numberofpoints; ++i) {

                newV = VertexR2(tri.pointlist[i * 2], tri.pointlist[i * 2 + 1]);
                v.push_back(newV);
            }

            Face newF;
            for (int i = 0; i < tri.numberoftriangles; ++i) {

                newF.v[0] = tri.trianglelist[i * 3];
                newF.v[1] = tri.trianglelist[i * 3 + 1];
                newF.v[2] = tri.trianglelist[i * 3 + 2];

                f.push_back(newF);
            }   
        }

        // Free memory.
        void freeIn(Triangulateio &in) const {

            free(in.pointlist);
            free(in.pointmarkerlist);

            free(in.segmentlist);
            free(in.segmentmarkerlist);

            if (_recordOrder) free(in.pointattributelist);    
        }

        void freeOut(Triangulateio &out) const {

            free(out.pointlist);
            free(out.pointmarkerlist);

            if (_recordOrder) free(out.pointattributelist);

            free(out.trianglelist);
            free(out.triangleattributelist);

            free(out.segmentlist);
            free(out.segmentmarkerlist);

            free(out.edgelist);
            free(out.edgemarkerlist);     
        }
    };

} // namespace gbc

#endif // GBC_TRIANGULATOR2_HPP
