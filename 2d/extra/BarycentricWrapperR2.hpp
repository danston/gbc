// Author: Dmitry Anisimov, danston@ymail.com.
// Copyright Dmitry Anisimov (c) 2016.

// README:
/*

    The classes in this file depend on two other classes that can be found in the extra folder:
    1. VertexExpressionsR2.hpp
    2. VertexR2.hpp

    This is a wrapper for my data structures that takes your point class with two functions implemented (see example in the main.cpp):
    1. x() - returns x coordinate of the point,
    2. y() - returns y coordinate of the point
    and gives you back the barycentric coordinates in the standard std::vector.

*/

#ifndef GBC_BARYCENTRICWRAPPERR2_HPP
#define GBC_BARYCENTRICWRAPPERR2_HPP

// STL includes.
#include <vector>
#include <cassert>

// Local includes.
#include "../extra/VertexR2.hpp"

namespace gbc {

    // This is a wrapper for all pointwise coordinates.
    template<class YourPointClass, class MyBarycentricClass>
    class PointwiseWrapperR2 {

    private:
        // Some typedefs.
        typedef YourPointClass Point;
        typedef MyBarycentricClass Coordinates;

    public:
        // Constructor.
        PointwiseWrapperR2(const std::vector<Point> &v, const double tol = 1.0e-10) : _tol(tol) { 

            _v.resize(v.size());

            for (size_t i = 0; i < v.size(); ++i) 
                _v[i] = VertexR2(v[i].x(), v[i].y());
        }

        // Function that computes coordinates b at a point p.
        inline void compute(const Point &p, std::vector<double> &b) {

            Coordinates coords(_v, _tol);
            coords.compute(VertexR2(p.x(), p.y()), b);
        }

        // Compute coordinates b at all points p in the vector.
        inline void compute(const std::vector<Point> &p, std::vector<std::vector<double> > &b) const {
            
            b.clear();
            b.resize(p.size());
            for (size_t i = 0; i < p.size(); ++i) compute(p[i], b[i]);
        }

    private:
        // Tolerance.
        const double _tol;

        // Vertices of the polygon.
        std::vector<VertexR2> _v;
    };

    // This is a wrapper for all mesh-based coordinates.
    template<class YourPointClass, class MyBarycentricClass>
    class MeshbasedWrapperR2 {

    private:
        // Some typedefs.
        typedef YourPointClass Point;
        typedef MyBarycentricClass Coordinates;

    public:
        // Constructor.
        MeshbasedWrapperR2(const std::vector<Point> &v, const double tol = 1.0e-10) : _tol(tol) { 

            _v.resize(v.size());

            for (size_t i = 0; i < v.size(); ++i) 
                _v[i] = VertexR2(v[i].x(), v[i].y());
        }

        // Function that computes coordinates bb at all points p.
        void compute(const std::vector<Point> &p, std::vector<std::vector<double> > &bb) {

            std::vector<VertexR2> tmp(p.size());
            for (size_t i = 0; i < p.size(); ++i) tmp[i] = VertexR2(p[i].x(), p[i].y());

            Coordinates coords(_v, _tol);
            coords.compute(tmp, bb);
        }

        // Compute coordinates bb at the vertices of the internal mesh with the given edgeLength
        // of the average triangle in this mesh.
        void compute(const double edgeLength, std::vector<std::vector<double> > &bb) {

            Coordinates coords(_v, _tol);
            coords.compute(edgeLength, bb);
        }

    private:
        // Tolerance.
        const double _tol;

        // Vertices of the polygon.
        std::vector<VertexR2> _v;
    };

} // namespace gbc

#endif // GBC_BARYCENTRICWRAPPERR2_HPP
