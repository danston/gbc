// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    This is the vertex data structure.

    This class depends:
    1. VertexExpressionsR2.hpp

*/

#ifndef GBC_VERTEXR2_HPP
#define GBC_VERTEXR2_HPP

// STL includes.
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// Local includes.
#include "VertexExpressionsR2.hpp"

namespace gbc {

    // Vertex class in R2.
    class VertexR2 : public VertexExpressionR2<VertexR2> {

    private:
        // Superclass.
        typedef VertexExpressionR2<VertexR2> super;

        // Fields.
        double _x; // x coordinate of the vertex
        double _y; // y coordinate of the vertex

        std::vector<double> _b; // barycentric coordinates attached to each vertex

        double _tol; // tolerance

    public:
        // Typedefs.
        typedef double EntryType;

        // Constructors.
        VertexR2() : super(), _x(0.0), _y(0.0), _b(), _tol(1.0e-10) { }

        VertexR2(const std::size_t bSize) : super(), _x(0.0), _y(0.0), _b(bSize, 0.0), _tol(1.0e-10) { }

        VertexR2(const VertexR2 &v) : super(v), _x(v.x()), _y(v.y()), _b(v.b()), _tol(1.0e-10) { }

        VertexR2(const double x, const double y, const VertexType vertexType = INTERIOR) : super(vertexType), _x(x), _y(y), _b(), _tol(1.0e-10) { }

        template<typename E>
        VertexR2(VertexExpressionR2<E> const &v) {

            _b.resize(v.size());
            for (size_t i = 0; i != v.size(); ++i) _b[i] = v.b((int) i);

            _x = v.x();
            _y = v.y();

            out   = v.out;
            val   = v.val;
            type  = v.type;
            alpha = v.alpha;

            assert(val >= 0);
            assert(type >= 0 && type <= 3);

            _tol = 1.0e-10;
        }

        // Return size of the barycentric coordinate vector.
        size_t size() const {
            return _b.size();
        }

        // Return x coordinate.
        inline double &x() {
            return _x;
        }

        // Return const x coordinate.
        inline double x() const {
            return _x;
        }

        // Return y coordinate.
        inline double &y() {
            return _y;
        }

        // Return const y coordinate.
        inline double y() const {
            return _y;
        }

        // Return barycentric coordinates attached to the vertex.
        inline std::vector<double> &b() {
            return _b;
        }

        // Return const barycentric coordinates attached to the vertex.
        inline const std::vector<double> &b() const {
            return _b;
        }

        // Return barycentric coordinate value.
        inline double b(const int i) const {
            
            assert(i >= 0 && i < (int) _b.size());
            return _b[i];
        }

        // Overload some basic operators.

        // Square brakets operator by value.
        inline double operator[](const int i) const {
            
            assert(i >= 0 && i < 2);
            return (&_x)[i];
        }

        // Square brakets operator by reference.
        inline double &operator[](const int i) {
            
            assert(i >= 0 && i < 2);
            return (&_x)[i];
        }

        // Addition of two vertices without creating a new vertex.
        void operator+=(const VertexR2 &v) {
            const std::size_t bSize = _b.size();

            _x += v.x();
            _y += v.y();

            for (size_t i = 0; i < bSize; ++i) _b[i] += v.b()[i];
        }

        // Subtraction of two vertices without creating a new vertex.
        void operator-=(const VertexR2 &v) {
            const std::size_t bSize = _b.size();

            _x -= v.x();
            _y -= v.y();

            for (size_t i = 0; i < bSize; ++i) _b[i] -= v.b()[i];
        }

        // Multiplication by a constant from the right without creating a new vertex.
        void operator*=(const double scalar) {
            const std::size_t bSize = _b.size();

            _x *= scalar;
            _y *= scalar;

            for (size_t i = 0; i < bSize; ++i) _b[i] *= scalar;
        }

        // Equal equal operator.
        inline bool operator==(const VertexR2 &v) const {
            return fabs(_x - v.x()) < _tol && fabs(_y - v.y()) < _tol;
        }

        // Not equal operator.
        inline bool operator!=(const VertexR2 &v) const {
            return !(this->operator==(v));
        }

        // Stream operator.
        friend inline std::ostream &operator<<(std::ostream &ostr, const VertexR2 &v) {
            return ostr << v.x() << " " << v.y();
        }
    };

} // namespace gbc

#endif // GBC_VERTEXR2_HPP
