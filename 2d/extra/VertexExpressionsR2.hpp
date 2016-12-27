// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    This is a set of stand-alone classes that are used in the VertexR2.hpp class.

*/

#ifndef GBC_VERTEXEXPRESSIONSR2_HPP
#define GBC_VERTEXEXPRESSIONSR2_HPP

// STL includes.
#include <cmath>
#include <cassert>
#include <algorithm>

namespace gbc {

    // Some forward definitions.
    template<typename E1, typename E2>
    class VertexTranslatedR2;

    template<typename E>
    class VertexRotatedR2;

    template<typename E>
    class VertexScaledNotUnifR2;

    enum VertexType { CONVEX, CONCAVE, FLAT, INTERIOR };

    // Class with vertex expressions in R2.
    template<typename E>
    class VertexExpressionR2 {

    public:
        // Constructors.
        VertexExpressionR2(const VertexType vertexType = INTERIOR) : out(-1), val(0), alpha(-1.0), type(vertexType) { }

        template<typename E1>
        VertexExpressionR2(VertexExpressionR2<E1> const &v) : out(v.out), val(v.val), alpha(v.alpha), type(v.type) { }

        // Return const x coordinate.
        virtual double x() const = 0;

        // Return const y coordinate.
        virtual double y() const = 0;

        // Return barycentric coordinate value.
        virtual double b(const int i) const = 0;

        // Return size of the barycentric coordinate vector.
        size_t size() const {
            return static_cast<E const &>(*this).size();
        }

        // Evaluation operators.
        operator E &() {
            return static_cast<E &>(*this);
        }

        operator E const &() const {
            return static_cast<const E &>(*this);
        }

        // Flags.
        int out;         // index of the outgoing edge (used in the halfedge datastructure)
        int val;         // valency of the vertex (used in the mesh class)
        double alpha;    // external angle if it is a corner vertex (used in the polygon class)
        VertexType type; // type of the vertex: 0 - CONVEX; 1 - CONCAVE; 2 - FLAT; 3 - INTERIOR (used in the mesh class)

        // Squared length.
        inline double squaredLength() const {
            return x() * x() + y() * y();
        }

        // Length - Euclidean 2-norm.
        inline double length() const {
            return sqrt(squaredLength());
        }

        // Scalar product.
        template<typename E1>
        inline double scalarProduct(const VertexExpressionR2<E1> &v) const {
            return x() * v.x() + y() * v.y();
        }

        // Cross product.
        template<typename E1>
        inline double crossProduct(const VertexExpressionR2<E1> &v) const {
            return x() * v.y() - y() * v.x();
        }

        // Signed angle.
        template<typename E1>
        inline double angle(const VertexExpressionR2<E1> &v) const {
            return atan2(this->crossProduct(v), this->scalarProduct(v));
        }

        // Translate vertex.
        template<typename E1>
        VertexTranslatedR2<E, E1> translated(const VertexExpressionR2<E1> &v) const {
            return VertexTranslatedR2<E, E1>(*this, v);
        }

        // Rotate vertex.
        VertexRotatedR2<E> rotated(const double alpha) const {
            return VertexRotatedR2<E>(alpha, *this);
        }

        // Non uniform scaling of the vertex.
        VertexScaledNotUnifR2<E> scaled(const double scaleX, const double scaleY) const {
            return VertexScaledNotUnifR2<E>(scaleX, scaleY, *this);
        }
    };

    // Multiplication of the vertex by a scalar.
    template<typename E>
    class VertexScaledR2 : public VertexExpressionR2<VertexScaledR2<E> > {

    private:
        // Superclass.
        typedef VertexExpressionR2<VertexScaledR2<E> > super;

        // Fields.
        const double _scale; // scaling factor
        E const &_v;

    public:
        // Constructor.
        VertexScaledR2(const double scale, VertexExpressionR2<E> const &v) : super(v), _scale(scale), _v(v) { }

        // Return size.
        size_t size() const {
            return _v.size();
        }

        // Return const x coordinate.
        double x() const {
            return _v.x() * _scale;
        }

        // Return const y coordinate.
        double y() const {
            return _v.y() * _scale;
        }

        // Return barycentric coordinate value.
        double b(const int i) const {

            assert(i >= 0 && i < (int) _v.size());
            return _v.b(i) * _scale;
        }
    };

    // Addition of two vertices.
    template<typename E1, typename E2>
    class VertexAdditionR2 : public VertexExpressionR2<VertexAdditionR2<E1, E2> > {

    private:
        // Superclass.
        typedef VertexExpressionR2<VertexAdditionR2<E1, E2> > super;

        // Fields.
        E1 const &_u;
        E2 const &_v;

        size_t _size;

    public:
        // Constructor.
        VertexAdditionR2(VertexExpressionR2<E1> const &u, VertexExpressionR2<E2> const &v) : super(), _u(u), _v(v) {
            _size = std::min(u.size(), v.size());
        }

        // Return size.
        size_t size() const {
            return _size;
        }

        // Return const x coordinate.
        double x() const {
            return _u.x() + _v.x();
        }

        // Return const y coordinate.
        double y() const {
            return _u.y() + _v.y();
        }

        // Return barycentric coordinate value.
        double b(const int i) const {

            assert(i >= 0 && i < (int) _u.size());
            assert(i >= 0 && i < (int) _v.size());

            return _u.b(i) + _v.b(i);
        }
    };

    // Subtraction of two vertices.
    template<typename E1, typename E2>
    class VertexDifferenceR2 : public VertexExpressionR2<VertexDifferenceR2<E1, E2> > {

    private:
        // Superclass.
        typedef VertexExpressionR2<VertexDifferenceR2<E1, E2> > super;

        // Fields.
        E1 const &_u;
        E2 const &_v;

        size_t _size;

    public:
        // Constructor.
        VertexDifferenceR2(VertexExpressionR2<E1> const &u, VertexExpressionR2<E2> const &v) : super(), _u(u), _v(v) {
            _size = std::min(u.size(), v.size());
        }

        // Return size.
        size_t size() const {
            return _size;
        }

        // Return const x coordinate.
        double x() const {
            return _u.x() - _v.x();
        }

        // Return const y coordinate.
        double y() const {
            return _u.y() - _v.y();
        }

        // Return barycentric coordinate value.
        double b(const int i) const {

            assert(i >= 0 && i < (int) _u.size());
            assert(i >= 0 && i < (int) _v.size());

            return _u.b(i) - _v.b(i);
        }
    };

    // Negation of the vertex.
    template<typename E>
    class VertexNegatedR2 : public VertexExpressionR2<VertexNegatedR2<E> > {

    private:
        // Superclass.
        typedef VertexExpressionR2<VertexNegatedR2<E> > super;

        // Fields.
        E const &_v;

    public:
        // Constructor.
        VertexNegatedR2(VertexExpressionR2<E> const &v) : super(), _v(v) { }

        // Return size.
        size_t size() const {
            return _v.size();
        }

        // Return const x coordinate.
        double x() const {
            return -_v.x();
        }

        // Return const y coordinate.
        double y() const {
            return -_v.y();
        }

        // Return barycentric coordinate value.
        double b(const int i) const {

            assert(i >= 0 && i < (int) _v.size());
            return -_v.b(i);
        }
    };

    // Translated vertex.
    template<typename E1, typename E2>
    class VertexTranslatedR2 : public VertexExpressionR2<VertexTranslatedR2<E1, E2> > {

    private:
        // Superclass.
        typedef VertexExpressionR2<VertexTranslatedR2<E1, E2> > super;

        // Fields.
        E1 const &_u;
        E2 const &_v;

    public:
        // Constructor.
        VertexTranslatedR2(VertexExpressionR2<E1> const &u, VertexExpressionR2<E2> const &v) : super(), _u(u), _v(v) { }

        // Return size.
        size_t size() const {
            return 0;
        }

        // Return const x coordinate.
        double x() const {
            return _u.x() + _v.x();
        }

        // Return const y coordinate.
        double y() const {
            return _u.y() + _v.y();
        }

        // Return barycentric coordinate value.
        double b(const int i) const {

            assert(i >= 0 && i < (int) _u.size());
            return _u.b(i);
        }
    };

    // Rotated vertex.
    template<typename E>
    class VertexRotatedR2 : public VertexExpressionR2<VertexRotatedR2<E> > {

    private:
        // Superclass.
        typedef VertexExpressionR2<VertexRotatedR2<E> > super;

        // Fields.
        const double _alpha; // rotation angle
        E const &_v;

    public:
        // Constructor.
        VertexRotatedR2(const double alpha, VertexExpressionR2<E> const &v) : super(), _alpha(alpha), _v(v) { }

        // Return size.
        size_t size() const {
            return 0;
        }

        // Return const x coordinate.
        double x() const {
            return cos(_alpha) * _v.x() - sin(_alpha) * _v.y();
        }

        // Return const y coordinate.
        double y() const {
            return sin(_alpha) * _v.x() + cos(_alpha) * _v.y();
        }

        // Return barycentric coordinate value.
        double b(const int i) const {
            assert(i >= 0 && i < (int) _v.size());
            return _v.b(i);
        }
    };

    // Not uniformly scaled vertex.
    template<typename E>
    class VertexScaledNotUnifR2 : public VertexExpressionR2<VertexScaledNotUnifR2<E> > {

    private:
        // Superclass.
        typedef VertexExpressionR2<VertexScaledNotUnifR2<E> > super;

        // Fields.
        const double _scaleX, _scaleY;
        E const &_v;

    public:
        // Constructor.
        VertexScaledNotUnifR2(const double scaleX, const double scaleY, VertexExpressionR2<E> const &v)
                : super(), _scaleX(scaleX), _scaleY(scaleY), _v(v) { }

        // Return size.
        size_t size() const {
            return 0;
        }

        // Return const x coordinate.
        double x() const {
            return _v.x() * _scaleX;
        }

        // Return const y coordinate.
        double y() const {
            return _v.y() * _scaleY;
        }

        // Return barycentric coordinate value.
        double b(const int i) const {
            assert(i >= 0 && i < (int) _v.size());
            return _v.b(i);
        }
    };

    // Overload some basic operators.
    template<typename E>
    VertexScaledR2<E> const
    operator*(const double scale, VertexExpressionR2<E> const &v) {
        return VertexScaledR2<E>(scale, v);
    }

    template<typename E>
    VertexScaledR2<E> const
    operator*(VertexExpressionR2<E> const &v, const double scale) {
        return VertexScaledR2<E>(scale, v);
    }

    template<typename E1, typename E2>
    VertexAdditionR2<E1, E2> const
    operator+(VertexExpressionR2<E1> const &u, VertexExpressionR2<E2> const &v) {
        return VertexAdditionR2<E1, E2>(u, v);
    }

    template<typename E1, typename E2>
    VertexDifferenceR2<E1, E2> const
    operator-(VertexExpressionR2<E1> const &u, VertexExpressionR2<E2> const &v) {
        return VertexDifferenceR2<E1, E2>(u, v);
    }

    template<typename E>
    VertexNegatedR2<E> const
    operator-(VertexExpressionR2<E> const &v) {
        return VertexNegatedR2<E>(v);
    }

} // namespace gbc

#endif // GBC_VERTEXEXPRESSIONSR2_HPP
