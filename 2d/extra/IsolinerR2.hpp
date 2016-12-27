// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// README:
/*

    Here we extract contours from the given basis function using the marching triangle approach.

    This class depends on:
    1. VertexExpressionsR2.hpp
    2. VertexR2.hpp
    3. Face.hpp

*/

#ifndef GBC_ISOLINER2_HPP
#define GBC_ISOLINER2_HPP

// STL includes.
#include <list>
#include <vector>
#include <cassert>

// Local includes.
#include "Face.hpp"
#include "VertexR2.hpp"

namespace gbc {

    // Look up support class.
    class LookUpTable {

    public:
        // Get middle edge.
        inline bool getMiddleEdge(const long ind,
                                  const VertexR2 &v0, const VertexR2 &v1, const VertexR2 &v2,
                                  VertexR2 &a, VertexR2 &b) {

            if (ind == 0 || ind == 7) return false;

            if (ind == 1 || ind == 6) {
                a = 0.5 * (v1 + v2);
                b = 0.5 * (v2 + v0);
                return true;
            }

            if (ind == 2 || ind == 5) {
                a = 0.5 * (v0 + v1);
                b = 0.5 * (v1 + v2);
                return true;
            }

            if (ind == 3 || ind == 4) {
                a = 0.5 * (v0 + v1);
                b = 0.5 * (v2 + v0);
                return true;
            }

            return false;
        }

        // Get edges.
        inline bool getEdges(const long ind, std::pair<size_t, size_t> &e1, std::pair<size_t, size_t> &e2) {

            if (ind == 0 || ind == 7) return false;

            if (ind == 1 || ind == 6) {
                e1.first = 1;
                e1.second = 2;
                e2.first = 2;
                e2.second = 0;
                return true;
            }

            if (ind == 2 || ind == 5) {
                e1.first = 0;
                e1.second = 1;
                e2.first = 1;
                e2.second = 2;
                return true;
            }

            if (ind == 3 || ind == 4) {
                e1.first = 0;
                e1.second = 1;
                e2.first = 2;
                e2.second = 0;
                return true;
            }

            return false;
        }
    };

    // This class extracts isolines from the given coordinate basis function given the corresponding mesh.
    class IsolinerR2 {

    public:
        // Constructor.
        IsolinerR2(const std::vector<VertexR2> &v, const std::vector<Face> &f, const size_t coordInd = 0)
                : _v(v), _f(f), _coordInd(coordInd), _tol(1.0e-10) { }

        // Function that extracts isolines.
        void getContours(const std::vector<double> &isoValues, std::vector<std::vector<std::list<VertexR2> > > &contours) {

            const size_t numC = isoValues.size();
            contours.resize(numC);

            for (size_t i = 0; i < numC; ++i)
                extractContour(isoValues[i], contours[i]);
        }

    private:
        // Vertices of the given mesh.
        const std::vector<VertexR2> &_v;

        // Faces of the given mesh.
        const std::vector<Face> &_f;

        // Index of the coordinate function.
        const size_t _coordInd;

        // Tolerance.
        const double _tol;

        // Extract one contour corresponding to the given isoValue.
        // I use the marching triangles approach here: https://en.wikipedia.org/wiki/Marching_squares
        void extractContour(const double isoValue, std::vector<std::list<VertexR2> > &result) {

            std::vector<VertexR2> p;

            const size_t numF = _f.size();

            std::vector<bool> binInd;

            createBinaryIndex(isoValue, binInd);

            LookUpTable lut;

            std::vector<bool> tmp(3);

            std::pair<size_t, size_t> e1;
            std::pair<size_t, size_t> e2;

            for (size_t i = 0; i < numF; ++i) {

                tmp[0] = binInd[_f[i].v[0]];
                tmp[1] = binInd[_f[i].v[1]];
                tmp[2] = binInd[_f[i].v[2]];

                const long uniqueIndex = binary2Decima(vector2Long(tmp));

                const bool found = lut.getEdges(uniqueIndex, e1, e2);

                if (found) getInterpolatedEdge(isoValue, e1, e2, i, p);
            }

            connectIsopieces(p, result);
        }

        // Convert a binary number to decimal.
        inline long binary2Decima(long bin) const {

            long dec = 0, base = 1, rem;
            while (bin > 0) {

                rem = bin % 10;
                dec += rem * base;
                base *= 2;
                bin /= 10;
            }

            return dec;
        }

        // Convert values from the vector with bools to a long number.
        inline long vector2Long(const std::vector<bool> &vec) {

            long result;
            std::stringstream ss;

            const size_t size = vec.size();
            for (size_t i = 0; i < size; ++i) ss << vec[i];

            ss >> result;
            return result;
        }

        // Create binary index.
        void createBinaryIndex(const double isoValue, std::vector<bool> &binInd) const {
            
            const size_t numV = _v.size();
            binInd.resize(numV);

            for (size_t i = 0; i < numV; ++i) {

                if (_v[i].b()[_coordInd] < isoValue) binInd[i] = 0;
                else binInd[i] = 1;
            }
        }

        // Get interpolated edge.
        inline void getInterpolatedEdge(const double isoValue,
                                        const std::pair<size_t, size_t> &e1,
                                        const std::pair<size_t, size_t> &e2,
                                        const size_t faceInd,
                                        std::vector<VertexR2> &p) const {

            getInterpolatedVertex(isoValue, e1, faceInd, p);
            getInterpolatedVertex(isoValue, e2, faceInd, p);
        }

        // Get interpolated vertex.
        void getInterpolatedVertex(const double isoValue,
                                   const std::pair<size_t, size_t> &e,
                                   const size_t faceInd,
                                   std::vector<VertexR2> &p) const {

            const VertexR2 &v1 = _v[_f[faceInd].v[e.first]];
            const VertexR2 &v2 = _v[_f[faceInd].v[e.second]];

            const double f1 = v1.b()[_coordInd];
            const double f2 = v2.b()[_coordInd];

            const double denom = f2 - f1;

            assert(denom != 0.0);

            double t = (isoValue - f1) / denom;

            p.push_back((1.0 - t) * v1 + t * v2);
        }

        // Connect isopieces. This function can be avoided, but then all the computed pieces will be unordered
        // and it wont be possible to simplify the final contours in Adobe Illustrator e.g.
        void connectIsopieces(const std::vector<VertexR2> &p, std::vector<std::list<VertexR2> > &result) const {

            result.clear();

            std::list<std::pair<VertexR2, VertexR2> > pool;

            for (size_t i = 0; i < p.size(); i += 2) pool.push_back(std::make_pair(p[i], p[i + 1]));

            typedef std::list<std::pair<VertexR2, VertexR2> >::iterator iter;

            while (!pool.empty()) {

                const std::pair<VertexR2, VertexR2> edge = pool.front();
                pool.pop_front();

                bool first = true;
                VertexR2 curr = edge.first;

                std::list<VertexR2> curve;
                curve.push_back(edge.second);
                curve.push_back(edge.first);

                while (true) {
                    bool found = false;
                    for (iter it = pool.begin(); it != pool.end(); ++it) {

                        if ((it->first - curr).length() < _tol) {
                            curr = it->second;
                            pool.erase(it);
                            found = true;

                            if (first) curve.push_back(curr);
                            else curve.push_front(curr);

                            break;

                        } else if ((it->second - curr).length() < _tol) {
                            curr = it->first;
                            pool.erase(it);
                            found = true;

                            if (first) curve.push_back(curr);
                            else curve.push_front(curr);

                            break;
                        }
                    }

                    if (!found && first) {
                        curr = edge.second;
                        first = false;
                    }
                    else if (!found) break;
                }

                result.push_back(curve);
            }
        }
    };

} // namespace gbc

#endif // GBC_ISOLINER2_HPP
