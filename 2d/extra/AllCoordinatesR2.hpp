// Author: Dmitry Anisimov, danston@ymail.com.
// Copyright Dmitry Anisimov (c) 2016.

// README:
/*

    This is an example of how to run my code for all coordinate classes at the same time.

*/

#ifndef GBC_ALLCOORDINATESR2_HPP
#define GBC_ALLCOORDINATESR2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <list>
#include <string>

// Local includes.
#include "IsolinerR2.hpp"

#include "../coords/WachspressR2.hpp"
#include "../coords/DiscreteHarmonicR2.hpp"
#include "../coords/MeanValueR2.hpp"
#include "../coords/ThreePointsR2.hpp"
#include "../coords/MetricR2.hpp"
#include "../coords/BilinearR2.hpp"
#include "../coords/PoissonR2.hpp"
#include "../coords/CubicMeanValueR2.hpp"
#include "../coords/GordonWixomR2.hpp"
#include "../coords/PositiveGordonWixomR2.hpp"
#include "../coords/PositiveMeanValueR2.hpp"
#include "../coords/MaximumEntropyR2.hpp"
#include "../coords/AffineR2.hpp"
#include "../coords/LaplaceR2.hpp"
#include "../coords/SibsonR2.hpp"
#include "../coords/HarmonicR2.hpp"
#include "../coords/LocalR2.hpp"

namespace gbc {

    // Compute all coordinates in R2.
    class AllCoordinatesR2 {

    public:
        // Constructor.
        AllCoordinatesR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) 
        : _v(v), _tol(tol), _power(1.0), _numSamples(1024), _edgeLength(-1.0) { }

        // Print coordinates.
        void print(const std::string &path,
                   const std::vector<VertexR2> &p, 
                   const size_t coordInd = 0) {

            std::vector<VertexR2> tmp(p);

            const double scale = 1000.0;

            std::vector<BarycentricCoordinatesR2 *> all;

            all.push_back(new WachspressR2(_v, _tol));
            all.push_back(new DiscreteHarmonicR2(_v, _tol));
            all.push_back(new MeanValueR2(_v, _tol));

            ThreePointsR2* tpc = new ThreePointsR2(_v, _tol);
            tpc->setPower(_power);
            all.push_back(tpc);

            all.push_back(new MetricR2(_v, _tol));
            all.push_back(new BilinearR2(_v, _tol));
            all.push_back(new PoissonR2(_v, _tol));
            all.push_back(new CubicMeanValueR2(_v, _tol));
            
            GordonWixomR2* gwc = new GordonWixomR2(_v, _tol);
            gwc->setIntegralPrecision(_numSamples);
            all.push_back(gwc);

            PositiveGordonWixomR2* pgwc = new PositiveGordonWixomR2(_v, _tol);
            pgwc->setIntegralPrecision(_numSamples);
            all.push_back(pgwc);

            all.push_back(new PositiveMeanValueR2(_v, _tol));
            all.push_back(new MaximumEntropyR2(_v, _tol));
            all.push_back(new AffineR2(_v, _tol));

            // all.push_back(new LaplaceR2(_v, _tol)); // does not work
            // all.push_back(new SibsonR2(_v, _tol)); // does not work
            
            // all.push_back(new HarmonicR2(_v, _tol)); // does not work

            // LocalR2* lbc = new LocalR2(_v, _tol); // does not work
            // lbc->setEdgeLength(_edgeLength);
            // all.push_back(lbc);

            std::vector<double> isoValues(9);
            for (size_t i = 0; i < 9; ++i)
                isoValues[i] = double(i + 1) * 0.1;

            std::vector<std::vector<std::list<VertexR2> > > contours;

            // IsolinerR2 iso(p, coordInd);

            const size_t k = all.size();
            for (size_t i = 0; i < k; ++i) {

                all[i]->bc(tmp);

                // iso.getContours(isoValues, contours);
                saveEps(path, all[i]->name(), coordInd, contours, scale);
            }
            std::cout << "\n";
        }

        // Set some parameters used in the coordinate classes.
        inline void setPower(const double power) {
            _power = power;
        }

        inline void setIntegralPrecision(const size_t numSamples) {

            assert(numSamples > 0);
            _numSamples = numSamples;
        }

        inline void setEdgeLength(const double newEdgeLength) {

            assert(newEdgeLength > 0.0);
            _edgeLength = newEdgeLength;
        }

    private:
        // Vertices of the polygon.
        const std::vector<VertexR2> &_v;

        // Tolerance.
        const double _tol;

        // Some parameters used in the coordinate classes.
        double _power;
        size_t _numSamples;
        double _edgeLength;

        // Save eps with contours.
        void saveEps(const std::string &path, const std::string &name, const size_t coordInd, 
                     const std::vector<std::vector<std::list<VertexR2> > > &contours, const double scale) const {

            std::cout << "\n" << name << "\n";

            assert(_v.size() != 0);

            if (path.empty()) return;

            std::string fileName = path;

            fileName.append(name + ".eps");

            std::ofstream saver(fileName.c_str(), std::ios_base::out);

            if (!saver) {
                std::cout << "\nError saving eps file with contours!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Compute bounding box.
            VertexR2 minB, maxB;
            boundingBox(minB, maxB);

            if ((maxB - minB).length() < 10.0 && scale == 1.0) scale *= 1000;

            // Header.
            setHeader(minB.x() * scale, minB.y() * scale, maxB.x() * scale, maxB.y() * scale, saver, "contours");

            // Start private namespace.
            saver << "0 dict begin gsave\n\n";

            // Save polygon with contours.

            // Save polygon.
            drawConnectedPath(saver, true, scale, Color(0, 0, 0));
            drawScatteredPoints(saver, scale, Color(0, 0, 0));

            // Save contours.
            const size_t numC = contours.size();
            for (size_t i = 0; i < numC; ++i) drawContours(contours[i], saver, scale, Color(0, 0, 0));

            // Finish private namespace.
            saver << "grestore end\n\n";
            saver << "%%EOF\n";

            // Close the file.
            saver.close();
        }
    };

} // namespace gbc

#endif // GBC_ALLCOORDINATESR2_HPP
