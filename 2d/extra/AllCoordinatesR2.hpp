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
                   std::vector<VertexR2> &tp,
                   const std::vector<Face> &tf, 
                   const size_t coordInd = 0) {

            // Set all coordinate classes.
            std::vector<BarycentricCoordinatesR2 *> all;

            all.push_back(new WachspressR2(_v, _tol));

            all.push_back(new DiscreteHarmonicR2(_v, _tol));
            all.push_back(new MeanValueR2(_v, _tol));

            ThreePointsR2* tpc = new ThreePointsR2(_v, _tol);
            tpc->setPower(_power);
            all.push_back(tpc);

            all.push_back(new MetricR2(_v, _tol));
            
            if (_v.size() == 4) all.push_back(new BilinearR2(_v, _tol));
            
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

            all.push_back(new LaplaceR2(_v, _tol));
            all.push_back(new SibsonR2(_v, _tol));
            
            HarmonicR2* hmc = new HarmonicR2(_v, _tol);
            hmc->setMesh(tp, tf);
            all.push_back(hmc);

            LocalR2* lbc = new LocalR2(_v, _tol);
            lbc->setEdgeLength(_edgeLength);
            lbc->setMesh(tp, tf);
            all.push_back(lbc);

            // Set iso values for the contour extraction from the coordinate basis function.
            std::vector<double> isoValues(9);
            for (size_t i = 0; i < 9; ++i)
                isoValues[i] = double(i + 1) * 0.1;

            std::vector<std::vector<std::list<VertexR2> > > contours;

            const double scale = 300.0;

            // Compute coordinates, create contours and plot them.
            std::cout << "Illustrating coordinates...\n";
            const size_t k = all.size();
            for (size_t i = 0; i < k; ++i) {

                all[i]->bc(tp);

                IsolinerR2 iso(tp, tf, coordInd);
                iso.getContours(isoValues, contours);

                saveEps(path, all[i]->name(), contours, scale);
            }
            std::cout << "...finished!\n\n";
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
        void saveEps(const std::string &path, const std::string &name,
                     const std::vector<std::vector<std::list<VertexR2> > > &contours, double scale) const {

            std::cout << name << "\n";

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

            if ((maxB - minB).length() < 10.0 && scale == 1.0) scale *= 1000.0;

            // Header.
            setHeader(minB.x() * scale, minB.y() * scale, maxB.x() * scale, maxB.y() * scale, saver, "contours");

            // Start private namespace.
            saver << "0 dict begin gsave\n\n";

            // Save polygon with contours.

            // Save polygon.
            drawPolygon(saver, scale);
            drawPolygonVertices(saver, 1.0, scale);

            // Save contours.
            const size_t numC = contours.size();
            for (size_t i = 0; i < numC; ++i) drawContours(saver, contours[i], scale);

            // Finish private namespace.
            saver << "grestore end\n\n";
            saver << "%%EOF\n";

            // Close the file.
            saver.close();
        }

        // Set eps header.
        void setHeader(const double llx, const double lly, const double urx, const double ury,
                       std::ofstream &ostr, const std::string &title) const {

            ostr << "%!PS-Adobe-3.0 EPSF-3.0\n";

            ostr << "%%BoundingBox: " << llx << " " << lly << " " << urx << " " << ury << "\n";

            ostr << "%%Pages: 1\n";

            ostr << "%%Creator: Dmitry Anisimov, danston@ymail.com\n";

            ostr << "%%Title: " << title.c_str() << "\n";

            ostr << "%%EndComments\n";

            ostr << "%%EndProlog\n\n";

            ostr << "%%Page: 1 1\n\n";
        }

        // Draw a polygon.
        void drawPolygon(std::ofstream &ostr, const double scale = 1.0) const {

            const size_t numV = _v.size();

            ostr << _v[0].x() * scale << " " << _v[0].y() * scale << " moveto\n";
            for (size_t i = 1; i < numV; ++i) ostr << _v[i].x() * scale << " " << _v[i].y() * scale << " lineto\n";

            ostr << _v[0].x() * scale << " " << _v[0].y() * scale << " lineto\n\n";
            ostr << "closepath\n\n";

            ostr << 0.0 << " setgray\n";
            ostr << "1 setlinewidth\n";
            ostr << "stroke\n\n";
        }

        // Draw the polygon's vertices.
        void drawPolygonVertices(std::ofstream &ostr, const double radius = 1.0, const double scale = 1.0) const {

            const size_t numV = _v.size();
            for (size_t i = 0; i < numV; ++i) drawDisc(ostr, _v[i], radius, scale);
        }

        // Draw a filled disc.
        void drawDisc(std::ofstream &ostr,
                      const VertexR2 &center,
                      const double radius,
                      const double scale = 1.0) const {

            ostr << 0.0 << " setgray\n";
            ostr << "0 setlinewidth\n\n";

            ostr << center.x() * scale << " " << center.y() * scale << " " << radius << " 0 360 arc closepath\n\n";
            ostr << "gsave\n";

            ostr << 0.0 << " setgray fill\n";

            ostr << "grestore\n";
            ostr << "stroke\n\n";
        }

        // Draw contours.
        void drawContours(std::ofstream &ostr,
                          const std::vector<std::list<VertexR2> > &p,
                          const double scale = 1.0,
                          const double lineWidth = 1.0,
                          const bool dashed = false) const {

            const size_t numP = p.size();

            for (size_t i = 0; i < numP; ++i) {

                typename std::list<VertexR2>::const_iterator it = p[i].begin();

                ostr << it->x() * scale << " " << it->y() * scale << " moveto\n";
                ++it;

                for (; it != p[i].end(); ++it) ostr << it->x() * scale << " " << it->y() * scale << " lineto\n";

                ostr << 0.0 << " setgray\n";

                if (dashed) ostr << "[4 1] 0 setdash\n";
                else ostr << "[] 0 setdash\n";

                ostr << lineWidth << " setlinewidth\n";
                ostr << "stroke\n\n";
            }
        }

        // Compute bounding box of the set of points.
        void boundingBox(VertexR2 &minB, VertexR2 &maxB) const {

            const size_t n = _v.size();

            double minX = std::numeric_limits<double>::max();
            double minY = std::numeric_limits<double>::max();

            double maxX = -std::numeric_limits<double>::max();
            double maxY = -std::numeric_limits<double>::max();

            for (size_t i = 0; i < n; ++i) {

                minX = std::min(minX, _v[i].x());
                maxX = std::max(maxX, _v[i].x());
                minY = std::min(minY, _v[i].y());
                maxY = std::max(maxY, _v[i].y());
            }

            minB = VertexR2(minX, minY);
            maxB = VertexR2(maxX, maxY);
        }
    };

} // namespace gbc

#endif // GBC_ALLCOORDINATESR2_HPP
