// Author: Dmitry Anisimov, danston@ymail.com.
// Copyright Dmitry Anisimov (c) 2016.

// README:
/*

    This is a test class for testing generalized barycentric coordinates.

*/

#ifndef GBC_TESTCOORDINATESR2_HPP
#define GBC_TESTCOORDINATESR2_HPP

// STL includes.
#include <vector>
#include <iostream>
#include <cmath>

// Local includes.
#include "Face.hpp"
#include "VertexR2.hpp"
#include "TriangulatorR2.hpp"
#include "BarycentricCoordinatesR2.hpp"

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

    // Different theoretical properties of generalized barycentric coordinates.
    class BarycentricPropertiesR2 {

    public:
        // Constructor.
        BarycentricPropertiesR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10, const bool verbose = true)
                : _v(v), _tol(tol), _verbose(verbose) { }

        // Properties.
        bool partitionOfUnity(const std::vector<double> &b) const {

            const size_t n = _v.size();

            double sum = 0.0;
            for (size_t i = 0; i < n; ++i) sum += b[i];

            const double diff = fabs(sum - 1.0);

            if (diff > _tol) {
                if (_verbose)
                    std::cerr << "ERROR: Partition of unity difference is " << diff << ";" << std::endl;
                return false;
            }

            if (_verbose) std::cout << "Partition of unity: SUCCESS" << std::endl;
            return true;
        }

        bool linearReproduction(const VertexR2 &p, const std::vector<double> &b) const {

            const size_t n = _v.size();

            VertexR2 sum;
            for (size_t i = 0; i < n; ++i) sum += b[i] * _v[i];

            const VertexR2 diff = sum - p;

            if (fabs(diff.x()) > _tol || fabs(diff.y()) > _tol) {
                if (_verbose)
                    std::cerr << "ERROR: Linear reproduction difference is " << diff << ";" << std::endl;
                return false;
            }

            if (_verbose) std::cout << "Linear reproduction: SUCCESS" << std::endl;
            return true;
        }

        bool boundedness(const std::vector<double> &b) const {

            const size_t n = _v.size();

            for (size_t i = 0; i < n; ++i) {
                if (b[i] < 0.0 || b[i] > 1.0) {
                    if (_verbose)
                        std::cerr << "ERROR: Value out of range [0,1] with index " << i << ": " << b[i] << ";" << std::endl;
                    return false;
                }
            }

            if (_verbose) std::cout << "Boundedness: SUCCESS" << std::endl;
            return true;
        }

        void checkAll(const VertexR2 &p, const std::vector<double> &b, bool &pu, bool &lr, bool &bd) const {

            pu = partitionOfUnity(b);
            lr = linearReproduction(p, b);
            bd = boundedness(b);

            if (_verbose) {

                std::cout << "\nPartition of unity: " << (pu ? "SUCCESS\n" : "FAILED\n");
                std::cout << "Linear reproduction: " << (lr ? "SUCCESS\n" : "FAILED\n");

                std::cout << "Boundedness: " << (bd ? "SUCCESS\n\n" : "FAILED\n");
            }
        }

        void checkAll(const std::vector<VertexR2> &p) {

            bool puRes = true, lrRes = true, bdRes = true;
            for (size_t i = 0; i < p.size(); ++i) {

                bool pu, lr, bd;
                checkAll(p[i], p[i].b(), pu, lr, bd);

                if (!pu) puRes = false;
                if (!lr) lrRes = false;
                if (!bd) bdRes = false;
            }

            if (puRes) std::cout << "FINAL: Partition of unity: SUCCESS!\n";
            else std::cout << "FINAL: Partition of unity: FAILED!\n";

            if (lrRes) std::cout << "FINAL: Linear reproduction: SUCCESS!\n";
            else std::cout << "FINAL: Linear reproduction: FAILED!\n";

            if (bdRes) std::cout << "FINAL: Boundedness: SUCCESS!\n";
            else std::cout << "FINAL: Boundedness: FAILED!\n";
        }

    private:
        // Vertices of the polygon.
        const std::vector<VertexR2> &_v;

        // Tolerance.
        const double _tol;

        // Print data.
        const bool _verbose;
    };

    // Class that verifies different properties and behaviour of generalzied barycentric coordinates in R2.
    class TestCoordinatesR2 {

    public:
        // Constructor.
        TestCoordinatesR2() : _tol(1.0e-10) { }

        // Function with different tests.
        void make() {

            // Edge length.
            const double edgeLength = 0.05;

            // Check barycentric properties on a irregular quad.
            std::vector<VertexR2> quad(4);

            quad[0] = VertexR2(0.1, 0.1);
            quad[1] = VertexR2(1.0, 0.0);
            quad[2] = VertexR2(0.9, 0.9);
            quad[3] = VertexR2(0.2, 1.0);

            std::vector<VertexR2> tp;
            std::vector<Face> tf;

            createMesh(quad, tp, tf, edgeLength);
            checkProperties(quad, tp, tf, edgeLength);

            // Compare Wachspress (WP), discrete harmonic (DH), Sibson (SB), and
            // Laplace (LP) coordinates on a circular polygon.
            std::vector<VertexR2> circular(7);

            circular[0] = VertexR2(0.745481549466664, 0.935590184544425);
            circular[1] = VertexR2(0.15182907469161, 0.858855133403299);
            circular[2] = VertexR2(0.000017768559782, 0.504215239502086);
            circular[3] = VertexR2(0.170481476091579, 0.123944761502761);
            circular[4] = VertexR2(0.767263215929374, 0.077424121120137);
            circular[5] = VertexR2(0.966378029661621, 0.319745919743979);
            circular[6] = VertexR2(0.877780665791445, 0.827538957307635);

            createMesh(circular, tp, tf, edgeLength);
            compare_WP_DH_SB_LP(circular, tp);

            // Compare mean value (MV) and positive mean value (PM) 
            // coordinates on a convex polygon.
            std::vector<VertexR2> convex(6);

            convex[0] = VertexR2(0.087272832224228, 0.554398725870316);
            convex[1] = VertexR2(0.199047453819107, 0.137461645317987);
            convex[2] = VertexR2(0.681629947054142, 0.066493631606953);
            convex[3] = VertexR2(0.857275780988953, 0.210203859371798);
            convex[4] = VertexR2(0.9, 0.4);
            convex[5] = VertexR2(0.784533566935143, 0.776173768717299);

            createMesh(convex, tp, tf, edgeLength);
            compare_MV_PM(convex, tp);
        }

    private:
        // Tolerance.
        const double _tol;

        // Create triangle mesh for a given polygon.
        void createMesh(const std::vector<VertexR2> &poly,
                        std::vector<VertexR2> &tp, 
                        std::vector<Face> &tf, 
                        const double edgeLength) {

            // Refine the polygon to create regular mesh.
            std::vector<VertexR2> refined;
            const size_t n = poly.size();

            for (size_t i = 0; i < n; ++i) {
                refined.push_back(poly[i]);

                const size_t ip = (i + 1) % n;
                const size_t numS = ceil((poly[ip] - poly[i]).length() / edgeLength);

                for (size_t j = 1; j < numS; ++j) {

                    VertexR2 vert = poly[i] + (double(j) / double(numS)) * (poly[ip] - poly[i]);
                    refined.push_back(vert);
                }
            }

            // Create mesh.
            TriangulatorR2 tri(refined, edgeLength, true);
            tri.setPlanarGraph(true);

            tp.clear();
            tf.clear();

            tri.generate(tp, tf);
        }

        // Check different barycentric properties.
        void checkProperties(const std::vector<VertexR2> &poly, 
                             std::vector<VertexR2> &tp, 
                             const std::vector<Face> &tf, 
                             const double edgeLength) const {

            std::vector<BarycentricCoordinatesR2 *> all;

            all.push_back(new WachspressR2(poly, _tol));
            all.push_back(new DiscreteHarmonicR2(poly, _tol));
            all.push_back(new MeanValueR2(poly, _tol));
            all.push_back(new ThreePointsR2(poly, _tol));
            all.push_back(new MetricR2(poly, _tol));
            all.push_back(new BilinearR2(poly, _tol));
            all.push_back(new PoissonR2(poly, _tol));
            all.push_back(new CubicMeanValueR2(poly, _tol));
            all.push_back(new GordonWixomR2(poly, _tol));
            all.push_back(new PositiveGordonWixomR2(poly, _tol));
            all.push_back(new PositiveMeanValueR2(poly, _tol));
            all.push_back(new MaximumEntropyR2(poly, _tol));
            all.push_back(new AffineR2(poly, _tol));
            all.push_back(new LaplaceR2(poly, _tol));
            all.push_back(new SibsonR2(poly, _tol));
            
            HarmonicR2* hmc = new HarmonicR2(poly, _tol);
            hmc->setMesh(tp, tf);
            all.push_back(hmc);

            LocalR2* lbc = new LocalR2(poly, _tol);
            lbc->setEdgeLength(edgeLength);
            lbc->setMesh(tp, tf);
            all.push_back(lbc);

            BarycentricPropertiesR2 bp(poly, _tol, false);
            for (size_t i = 0; i < all.size(); ++i) {

                all[i]->bc(tp);

                bool puRes = true, lrRes = true, bdRes = true;
                for (size_t j = 0; j < tp.size(); ++j) {

                    bool pu, lr, bd;
                    bp.checkAll(tp[j], tp[j].b(), pu, lr, bd);

                    if (!pu) puRes = false;
                    if (!lr) lrRes = false;
                    if (!bd) bdRes = false;
                }

                if (puRes && lrRes && bdRes) 
                    std::cout << "Properties of " << all[i]-> name() << ": SUCCESS!\n";
                else 
                    std::cout << "Properties " << puRes << " "  << lrRes << " " << bdRes
                                               <<  " of " << all[i]-> name() << ": FAILED!\n";
            }
        }

        // Compare Wachspress, discrete harmonic, Sibson, and Laplace coordinates on a regular polygon.
        // They all should give the same result.
        void compare_WP_DH_SB_LP(const std::vector<VertexR2> &circular, 
                                 std::vector<VertexR2> &tp) const {

            std::vector<std::vector<double> > wpBB;
            std::vector<std::vector<double> > dhBB;
            std::vector<std::vector<double> > sbBB;
            std::vector<std::vector<double> > lpBB;

            WachspressR2 wp(circular, _tol);
            wp.compute(tp, wpBB);

            DiscreteHarmonicR2 dh(circular, _tol);
            dh.compute(tp, dhBB);

            WachspressR2 sb(circular, _tol);
            sb.compute(tp, sbBB);

            DiscreteHarmonicR2 lp(circular, _tol);
            lp.compute(tp, lpBB);

            assert(wpBB.size() == tp.size());
            assert(wpBB.size() == dhBB.size() && 
                   dhBB.size() == sbBB.size() &&
                   sbBB.size() == lpBB.size());

            bool res = true; double value = -100.0;
            for (size_t i = 0; i < tp.size(); ++i) {
                for (size_t j = 0; j < circular.size(); ++j) {

                    if (fabs(wpBB[i][j] - dhBB[i][j]) > _tol ||
                        fabs(dhBB[i][j] - sbBB[i][j]) > _tol ||
                        fabs(sbBB[i][j] - lpBB[i][j]) > _tol) {
                    
                        value = std::max(fabs(wpBB[i][j] - dhBB[i][j]), value);
                        value = std::max(fabs(dhBB[i][j] - sbBB[i][j]), value);
                        value = std::max(fabs(sbBB[i][j] - lpBB[i][j]), value);

                        res = false;
                        break;
                    }
                }
                if (!res) break;
            }

            if (res) std::cout << "Comparison of WP, DH, SB, and LP: SUCCESS!\n";
            else std::cout << "Comparison of WP, DH, SB, and LP: " << value << " FAILED!\n";
        }

        // Compare mean value and positive mean value coordinates on a convex polygon.
        // They should give the same result.
        void compare_MV_PM(const std::vector<VertexR2> &convex, 
                            std::vector<VertexR2> &tp) const {

            std::vector<std::vector<double> > mvBB;
            std::vector<std::vector<double> > pmBB;

            MeanValueR2 mv(convex, _tol);
            mv.compute(tp, mvBB);

            PositiveMeanValueR2 pm(convex, _tol);
            pm.compute(tp, pmBB);

            assert(mvBB.size() == tp.size());
            assert(mvBB.size() == pmBB.size());

            bool res = true; double value = 0.0;
            for (size_t i = 0; i < tp.size(); ++i) {
                for (size_t j = 0; j < convex.size(); ++j) {

                    if (fabs(mvBB[i][j] - pmBB[i][j]) > _tol) {
                    
                        value = fabs(mvBB[i][j] - pmBB[i][j]);
                        res = false;

                        break;
                    }
                }
                if (!res) break;
            }

            if (res) std::cout << "Comparison of MV and PM: SUCCESS!\n";
            else std::cout << "Comparison of MV and PM: " << value << " FAILED!\n";
        }
    };

} // namespace gbc

#endif // GBC_TESTCOORDINATESR2_HPP
