//
// Created by Dmitry Anisimov on 11/03/16.
//

#ifndef GBCBOOK_BARYCENTRICPROPERTIESR2_H
#define GBCBOOK_BARYCENTRICPROPERTIESR2_H

// STL includes.
#include <vector>
#include <iostream>

// Local includes.
#include "../structs/VertexR2.h"

class BarycentricPropertiesR2 {

private:
    // Fields.
    const std::vector<VertexR2> &_v;

    // Tolerance.
    const double _tol;

    // Print data.
    const bool _verbose;

public:
    // Constructor.
    BarycentricPropertiesR2(const std::vector<VertexR2> &v, const double tol = 1.0e-15, const bool verbose = true)
            : _v(v), _tol(tol), _verbose(verbose) { }

    // Properties.
    bool partitionOfUnity(const std::vector<double> &phi) const;

    bool linearReproduction(const VertexR2 &p, const std::vector<double> &phi) const;

    bool boundedness(const std::vector<double> &phi) const;

    void checkAll(const VertexR2 &p, const std::vector<double> &phi, bool &pu, bool &lr, bool &bd) const {

        pu = partitionOfUnity(phi);
        lr = linearReproduction(p, phi);
        bd = boundedness(phi);

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
            checkAll(p[i], p[i].phi(), pu, lr, bd);

            if (!pu) puRes = false;
            if (!lr) lrRes = false;
            if (!bd) bdRes = false;
        }

        if (puRes) Utils::print("FINAL: Partition of unity: SUCCESS!\n");
        else Utils::print("FINAL: Partition of unity: FAILED!\n");

        if (lrRes) Utils::print("FINAL: Linear reproduction: SUCCESS!\n");
        else Utils::print("FINAL: Linear reproduction: FAILED!\n");

        if (bdRes) Utils::print("FINAL: Boundedness: SUCCESS!\n");
        else Utils::print("FINAL: Boundedness: FAILED!\n");
    }
};

bool BarycentricPropertiesR2::partitionOfUnity(const std::vector<double> &phi) const {

    const size_t n = _v.size();

    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) sum += phi[i];

    const double diff = fabs(sum - 1.0);

    if (diff > _tol) {
        if (_verbose)
            std::cerr << "ERROR: Partition of unity difference is "
            << diff << ";" << std::endl;
        return false;
    }

    if (_verbose) std::cout << "Partition of unity: SUCCESS" << std::endl;
    return true;
}

bool BarycentricPropertiesR2::linearReproduction(const VertexR2 &x, const std::vector<double> &phi) const {

    const size_t n = _v.size();

    VertexR2 sum;
    for (size_t i = 0; i < n; ++i) sum += phi[i] * _v[i];

    const VertexR2 diff = sum - x;

    if (fabs(diff.x()) > _tol || fabs(diff.y()) > _tol) {
        if (_verbose)
            std::cerr << "ERROR: Linear reproduction difference is "
            << diff << ";" << std::endl;
        return false;
    }

    if (_verbose) std::cout << "Linear reproduction: SUCCESS" << std::endl;
    return true;
}

bool BarycentricPropertiesR2::boundedness(const std::vector<double> &phi) const {

    const size_t n = _v.size();

    for (size_t i = 0; i < n; ++i) {
        if (phi[i] < 0.0 || phi[i] > 1.0) {
            if (_verbose)
                std::cerr << "ERROR: Value out of range [0,1] with index " << i << ": "
                << phi[i] << ";" << std::endl;
            return false;
        }
    }

    if (_verbose) std::cout << "Boundedness: SUCCESS" << std::endl;
    return true;
}

#endif // GBCBOOK_BARYCENTRICPROPERTIESR2_H
