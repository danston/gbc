// Authors: Juyong Zhang, Bailin Deng, Zishun Liu, Giuseppe Patane, Sofien Bouaziz, Kai Hormann, and Ligang Liu, https://github.com/bldeng/LBC

#ifndef GBC_DATATYPES_HPP
#define GBC_DATATYPES_HPP

// Eigen includes.
#include "../eigen/Eigen/Dense"
#include "../eigen/Eigen/Sparse"

namespace LBC {

    typedef Eigen::MatrixXd DenseMatrix;
    typedef Eigen::Matrix3Xd DenseMatrix3X;
    typedef Eigen::MatrixXi DenseIndexMatrix;

    typedef Eigen::VectorXd DenseVector;
    typedef Eigen::VectorXi IndexVector;

    typedef Eigen::Vector2d Vector2d;
    typedef Eigen::Vector3d Vector3d;

    typedef Eigen::Triplet<double> TripletD;
    typedef Eigen::SparseMatrix<double> SparseMatrix;

} // namespace LBC

#endif // GBC_DATATYPES_HPP
