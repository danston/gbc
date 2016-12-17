# Open terminal
# cd your_path/2d/bin/
# cmake -DCMAKE_BUILD_TYPE=Debug ..
# cmake -DCMAKE_BUILD_TYPE=Release ..
# make
# ./gbc

cmake_minimum_required(VERSION 3.3)

project(gbc)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -std=c++11 ")

set(SRC 
        extra/BarycentricCoordinatesR2.hpp
        extra/SegmentCoordinatesR2.hpp
        extra/VertexExpressionsR2.hpp
        extra/VertexR2.hpp
        extra/Halfedge.hpp
        extra/Face.hpp
        extra/MeshR2.hpp
        extra/TriangulatorR2.hpp
        extra/IsolinerR2.hpp
        extra/TestCoordinatesR2.hpp
        extra/libs/gaussinteg/gauss_legendre.hpp
        extra/libs/gaussinteg/gauss_legendre.cpp
        extra/libs/triangle/triangle.hpp
        extra/libs/triangle/triangle.cpp
        extra/libs/lbcsolver/DataSetup.hpp
        extra/libs/lbcsolver/DataTypes.hpp
        extra/libs/lbcsolver/LbcSolver.hpp
        coords/WachspressR2.hpp
        coords/DiscreteHarmonicR2.hpp
        coords/MeanValueR2.hpp
        coords/ThreePointsR2.hpp
        coords/MetricR2.hpp
        coords/BilinearR2.hpp
        coords/PoissonR2.hpp
        coords/CubicMeanValueR2.hpp
        coords/GordonWixomR2.hpp
        coords/PositiveGordonWixomR2.hpp
        coords/PositiveMeanValueR2.hpp
        coords/MaximumEntropyR2.hpp
        coords/AffineR2.hpp
        coords/LaplaceR2.hpp
        coords/SibsonR2.hpp
        coords/HarmonicR2.hpp
        coords/LocalR2.hpp
        main.cpp
)

add_executable(gbc ${SRC})