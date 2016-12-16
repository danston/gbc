// Author: Dmitry Anisimov, danston@ymail.com.
// Copyright Dmitry Anisimov (c) 2016.

// Some examples on how to use the coordinate classes.

// README:
/*

    A quick way to incorporate my code in your code is to implement a wrapper that converts data
    stored in your data structures to my data structures. As a quick example:

    std::vector<VertexR2> square(4);

    square[0] = VertexR2(YourPointClass.x(), YourPointClass.y());
    square[1] = VertexR2(YourPointClass.x(), YourPointClass.y());
    square[2] = VertexR2(YourPointClass.x(), YourPointClass.y());
    square[3] = VertexR2(YourPointClass.x(), YourPointClass.y());

    VertexR2 centre(YourPointClass.x(), YourPointClass.y());

    std::vector<double> res;

    MyBarycentricClassR2 bary(square);
    bary.compute(centre, res);

*/

#include "./coords/MeanValueR2.hpp"
#include "./coords/HarmonicR2.hpp"

#include "./extra/TriangulatorR2.hpp"
#include "./extra/AllCoordinatesR2.hpp"

// Example.
int main() {

    using namespace gbc;


    // Polygon.
    std::vector<VertexR2> poly(4);

    poly[0] = VertexR2(0.0, 0.0);
    poly[1] = VertexR2(1.0, 0.0);
    poly[2] = VertexR2(1.0, 1.0);
    poly[3] = VertexR2(0.0, 1.0);


    // std::vector<VertexR2> poly(6);

    // poly[0] = VertexR2(0.087272832224228, 0.554398725870316);
    // poly[1] = VertexR2(0.199047453819107, 0.137461645317987);
    // poly[2] = VertexR2(0.681629947054142, 0.066493631606953);
    // poly[3] = VertexR2(0.857275780988953, 0.210203859371798);
    // poly[4] = VertexR2(0.9, 0.4);
    // poly[5] = VertexR2(0.784533566935143, 0.776173768717299);


    // Pointwise example.

    // Evaluation point.
    VertexR2 centre(0.5, 0.5);

    // Storage for the computed barycentric coordinates.
    std::vector<double> b;

    // Compute barycentric coordinates.
    MeanValueR2 pbc(poly);
    pbc.compute(centre, b);

    // Output the resulting coordinates.
    std::cout << "\nPointwise result: ";
    for (size_t i = 0; i < b.size(); ++i) std::cout << b[i] << " ";
    std::cout << "\n\n";


    // Mesh-based example.

    // Evaluation points.
    const double edgeLength = 0.05;

    TriangulatorR2 tri(poly, edgeLength, true);
    tri.setPlanarGraph(true);

    std::vector<VertexR2> queries;
    std::vector<Face> faces;

    tri.generate(queries, faces);

    // Storage for the computed barycentric coordinates.
    std::vector<std::vector<double> > bb;

    // Compute barycentric coordinates.
    HarmonicR2 mbc(poly);
    mbc.compute(queries, bb);

    // You can also compute the coordinates given some edgeLength
    // of the average triangle in the internal triangle mesh.
    // mbc.compute(edgeLength, bb);

    // Output the resulting coordinates.
    std::cout << "Mesh-based result: \n\n";
    for (size_t i = 0; i < bb.size(); ++i) {

        for (size_t j = 0; j < bb[i].size(); ++j)
            std::cout << bb[i][j] << " ";

        std::cout << std::endl;
    }
    std::cout << "\n";


    // Example with all coordinates.
    // Eps print all coordinate basis functions with the index coordInd in the path.
    const size_t coordInd = 0;
    const std::string path = "/Users/danston/Documents/github/gbc/2d/out/";

    AllCoordinatesR2 all(poly);

    all.setPower(0.5); // for three-point coordinates
    all.setIntegralPrecision(100); // for Gordon-Wixom and positive Gordon-Wixom coordinates
    all.setEdgeLength(edgeLength); // for local coordinates

    all.print(path, queries, coordInd);
}
