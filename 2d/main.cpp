// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

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

    std::vector<double> result;

    MyBarycentricClassR2 bary(square);
    bary.compute(centre, result);

    Below you can find other examples. For more information see README.md.

*/

// Local includes
#include "./coords/MeanValueR2.hpp"
#include "./coords/HarmonicR2.hpp"

#include "./extra/MeshR2.hpp"
#include "./extra/TriangulatorR2.hpp"
#include "./extra/AllCoordinatesR2.hpp"
#include "./extra/TestCoordinatesR2.hpp"

// Examples.
int main() {

    using namespace gbc;


    // Polygon.
    std::vector<VertexR2> poly(4);

    poly[0] = VertexR2(0.1, 0.1);
    poly[1] = VertexR2(1.0, 0.0);
    poly[2] = VertexR2(0.9, 0.9);
    poly[3] = VertexR2(0.2, 1.0);


    // 1. Pointwise example.

    // Evaluation point.
    VertexR2 centre(0.5, 0.5);

    // Storage for the computed barycentric coordinates.
    std::vector<double> b;

    // Compute mean value coordinates.
    MeanValueR2 pbc(poly);
    pbc.compute(centre, b);

    // Output the resulting coordinates.
    std::cout << "\nPointwise result: ";
    for (size_t i = 0; i < b.size(); ++i) std::cout << b[i] << " ";
    std::cout << "\n\n";


    // 2. Mesh-based example.

    // Create evaluation points.
    const double edgeLength = 0.05;

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

    std::vector<VertexR2> queries;
    std::vector<Face> faces;

    tri.generate(queries, faces);

    MeshR2 mesh;
    mesh.initialize(queries, faces);

    // Clean mesh from the polygon vertices.
    std::vector<VertexR2> cleaned;
    for (size_t i = 0; i < mesh.numVertices(); ++i) {
        if (mesh.vertices()[i].type == INTERIOR || mesh.vertices()[i].type == FLAT)
            cleaned.push_back(mesh.vertices()[i]);
    }

    // Storage for the computed barycentric coordinates.
    std::vector<std::vector<double> > bb;

    // Compute harmonic coordinates.
    // Here the set of points must exclude the polygon's vertices!
    HarmonicR2 mbc(poly);
    mbc.compute(cleaned, bb);

    // You can also compute the coordinates given some edgeLength
    // of the average triangle in the internal triangle mesh.
    // mbc.compute(edgeLength, bb);

    // Output the resulting coordinates.
    std::cout << "Mesh-based result: \n";
    for (size_t i = 0; i < bb.size(); ++i) {

        for (size_t j = 0; j < bb[i].size(); ++j)
            std::cout << bb[i][j] << " ";

        std::cout << std::endl;
    }
    std::cout << "\n";


    // 3. Example with all coordinates, where we print them in eps.
    
    // Index of the printed basis function.
    const size_t coordInd = 0;

    // Change to your path here!
    const std::string path = "/Users/path_to_the_folder/gbc/2d/out/";

    // Compute and draw all coordinates.
    AllCoordinatesR2 all(poly);

    all.setPower(0.5);             // for three-point coordinates
    all.setIntegralPrecision(100); // for Gordon-Wixom and positive Gordon-Wixom coordinates
    all.setEdgeLength(edgeLength); // for local coordinates

    all.print(path, queries, faces, coordInd);


    // 4. Example, where we test different properties of barycentric coordinates.
    std::cout << "Test coordinates:\n";

    TestCoordinatesR2 test;
    test.make();
    
    std::cout << "\n";
}
