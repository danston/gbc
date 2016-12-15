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

    See also an example below with my wrapper for your point class.

*/

#include "./coords/MeanValueR2.hpp"
#include "./coords/HarmonicR2.hpp"
#include "./extra/BarycentricWrapperR2.hpp"
#include "./extra/TriangulatorR2.hpp"
#include "./extra/MeshR2.hpp"

// Example of a custom point class for the wrapper used below.
class YourPointClass {

public:
    // Constructor.
    YourPointClass(const double x = 0.0, const double y = 0.0) : _x(x), _y(y) { }

    // Return x coordinate.
    inline double &x() {
        return _x;
    }

    // Return const x coordinate.
    inline double x() const {
        return _x;
    }

    // Return y coordinate.
    inline double &y() {
        return _y;
    }

    // Return const y coordinate.
    inline double y() const {
        return _y;
    }

private:
    // Coordinates x and y of the point.
    double _x, _y;
};

// Example.
int main() {

    using namespace gbc;


    // Polygon.
    std::vector<YourPointClass> square(4);

    square[0] = YourPointClass(0.0, 0.0);
    square[1] = YourPointClass(1.0, 0.0);
    square[2] = YourPointClass(1.0, 1.0);
    square[3] = YourPointClass(0.0, 1.0);


    // Pointwise example.

    // Evaluation point.
    YourPointClass centre(0.5001, 0.5001);

    // Storage for the computed barycentric coordinates.
    std::vector<double> b;

    // Compute barycentric coordinates.
    PointwiseWrapperR2<YourPointClass, MeanValueR2> pbc(square);
    pbc.compute(centre, b);

    // Output the resulting coordinates.
    std::cout << "\nPointwise result: ";
    for (size_t i = 0; i < b.size(); ++i) std::cout << b[i] << " ";
    std::cout << "\n\n";


    // Mesh-based example.

    // Evaluation points.
    std::vector<VertexR2> sq(4);
    for (size_t i = 0; i < 4; ++i)
        sq[i] = VertexR2(square[i].x(), square[i].y());

    const double edgeLength = 0.05;

    TriangulatorR2 tri(sq, edgeLength, true);
    tri.setPlanarGraph(true);

    std::vector<VertexR2> tp;
    std::vector<Face> tf;

    tri.generate(tp, tf);

    MeshR2 mesh;
    mesh.initialize(tp, tf);

    // Here you can define your own set of the evaluation points.
    // This set MUST EXCLUDE the polygon's vertices!
    std::vector<YourPointClass> queries;

    const size_t tps = tp.size();
    assert(tps == mesh.numVertices());

    for (size_t i = 0; i < tps; ++i) {
        
        if (mesh.vertices()[i].type == INTERIOR || 
            mesh.vertices()[i].type == FLAT) {
        
            queries.push_back(YourPointClass(tp[i].x(), tp[i].y()));
        }
    }

    // Storage for the computed barycentric coordinates.
    std::vector<std::vector<double> > bb;

    // Compute barycentric coordinates.
    MeshbasedWrapperR2<YourPointClass, HarmonicR2> mbc(square);
    // mbc.compute(queries, bb);

    // You can also compute the coordinates given some edgeLength
    // of the average triangle in the internal triangle mesh.
    mbc.compute(edgeLength, bb);

    // Output the resulting coordinates.
    std::cout << "Mesh-based result: \n\n";
    for (size_t i = 0; i < bb.size(); ++i) {

        for (size_t j = 0; j < bb[i].size(); ++j)
            std::cout << bb[i][j] << " ";

        std::cout << std::endl;
    }
    std::cout << "\n";
}
