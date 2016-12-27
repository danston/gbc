# Generalized barycentric coordinates (GBC)

### Description

Implementation of some [generalized barycentric coordinates](http://www.inf.usi.ch/hormann/barycentric/index.html) in R2.

This package includes:

__1. Pointwise coordinates:__

A. Coordinates for polygons:
  * Bilinear coordinates
  * Wachspress coordinates
  * Discrete harmonic coordinates
  * Mean value coordinates
  * Positive mean value coordinates
  * Cubic mean value coordinates
  * Three-point coordinates
  * Metric coordinates
  * Poisson coordinates
  * Gordon-Wixom coordinates
  * Positive Gordon-Wixom coordinates
  * Maximum entropy coordinates
    
B. Coordinates for scattered points:
  * Affine coordinates
  * Sibson coordinates
  * Laplace coordinates

__2. Mesh-based coordinates:__

A. Coordinates for polygons:
  * Harmonic coordinates
  * Local barycentric coordinates

I tried to keep the code as stand alone as possible, but many classes still have quite a few dependencies. However, all these dependencies are included in the package and do not require any special treatment like building or compiling. At the beginning of each file I also indicate which dependencies the corresponding class has.

##### NOTE: This code has been tested only on Mac OS!

### Run the code

In order to run the code do the following:

1. Install [macports](https://www.macports.org/install.php)
2. Open main.cpp and edit the path to the folder with images on the line 142
3. Open terminal and type the following:

```bash
  sudo port install cmake
```
```bash
  cd path_to_your_folder_with_gbc_package/2d/
```
```bash
  mkdir bin
```
```bash
  cd bin
```
```bash
  cmake -DCMAKE_BUILD_TYPE=Debug ..
```
```bash
  make
```
```bash
  ./gbc
```

For the release version use instead: 

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

##### NOTE: If the output result shows FAILED for some tests, it is ok. It happens because the corresponding coordinate functions do not satisfy some theoretical properties of generalized barycentric coordinates.

### Simple example

```C++
// Polygon.
std::vector<VertexR2> poly(4);

poly[0] = VertexR2(0.1, 0.1);
poly[1] = VertexR2(1.0, 0.0);
poly[2] = VertexR2(0.9, 0.9);
poly[3] = VertexR2(0.2, 1.0);

// Evaluation point.
VertexR2 centre(0.5, 0.5);

// Storage for the computed barycentric coordinates.
std::vector<double> b;

// Compute mean value coordinates.
MeanValueR2 mvc(poly);
mvc.compute(centre, b);

// Output the resulting coordinates.
std::cout << "\nMean value coordinates: ";
for (size_t i = 0; i < b.size(); ++i) std::cout << b[i] << " ";
std::cout << "\n\n";
```

##### NOTE: For more examples see main.cpp.

### Visual output

This package allows you to save the coordinate functions in the eps format, for example:

**_Wachspress coordinates_**  
![Wachspress coordinates](2d/examples/WachspressR2.png)

**_Mean value coordinates_**  
![Mean value coordinates](2d/examples/MeanValueR2.png)

**_Metric coordinates_**  
![Metric coordinates](2d/examples/MetricR2.png)

**_Harmonic coordinates_**  
![Harmonic coordinates](2d/examples/HarmonicR2.png)
