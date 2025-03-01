// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <string>           // For std::string (id, name)
#include <vector>           // For std::vector (distanceList2dDual return, possible process() use)
#include "../observable.hpp" // Base Observable class
#include "../universe.hpp"   // Access to Universe triangulation data

// Comment: Hausdorff2dDual measures the 2D Hausdorff dimension in the dual graph of spatial slices (Sec. 3.4).
class Hausdorff2dDual : public Observable {  // Inherits from Observable
public:
    using Observable::Observable;  // Inherits Observable constructors
    // Comment: Allows base class constructor to be used directly (e.g., Observable(id)).

    Hausdorff2dDual(std::string id) : Observable(id) {  // Constructor
        name = "hausdorff2d_dual";  // Set observable name
    }
    // Comment: Initializes with identifier (e.g., from main.cpp); no averaging option (Sec. 3).

    void process();  // Computes 2D Hausdorff dimension in dual graph
    // Comment: Overrides Observable::process() to measure dual Hausdorff dimension (Sec. 3.4).
    // HPC Target [OpenMP #2]: Potential for parallel processing if loop-based.

private:
    int max_epsilon;  // Maximum distance for dual Hausdorff measurement
    // Comment: Defines upper bound for distance sampling in dual graph (Sec. 3.4).

    void initialize() {}  // Empty initialization
    // Comment: Overrides Observable::initialize(); no setup needed (Sec. 3.4).

    std::vector<int> distanceList2dDual(Triangle::Label origin);  // Computes dual 2D distances
    // Comment: Helper method; likely uses Observable::sphere2dDual (Sec. 3.4).
    // HPC Target [GPU #8]: GPU-acceleratable if BFS-based.
};

// HPC Targets Summary:
// [OpenMP #2]: Parallelize process() if it loops over triangles or slices (Sec. 3.4).
// [OpenMP #3]: Parallelize BFS in distanceList2dDual() if applicable (Sec. 3.4).
// [GPU #8]: GPU-accelerate BFS for large max_epsilon (Sec. 3.4).
// [General #10]: Pre-allocate vectors in process() or distanceList2dDual() (Sec. 3.1).