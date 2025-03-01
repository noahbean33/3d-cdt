// Copyright 2020 Joren Brunekreef and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard; earlier copyright (2020) vs. others (2021).

#include <string>           // For std::string (id, name)
#include <vector>           // For std::vector (epsilons, doneLr, triangleLr)
#include "../observable.hpp" // Base Observable class
#include "../universe.hpp"   // Access to Universe triangulation data

// Comment: Ricci2dDual measures 2D Ricci curvature in the dual graph of spatial slices (Sec. 3.4).
class Ricci2dDual : public Observable {  // Inherits from Observable
public:
    Ricci2dDual(std::string id) : Observable(id) {  // Default constructor
        name = "ricci2d_dual";  // Set observable name
        eps_max = 10;           // Default maximum epsilon (radius)
    }

    Ricci2dDual(std::string id, int eps_max_) : Observable(id) {  // Parameterized constructor
        name = "ricci2d_dual";  // Set observable name
        eps_max = eps_max_;     // Set custom maximum epsilon
    }
    // Comment: Constructors initialize with identifier (e.g., from main.cpp) and set epsilon range (Sec. 3).

    void process();  // Computes 2D Ricci curvature in dual graph
    // Comment: Overrides Observable::process() to calculate dual curvature per triangle (Sec. 3.4).
    // HPC Target [OpenMP #2]: Potential for parallel triangle processing.

private:
    int eps_max;  // Maximum radius for dual sphere measurements
    // Comment: Defines upper bound for epsilon in curvature calculation (Sec. 3.4).

    std::vector<int> epsilons;  // List of epsilon values computed
    // Comment: Stores radii or distances for dual curvature (Sec. 3.4).

    std::vector<bool> doneLr;  // Visitation flags for dual BFS
    // Comment: Tracks visited triangles or tetrahedra in BFS (Sec. 3.4).
    // HPC Target [OpenMP #3]: Needs thread-local copies for parallel BFS.

    std::vector<bool> triangleLr;  // Triangle-specific flags
    // Comment: Likely tracks processed triangles in dual graph (Sec. 3.4).
    // HPC Target [OpenMP #2]: Thread safety needed if updated in parallel.

    double averageSphereDistanceDual(Triangle::Label p1, int epsilon);  // Computes dual average distance
    // Comment: Helper method for dual curvature; uses Observable::sphere2dDual (Sec. 3.4).
    // HPC Target [GPU #8]: GPU-acceleratable if BFS-based.

    void initialize() {}  // Empty initialization
    // Comment: Overrides Observable::initialize(); no setup needed (Sec. 3.4).
};

// HPC Targets Summary:
// [OpenMP #2]: Parallelize process() for triangle loops (Sec. 3.4).
// [OpenMP #3]: Parallelize BFS in averageSphereDistanceDual() (Sec. 3.4).
// [GPU #8]: GPU-accelerate BFS for large epsilon (Sec. 3.4).
// [General #10]: Pre-allocate epsilons, doneLr, triangleLr (Sec. 3.1).