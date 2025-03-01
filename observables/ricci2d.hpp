// Copyright 2020 Joren Brunekreef and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard; note earlier copyright year (2020) vs. others (2021).

#include <string>           // For std::string (id, name)
#include <vector>           // For std::vector (epsilons, doneLr, vertexLr)
#include "../observable.hpp" // Base Observable class
#include "../universe.hpp"   // Access to Universe triangulation data

// Comment: Ricci2d measures 2D Ricci curvature on spatial slices in 3D CDT (Sec. 3.4).
class Ricci2d : public Observable {  // Inherits from Observable
public:
    Ricci2d(std::string id) : Observable(id) {  // Default constructor
        name = "ricci2d";  // Set observable name
        eps_max = 10;      // Default maximum epsilon (radius)
    }

    Ricci2d(std::string id, int eps_max_) : Observable(id) {  // Parameterized constructor
        name = "ricci2d";  // Set observable name
        eps_max = eps_max_; // Set custom maximum epsilon
    }
    // Comment: Constructors initialize with identifier (e.g., from main.cpp) and set epsilon range (Sec. 3).

    void process();  // Computes 2D Ricci curvature
    // Comment: Overrides Observable::process() to calculate curvature per vertex (Sec. 3.4).
    // HPC Target [OpenMP #2]: Potential for parallel vertex processing.

private:
    int eps_max;  // Maximum radius for sphere measurements
    // Comment: Defines upper bound for epsilon in curvature calculation (Sec. 3.4).

    std::vector<int> epsilons;  // List of epsilon values computed
    // Comment: Stores computed sphere distances or radii used (Sec. 3.4).

    std::vector<bool> doneLr;  // Visitation flags for BFS (likely per vertex or tetra)
    // Comment: Tracks visited nodes; specific use unclear without process().
    // HPC Target [OpenMP #3]: Needs thread-local copies for parallel BFS.

    std::vector<bool> vertexLr;  // Vertex-specific flags
    // Comment: Likely tracks processed vertices; purpose depends on process().
    // HPC Target [OpenMP #2]: Thread safety needed if updated in parallel.

    double averageSphereDistance(Vertex::Label p1, int epsilon);  // Computes average distance
    // Comment: Helper method for curvature; likely uses Observable::sphere2d (Sec. 3.4).
    // HPC Target [GPU #8]: GPU-acceleratable if BFS-based.

    void initialize() {}  // Empty initialization
    // Comment: Overrides Observable::initialize(); no setup needed (Sec. 3.4).
};

// HPC Targets Summary:
// [OpenMP #2]: Parallelize process() for vertex loops (Sec. 3.4).
// [OpenMP #3]: Parallelize BFS in averageSphereDistance() (Sec. 3.4).
// [GPU #8]: GPU-accelerate BFS for large epsilon (Sec. 3.4).
// [General #10]: Pre-allocate epsilons, doneLr, vertexLr (Sec. 3.1).