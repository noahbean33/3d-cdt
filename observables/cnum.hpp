// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <string>           // For std::string (id, name)
#include "../observable.hpp" // Base Observable class
#include "../universe.hpp"   // Access to Universe triangulation data

// Comment: CNum measures the coordination number (cnum) distribution in 3D CDT (Sec. 3.4).
class CNum : public Observable {  // Inherits from Observable
public:
    CNum(std::string id) : Observable(id) {  // Constructor
        name = "cnum";  // Set observable name
    }
    // Comment: Initializes with identifier (e.g., from main.cpp); sets name for output file (Sec. 3).

    void process();  // Computes coordination number distribution
    // Comment: Overrides Observable::process() to measure cnum per vertex (Sec. 3.4).
    // HPC Target [OpenMP #2]: Potential for parallel processing if loop-based.

private:
    void initialize() {}  // Empty initialization
    // Comment: Overrides Observable::initialize(); no setup needed (Sec. 3.4).
};

// HPC Targets Summary:
// [OpenMP #2]: Parallelize process() if it loops over vertices (Sec. 3.4).
// [General #10]: Pre-allocate vectors in process() for cache efficiency (Sec. 3.1).
// No direct BFS (OpenMP #3, GPU #8) unless process() uses sphere/distance methods.