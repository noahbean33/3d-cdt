// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <string>           // For std::string (id, name)
#include <vector>           // For std::vector (possibly used in process())
#include "../observable.hpp" // Base Observable class
#include "../universe.hpp"   // Access to Universe triangulation data

// Comment: Minbu measures the minimal bunching (Minbu) observable in 3D CDT (Sec. 3.4).
class Minbu : public Observable {  // Inherits from Observable
public:
    using Observable::Observable;  // Inherits Observable constructors
    // Comment: Allows base class constructor to be used directly (e.g., Observable(id)).

    Minbu(std::string id) : Observable(id) {  // Constructor
        name = "minbu";  // Set observable name
    }
    // Comment: Initializes with identifier (e.g., from main.cpp); sets name for output file (Sec. 3).

    void process();  // Computes the Minbu observable
    // Comment: Overrides Observable::process() to calculate minimal bunching (Sec. 3.4).
    // HPC Target [OpenMP #2]: Potential for parallel processing if loop-based.

private:
    void initialize() {}  // Empty initialization
    // Comment: Overrides Observable::initialize(); no setup needed (Sec. 3.4).
};

// HPC Targets Summary:
// [OpenMP #2]: Parallelize process() if it involves looping over slices or vertices (Sec. 3.4).
// [General #10]: Minimal impact here; cache efficiency depends on process() implementation.
// No direct BFS (OpenMP #3, GPU #8) unless process() uses sphere/distance methods.