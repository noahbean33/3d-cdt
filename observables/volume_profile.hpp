// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <string>           // For std::string (id, name)
#include "../observable.hpp" // Base Observable class
#include "../universe.hpp"   // Access to Universe triangulation data

// Comment: VolumeProfile measures the volume distribution across time slices in 3D CDT (Sec. 3.4).
class VolumeProfile : public Observable {  // Inherits from Observable
public:
    VolumeProfile(std::string id) : Observable(id) {  // Constructor
        name = "volume_profile";  // Set observable name
    }
    // Comment: Initializes with identifier (e.g., from main.cpp); sets name for output file (Sec. 3).

    void process();  // Computes volume profile
    // Comment: Overrides Observable::process() to calculate volume per time slice (Sec. 3.4).
    // HPC Target [OpenMP #2]: Potential for parallel measurement.

private:
    void initialize() {}  // Empty initialization
    // Comment: Overrides Observable::initialize(); no setup needed (Sec. 3.4).
};

// HPC Targets Summary:
// [OpenMP #2]: Parallelize process() if loop-based (Sec. 3.4).
// [General #10]: Minimal impact here; cache efficiency depends on process() implementation.
// No direct BFS (OpenMP #3, GPU #8) unless process() uses sphere/distance methods.