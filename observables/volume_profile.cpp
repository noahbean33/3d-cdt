// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "volume_profile.hpp"  // VolumeProfile class definition

// Comment: Implements the process() method to compute the volume profile across time slices (Sec. 3.4).
void VolumeProfile::process() {
    std::string tmp = "";  // Temporary string to build output
    for (auto l : Universe::sliceSizes) {  // Iterate over time slices
        tmp += std::to_string(l);  // Append volume (vertices per slice)
        tmp += " ";                // Space separator
    }
    tmp.pop_back();  // Remove trailing space

    output = tmp;  // Set output string for Observable::write()
    // Comment: Generates a space-separated list of slice volumes (Sec. 3.4).
    // HPC Target [OpenMP #2]: Parallelizable loop with thread-local string buffers.
    // HPC Target [General #10]: Pre-allocating tmp could improve cache efficiency.
}