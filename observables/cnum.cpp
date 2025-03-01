// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "cnum.hpp"  // CNum class definition

// Comment: Implements CNum to measure the spatial coordination number (scnum) distribution (Sec. 3.4).
void CNum::process() {
    std::string tmp = "";  // Output string

    std::array<int, 750> histogram;  // Histogram of scnum values
    std::fill(histogram.begin(), histogram.end(), 0);  // Initialize to zero
    // Comment: Fixed-size array for scnum counts, assuming max scnum < 750 (Sec. 3.2.2).

    for (auto v : Universe::vertices) {  // Iterate over vertices
        if (Universe::sliceSizes[v->time] != Simulation::target2Volume) continue;  // Filter by slice
        if (v->scnum > histogram.size() - 1) {  // Check for overflow
            printf("overflow. cnum: %d\n", v->scnum);  // Log overflow
            continue;  // Skip if scnum exceeds histogram size
        }
        histogram.at(v->scnum) += 1;  // Increment bin for spatial coordination number
    }
    // Comment: Builds histogram of scnum for vertices in target slice (Sec. 3.4).

    for (auto h : histogram) {  // Format output
        tmp += std::to_string(h) + " ";
    }
    tmp.pop_back();  // Remove trailing space

    output = tmp;  // Set for Observable::write()
    // Comment: Outputs space-separated histogram (e.g., "0 5 10 ...") (Sec. 3.4).
    // HPC Targets [OpenMP #2, General #10]: Parallelize loop; optimize tmp allocation.
}