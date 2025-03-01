// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "hausdorff2d_dual.hpp"  // Hausdorff2dDual class definition

// Comment: Implements Hausdorff2dDual to measure 2D Hausdorff dimension in the dual graph (Sec. 3.4).
void Hausdorff2dDual::process() {
    std::string tmp = "";  // Output string

    std::uniform_int_distribution<> triangleGen(0, Universe::triangles.size() - 1);
    Triangle::Label tr;
    do {  // Select a triangle from a slice matching target2Volume
        tr = Universe::triangles.at(triangleGen(rng));  // Random triangle (Sec. 3.2)
    } while (Universe::sliceSizes[tr->time] != Simulation::target2Volume);
    // Comment: Ensures triangle originates from a specific slice (Sec. 2.4).
    // Commented alternative: /* } while (tr->time != 1); */ suggests past hardcoded slice 1.

    auto dsts = distanceList2dDual(tr);  // Compute dual distance profile

    for (auto d : dsts) {  // Format output
        tmp += std::to_string(d) + " ";
    }
    tmp.pop_back();  // Remove trailing space
    output = tmp;    // Set for Observable::write()
    // HPC Target [General #10]: Pre-allocate tmp for cache efficiency.
}

std::vector<int> Hausdorff2dDual::distanceList2dDual(Triangle::Label origin) {
    std::vector<int> dsts;  // Distance profile (dual sphere sizes per depth)
    std::vector<Triangle::Label> done;  // Visited triangles
    std::vector<Triangle::Label> thisDepth;  // Current BFS depth
    std::vector<Triangle::Label> nextDepth;  // Next BFS depth

    done.push_back(origin);
    thisDepth.push_back(origin);

    int currentDepth = 0;  // Unused; could track depth explicitly
    do {
        for (auto t : thisDepth) {
            for (auto neighbor : t->trnbr) {  // Adjacent triangles via Triangle::trnbr
                if (std::find(done.begin(), done.end(), neighbor) == done.end()) {  // Unvisited
                    nextDepth.push_back(neighbor);
                    done.push_back(neighbor);
                }
            }
        }

        dsts.push_back(thisDepth.size());  // Record dual sphere size at this depth

        thisDepth = nextDepth;
        nextDepth.clear();
        currentDepth++;  // Increment depth (unused in output)
    } while (thisDepth.size() > 0);

    return dsts;  // Comment: Returns dual sphere sizes up to max distance (Sec. 3.4).
    // HPC Targets [OpenMP #3, GPU #8]: Parallelize BFS for speedup.
}