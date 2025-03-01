// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "hausdorff2d.hpp"  // Hausdorff2d class definition

// Comment: Implements Hausdorff2d to measure 2D Hausdorff dimension on spatial slices (Sec. 3.4).
void Hausdorff2d::process() {
    std::string tmp = "";  // Output string
    std::vector<int> profile = {};  // Profile of sphere sizes or distances
    int max_epsilon = 30;  // Maximum distance/epsilon for measurement
    profile.resize(max_epsilon, 0);  // Pre-allocate profile up to max_epsilon
    // Comment: max_epsilon shadows class member; likely meant to be configurable.

    int vmax = 0;  // Maximum vertex label
    for (auto v : Universe::vertices) {
        if (v > vmax) vmax = v;  // Find max vertex index
    }
    doneL.resize(vmax + 1, false);  // Resize BFS visitation flags (inherited from Observable)
    // Comment: Prepares doneL for BFS; size based on vertex count (Sec. 3.4).

    if (!average) {  // Non-averaged mode
        for (int i = 1; i <= max_epsilon; i++) {
            Vertex::Label v;
            do {  // Select vertex from slice matching target2Volume
                v = Universe::verticesAll.pick();  // Random vertex (Sec. 3.1.2)
            } while (Universe::sliceSizes[v->time] != Simulation::target2Volume);

            auto s1 = sphere2d(v, i);  // 2D sphere at radius i (Sec. 3.4)
            profile.at(i - 1) = s1.size();  // Store sphere size
        }
    } else if (average) {  // Averaged mode
        printf("avg\n");

        int counter = 0;  // Number of vertices processed
        for (auto v : Universe::verticesAll) {  // Iterate over all vertices
            if (Universe::sliceSizes[v->time] != Simulation::target2Volume) continue;  // Filter by slice
            counter++;

            auto singleProfile = distanceList2d(v);  // Compute distance profile for v
            if (singleProfile.size() > profile.size()) profile.resize(singleProfile.size(), 0);  // Expand if needed

            std::string tmp = "";  // Local tmp (shadows outer tmp; unused in final output)
            for (int i = 0; i < singleProfile.size(); i++) {
                profile.at(i) += singleProfile.at(i);  // Accumulate distances
                tmp += std::to_string(profile.at(i) / counter) + " ";  // Intermediate average (unused)
            }
        }

        for (int i = 0; i < profile.size(); i++) {
            profile.at(i) /= counter;  // Compute final averages
        }
    }
    // HPC Target [OpenMP #2]: Parallelize vertex loop in average mode.

    for (auto d : profile) {  // Format output
        tmp += std::to_string(d) + " ";
    }
    tmp.pop_back();  // Remove trailing space
    output = tmp;    // Set for Observable::write()
    // HPC Target [General #10]: Pre-allocate tmp for cache efficiency.
}

std::vector<int> Hausdorff2d::distanceList2d(Vertex::Label origin) {
    // TODO(JorenB): optimize BFS - Indicates potential improvement
    std::vector<int> dsts;  // Distance profile (sphere sizes per depth)
    std::vector<Vertex::Label> done;  // Visited vertices
    std::vector<Vertex::Label> thisDepth;  // Current BFS depth
    std::vector<Vertex::Label> nextDepth;  // Next BFS depth

    done.push_back(origin);
    thisDepth.push_back(origin);

    int currentDepth = 0;  // Unused; could track depth explicitly
    do {
        for (auto v : thisDepth) {
            for (auto neighbor : Universe::vertexNeighbors[v]) {
                if (neighbor->time != origin->time) continue;  // Restrict to slice
                if (std::find(done.begin(), done.end(), neighbor) == done.end()) {  // Unvisited
                    nextDepth.push_back(neighbor);
                    done.push_back(neighbor);
                }
            }
        }
        dsts.push_back(thisDepth.size());  // Record sphere size at this depth

        thisDepth = nextDepth;
        nextDepth.clear();
        currentDepth++;  // Increment depth (unused in output)
    } while (thisDepth.size() > 0);

    return dsts;  // Comment: Returns sphere sizes up to max distance (Sec. 3.4).
    // HPC Targets [OpenMP #3, GPU #8]: Parallelize BFS for speedup.
}