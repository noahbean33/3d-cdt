// Copyright 2020 Joren Brunekreef and Andrzej Görlich
#include <string>           // For std::string (tmp, output)
#include <vector>           // For std::vector (epsilons, distanceList, etc.)
#include <algorithm>        // For std::accumulate in averageSphereDistance()
#include <unordered_map>    // Unused; possibly vestigial
#include "ricci2d.hpp"      // Ricci2d class definition
#include <chrono>           // For timing (commented out)

// Comment: Implements Ricci2d to measure 2D Ricci curvature on spatial slices (Sec. 3.4).
void Ricci2d::process() {
    epsilons = {};  // Clear epsilons vector
    for (int i = 1; i <= eps_max; i++) {  // Populate epsilon range
        epsilons.push_back(i);
    }
    // Comment: Sets epsilons from 1 to eps_max (e.g., 10 by default).

    std::vector<double> epsilonDistanceList;  // Stores average distances per epsilon
    std::vector<Vertex::Label> origins;       // Origin vertices for each epsilon

    int vmax = 0;  // Maximum vertex label
    for (auto v : Universe::vertices) {
        if (v > vmax) vmax = v;  // Find max vertex index
    }
    doneL.resize(vmax + 1, false);    // Resize BFS visitation flags
    doneLr.resize(vmax + 1, false);   // Resize local BFS flags
    vertexLr.resize(vmax + 1, false); // Resize vertex flags
    // Comment: Prepares vectors for BFS; doneL inherited from Observable (Sec. 3.4).

    for (std::vector<int>::iterator it = epsilons.begin(); it != epsilons.end(); it++) {
        Vertex::Label v;
        do {  // Select vertex from slice matching target2Volume
            v = Universe::verticesAll.pick();  // Random vertex (Sec. 3.1.2)
        } while (Universe::sliceSizes[v->time] != Simulation::target2Volume);
        origins.push_back(v);  // Store origin for this epsilon
    }
    // HPC Target [OpenMP #2]: Could parallelize origin selection with thread-local RNG.

    for (int i = 0; i < epsilons.size(); i++) {  // Process each epsilon
        int epsilon = epsilons[i];
        auto origin = origins[i];
        double averageDistance = averageSphereDistance(origin, epsilon);  // Compute average
        epsilonDistanceList.push_back(averageDistance);
    }
    // HPC Target [OpenMP #2]: Parallelize this loop with private epsilonDistanceList.

    std::string tmp = "";
    for (double dst : epsilonDistanceList) {  // Format output
        tmp += std::to_string(dst) + " ";
    }
    tmp.pop_back();  // Remove trailing space
    output = tmp;    // Set for Observable::write()
    // HPC Target [General #10]: Pre-allocate tmp for cache efficiency.
}

double Ricci2d::averageSphereDistance(Vertex::Label p1, int epsilon) {
    auto s1 = sphere2d(p1, epsilon);  // 2D sphere around p1 (Sec. 3.4)
    if (s1.size() == 0) return 0.0;   // Early exit if empty
    int t1 = p1->time;                // Time slice of p1

    std::uniform_int_distribution<> rv(0, s1.size() - 1);
    auto p2 = s1.at(rv(rng));         // Random vertex from s1
    auto s2 = sphere2d(p2, epsilon);  // 2D sphere around p2
    if (s2.size() == 0) return 0.0;   // Early exit if empty
    int t2 = p2->time;                // Time slice of p2 (should match t1)

    if (s2.size() < s1.size()) {  // Ensure s1 is smaller for efficiency
        std::swap(s1, s2);
    }

    std::vector<int> distanceList;  // Distances from s1 vertices to s2

    for (auto b : s1) {  // BFS from each vertex in s1
        std::fill(doneLr.begin(), doneLr.end(), false);    // Reset BFS flags
        std::fill(vertexLr.begin(), vertexLr.end(), false); // Reset vertex flags
        for (auto v : s2) {
            vertexLr.at(v) = true;  // Mark s2 vertices as targets
        }

        int countdown = s2.size();  // Number of targets to find
        std::vector<Vertex::Label> thisDepth;  // Current BFS depth
        std::vector<Vertex::Label> nextDepth;  // Next BFS depth

        doneLr.at(b) = true;
        thisDepth.push_back(b);

        for (int currentDepth = 0; currentDepth < 3 * epsilon + 1; currentDepth++) {
            for (auto v : thisDepth) {
                if (vertexLr.at(v)) {  // Found a target
                    distanceList.push_back(0);  // Distance 0 (same vertex)
                    vertexLr.at(v) = false;
                    countdown--;
                }

                for (auto neighbor : Universe::vertexNeighbors[v]) {
                    if (neighbor->time != v->time) continue;  // Restrict to slice
                    if (!doneLr.at(neighbor)) {
                        nextDepth.push_back(neighbor);
                        doneLr.at(neighbor) = true;
                        if (vertexLr.at(neighbor)) {  // Found a target
                            distanceList.push_back(currentDepth + 1);
                            vertexLr.at(neighbor) = false;
                            countdown--;
                        }
                    }
                    if (countdown == 0) break;  // All targets found
                }
                if (countdown == 0) break;
            }
            thisDepth = nextDepth;
            nextDepth.clear();
            if (countdown == 0) break;  // Early exit
        }
        assert(countdown == 0);  // Ensure all s2 vertices found
    }

    int distanceSum = std::accumulate(distanceList.begin(), distanceList.end(), 0);
    double averageDistance = static_cast<double>(distanceSum) / (epsilon * distanceList.size());
    return averageDistance;
    // HPC Targets [OpenMP #3, GPU #8]: Parallelize BFS; optimize vector usage.
}