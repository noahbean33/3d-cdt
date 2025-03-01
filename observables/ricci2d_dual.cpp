// Copyright 2020 Joren Brunekreef and Andrzej Görlich
#include <string>           // For std::string (tmp, output)
#include <vector>           // For std::vector (epsilons, distanceList, etc.)
#include <algorithm>        // For std::accumulate in averageSphereDistanceDual()
#include <unordered_map>    // Unused; possibly vestigial
#include "ricci2d_dual.hpp" // Ricci2dDual class definition
#include <chrono>           // For timing (commented out)

// Comment: Implements Ricci2dDual to measure 2D Ricci curvature on the dual graph of spatial slices (Sec. 3.4).
void Ricci2dDual::process() {
    epsilons = {};  // Clear epsilons vector
    for (int i = 1; i <= eps_max; i++) {  // Populate epsilon range
        epsilons.push_back(i);
    }
    // Comment: Sets epsilons from 1 to eps_max (e.g., 10 by default).

    std::vector<double> epsilonDistanceList;  // Stores average distances per epsilon
    std::vector<Triangle::Label> origins;     // Origin triangles for each epsilon

    int tmax = 0;  // Maximum triangle label
    for (auto t : Universe::triangles) {
        if (t > tmax) tmax = t;  // Find max triangle index
    }
    doneL.resize(tmax + 1, false);      // Resize BFS visitation flags (inherited)
    doneLr.resize(tmax + 1, false);     // Resize local BFS flags
    triangleLr.resize(tmax + 1, false); // Resize triangle flags
    // Comment: Prepares vectors for dual BFS; size based on triangle count (Sec. 3.4).

    std::uniform_int_distribution<> rt(0, Universe::triangles.size() - 1);
    for (std::vector<int>::iterator it = epsilons.begin(); it != epsilons.end(); it++) {
        Triangle::Label t;
        do {  // Select triangle from slice matching target2Volume
            t = Universe::triangles.at(rt(rng));  // Random triangle (Sec. 3.2)
        } while (Universe::sliceSizes[t->time] != Simulation::target2Volume);
        origins.push_back(t);  // Store origin for this epsilon
    }
    // HPC Target [OpenMP #2]: Could parallelize origin selection with thread-local RNG.

    for (int i = 0; i < epsilons.size(); i++) {  // Process each epsilon
        int epsilon = epsilons[i];
        auto origin = origins[i];
        double averageDistance = averageSphereDistanceDual(origin, epsilon);  // Compute average
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

double Ricci2dDual::averageSphereDistanceDual(Triangle::Label p1, int epsilon) {
    auto s1 = sphere2dDual(p1, epsilon);  // Dual 2D sphere around p1 (Sec. 3.4)
    if (s1.size() == 0) return 0.0;       // Early exit if empty
    int t1 = p1->time;                    // Time slice of p1

    std::uniform_int_distribution<> rt(0, s1.size() - 1);
    auto p2 = s1.at(rt(rng));             // Random triangle from s1
    auto s2 = sphere2dDual(p2, epsilon);  // Dual 2D sphere around p2
    if (s2.size() == 0) return 0.0;       // Early exit if empty
    int t2 = p2->time;                    // Time slice of p2 (should match t1)

    if (s2.size() < s1.size()) {  // Ensure s1 is smaller for efficiency
        std::swap(s1, s2);
    }

    std::vector<int> distanceList;  // Distances from s1 triangles to s2

    for (auto b : s1) {  // BFS from each triangle in s1
        std::fill(doneLr.begin(), doneLr.end(), false);    // Reset BFS flags
        std::fill(triangleLr.begin(), triangleLr.end(), false); // Reset triangle flags
        for (auto t : s2) {
            triangleLr.at(t) = true;  // Mark s2 triangles as targets
        }

        int countdown = s2.size();  // Number of targets to find
        std::vector<Triangle::Label> thisDepth;  // Current BFS depth
        std::vector<Triangle::Label> nextDepth;  // Next BFS depth

        doneLr.at(b) = true;
        thisDepth.push_back(b);

        for (int currentDepth = 0; currentDepth < 3 * epsilon + 1; currentDepth++) {
            for (auto t : thisDepth) {
                if (triangleLr.at(t)) {  // Found a target
                    distanceList.push_back(0);  // Distance 0 (same triangle)
                    triangleLr.at(t) = false;
                    countdown--;
                }

                for (auto neighbor : Universe::triangleNeighbors[t]) {
                    if (neighbor->time != t->time) continue;  // Restrict to slice
                    if (!doneLr.at(neighbor)) {
                        nextDepth.push_back(neighbor);
                        doneLr.at(neighbor) = true;
                        if (triangleLr.at(neighbor)) {  // Found a target
                            distanceList.push_back(currentDepth + 1);
                            triangleLr.at(neighbor) = false;
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
        assert(countdown == 0);  // Ensure all s2 triangles found
    }

    int distanceSum = std::accumulate(distanceList.begin(), distanceList.end(), 0);
    double averageDistance = static_cast<double>(distanceSum) / (epsilon * distanceList.size());
    return averageDistance;
    // HPC Targets [OpenMP #3, GPU #8]: Parallelize BFS; optimize vector usage.
}