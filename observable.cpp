// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include <iostream>       // For std::cout (unused here)
#include <fstream>        // For std::ofstream in write(), clear()
#include <vector>         // For std::vector in sphere methods
#include "observable.hpp"  // Observable class definition

// Static member initializations (Sec. 3)
std::default_random_engine Observable::rng(0);  // RNG with default seed 0
// TODO(JorenB): seed properly - Not seeded here; likely set elsewhere (e.g., Simulation::start)
// HPC Target [MPI #4]: Static RNG needs per-rank instances.
std::string Observable::data_dir = "";  // Output directory, set by main.cpp
std::vector<bool> Observable::doneL;    // BFS visitation tracker (Sec. 3.4)
// HPC Target [General #10]: Should be pre-allocated; [OpenMP #3]: Needs thread-local copies.

void Observable::write() {  // Writes output to file (Sec. 3.4)
    std::string filename = data_dir + "/" + name + "-" + identifier + extension;  // e.g., "data/VolumeProfile-fID.dat"
    std::ifstream infile(filename);  // Check if file exists (unused result)

    std::ofstream file;
    file.open(filename, std::ios::app);  // Append mode
    assert(file.is_open());  // Ensure file opened

    file << output << "\n";  // Write output string
    file.close();
    // Comment: Appends computed data from process() to file.
}

void Observable::clear() {  // Resets observable data (Sec. 3.4)
    std::string filename = data_dir + "/" + name + "-" + identifier + extension;
    std::ofstream file;
    file.open(filename, std::ios::app);  // Append mode (no truncation)
    assert(file.is_open());
    file.close();  // No data written; just ensures file exists

    initialize();  // Reset derived class state
    // Comment: Called by Simulation::start() to clear prior runs; file not cleared, just appended.
    // Potential Bug: Should use std::ios::trunc to clear file if reset intended.
}

std::vector<Vertex::Label> Observable::sphere(Vertex::Label origin, int radius) {  // Vertex sphere (Sec. 3.4)
    std::vector<Vertex::Label> thisDepth;  // Current BFS depth
    std::vector<Vertex::Label> nextDepth;  // Next BFS depth
    std::vector<Vertex::Label> vertexList; // Resulting sphere vertices
    std::vector<Vertex::Label> flippedVertices;  // Track visited for reset

    doneL.at(origin) = true;  // Mark origin visited
    thisDepth.push_back(origin);
    flippedVertices.push_back(origin);

    for (int currentDepth = 0; currentDepth < radius; currentDepth++) {  // BFS loop
        for (auto v : thisDepth) {  // Process current depth
            for (auto neighbor : Universe::vertexNeighbors[v]) {  // Check neighbors
                if (!doneL.at(neighbor)) {  // If unvisited
                    flippedVertices.push_back(neighbor);
                    nextDepth.push_back(neighbor);
                    doneL.at(neighbor) = true;
                    if (currentDepth == radius - 1) vertexList.push_back(neighbor);  // Add to result at final depth
                }
            }
        }
        thisDepth = nextDepth;  // Move to next depth
        nextDepth.clear();
    }

    for (auto v : flippedVertices) doneL.at(v) = false;  // Reset visited flags
    return vertexList;
    // Comment: BFS for link distance sphere; returns vertices at radius.
    // HPC Targets [OpenMP #3, GPU #8]: Parallelizable BFS.
}

std::vector<Vertex::Label> Observable::sphere2d(Vertex::Label origin, int radius) {  // 2D vertex sphere (Sec. 3.4)
    std::vector<Vertex::Label> thisDepth;
    std::vector<Vertex::Label> nextDepth;
    std::vector<Vertex::Label> vertexList;
    std::vector<Vertex::Label> flippedVertices;

    doneL.at(origin) = true;
    thisDepth.push_back(origin);
    flippedVertices.push_back(origin);

    for (int currentDepth = 0; currentDepth < radius; currentDepth++) {
        for (auto v : thisDepth) {
            for (auto neighbor : Universe::vertexNeighbors[v]) {
                if (neighbor->time != origin->time) continue;  // Restrict to same time slice
                if (!doneL.at(neighbor)) {
                    flippedVertices.push_back(neighbor);
                    nextDepth.push_back(neighbor);
                    doneL.at(neighbor) = true;
                    if (currentDepth == radius - 1) vertexList.push_back(neighbor);
                }
            }
        }
        thisDepth = nextDepth;
        nextDepth.clear();
    }

    for (auto v : flippedVertices) doneL.at(v) = false;
    return vertexList;
    // Comment: 2D BFS within origin’s time slice; supports 2D observables (e.g., Ricci2d).
    // HPC Targets [OpenMP #3, GPU #8]: Parallelizable BFS.
}

std::vector<Tetra::Label> Observable::sphereDual(Tetra::Label origin, int radius) {  // Dual tetra sphere (Sec. 3.4)
    // TODO(JorenB): optimize BFS (see sphere()) - Indicates potential improvement
    std::vector<Tetra::Label> done;  // Visited tetrahedra
    std::vector<Tetra::Label> thisDepth;
    std::vector<Tetra::Label> nextDepth;

    done.push_back(origin);
    thisDepth.push_back(origin);
    std::vector<Tetra::Label> tetraList;

    for (int currentDepth = 0; currentDepth < radius; currentDepth++) {
        for (auto t : thisDepth) {
            for (auto neighbor : t->tnbr) {  // Check tetra neighbors
                if (std::find(done.begin(), done.end(), neighbor) == done.end()) {  // If unvisited
                    nextDepth.push_back(neighbor);
                    done.push_back(neighbor);
                    if (currentDepth == radius - 1) tetraList.push_back(neighbor);
                }
            }
        }
        thisDepth = nextDepth;
        nextDepth.clear();
    }

    return tetraList;
    // Comment: BFS for dual link distance; doesn’t use doneL (less efficient).
    // HPC Targets [OpenMP #3, GPU #8]: Parallelizable BFS.
}

std::vector<Triangle::Label> Observable::sphere2dDual(Triangle::Label origin, int radius) {  // 2D dual triangle sphere
    std::vector<Triangle::Label> done;
    std::vector<Triangle::Label> thisDepth;
    std::vector<Triangle::Label> nextDepth;

    done.push_back(origin);
    thisDepth.push_back(origin);
    std::vector<Triangle::Label> triangleList;

    for (int currentDepth = 0; currentDepth < radius; currentDepth++) {
        for (auto t : thisDepth) {
            for (auto neighbor : t->trnbr) {  // Check triangle neighbors
                if (std::find(done.begin(), done.end(), neighbor) == done.end()) {
                    nextDepth.push_back(neighbor);
                    done.push_back(neighbor);
                    if (currentDepth == radius - 1) triangleList.push_back(neighbor);
                }
            }
        }
        thisDepth = nextDepth;
        nextDepth.clear();
    }

    return triangleList;
    // Comment: 2D dual BFS; doesn’t use doneL.
    // HPC Targets [OpenMP #3, GPU #8]: Parallelizable BFS.
}