// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <string>        // For std::string (name, identifier, output)
#include <algorithm>     // For std::fill in sphere() methods
#include <vector>        // For std::vector in sphere and distance methods
#include "universe.hpp"  // Access to Universe triangulation data
#include "simulation.hpp" // Interaction with Simulation class

// Comment: Observable base class for measuring geometric properties in 3D CDT (Sec. 3.4).
class Observable {
public:
    std::string name;  // Observable name
    // Comment: Public identifier for the observable (e.g., "VolumeProfile").

    Observable(std::string identifier_) {
        identifier = identifier_;  // Set unique identifier
    }
    // Comment: Constructor sets identifier (e.g., file ID from main.cpp, Sec. 3).

    void measure() {
        process();  // Compute observable
        write();    // Output results
    }
    // Comment: Called during simulation to measure and record (Sec. 3.3.3).

    void clear();  // Resets observable data
    // Comment: Clears prior measurements (Sec. 3.4).

    static std::string data_dir;  // Output directory for data files
    // Comment: Set by main.cpp; shared across all observables (Sec. 3).

private:
    std::string identifier;  // Unique ID for this observable instance
    // Comment: Used for file naming (e.g., "<identifier>.dat").

protected:
    static std::default_random_engine rng;  // RNG for observable calculations
    // Comment: Shared RNG; usage unclear without process() implementations.
    // HPC Target [MPI #4]: Needs per-rank instances for distributed runs.

    virtual void process() = 0;  // Pure virtual method to compute observable
    // Comment: Implemented by derived classes (e.g., VolumeProfile, Sec. 3.4).
    virtual void initialize() = 0;  // Pure virtual method to initialize
    // Comment: Sets up observable state; called implicitly or in derived constructors.

    void write();  // Writes output to file
    // Comment: Saves computed data (Sec. 3.4).

    static std::vector<bool> doneL;  // Tracks visited nodes in BFS
    // Comment: Shared across instances for sphere/distance methods (Sec. 3.4).
    // HPC Target [OpenMP #3]: Needs thread-local copies for parallel BFS.

    // Toolbox methods for geometric measurements (Sec. 3.4)
    static std::vector<Vertex::Label> sphere(Vertex::Label origin, int radius);
    // Comment: Computes vertex sphere of radius from origin (link distance).
    static std::vector<Vertex::Label> sphere2d(Vertex::Label origin, int radius);
    // Comment: 2D vertex sphere (within a slice).
    static std::vector<Tetra::Label> sphereDual(Tetra::Label origin, int radius);
    // Comment: Dual tetrahedron sphere (dual link distance).
    static std::vector<Triangle::Label> sphere2dDual(Triangle::Label origin, int radius);
    // Comment: 2D dual triangle sphere.

    static int distance(Vertex::Label v1, Vertex::Label v2);
    // Comment: Shortest link distance between vertices.
    static int distanceDual(Tetra::Label t1, Tetra::Label t2);
    // Comment: Shortest dual link distance between tetrahedra.
    // HPC Target [OpenMP #3, GPU #8]: BFS methods are parallelizable.

    std::string extension = ".dat";  // File extension for output
    std::string output;             // Data to write
    // Comment: Output string built by process(), saved by write().
};

// HPC Targets Summary:
// [OpenMP #3]: Parallelize BFS in sphere/distance methods (Sec. 3.4).
// [GPU #8]: GPU-accelerate BFS for large radii (Sec. 3.4).
// [MPI #4]: Distribute rng for parallel runs (Sec. 2).
// [General #10]: Pre-allocate doneL and output buffers (Sec. 3.1).