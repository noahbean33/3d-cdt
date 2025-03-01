// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
// Comment: Indicates authorship and year of creation; part of the 3D CDT codebase.

#pragma once
// Comment: Prevents multiple inclusions of this header file, a standard C++ practice to avoid redefinition errors.

#include <string>
// Comment: Includes the standard string library for handling std::string, used for the id and name fields.

#include <vector>
// Comment: Includes the standard vector library for std::vector, used in the distanceList2d() return type and potentially in process().

#include "../observable.hpp"
// Comment: Includes the base Observable class header, located in a parent directory, providing inheritance functionality for measurement and output.

#include "../universe.hpp"
// Comment: Includes the Universe class header, providing access to the triangulation data (e.g., vertices, sliceSizes) used in calculations.

// Comment: Hausdorff2d is a class that measures the 2D Hausdorff dimension on spatial slices (S^2 topology) within the 3D CDT triangulation.
// The Hausdorff dimension quantifies the fractal scaling of distances within a 2D slice, often through sphere growth or pairwise distances (Sec. 3.4).
class Hausdorff2d : public Observable {  // Inherits publicly from Observable
public:
    using Observable::Observable;
    // Comment: Inherits the constructor(s) of the base Observable class (C++11 feature), allowing Hausdorff2d to use Observable(id) directly without redefining it.
    // This ensures the base class's initialization (e.g., setting identifier) is reused.

    Hausdorff2d(std::string id) : Observable(id) {  // Default constructor
        name = "hausdorff2d";  // Sets the observable's name, used for output file naming (e.g., "data/hausdorff2d-<id>.dat").
        average = false;       // Default behavior: no averaging of distances; likely outputs raw sphere sizes or distances.
    }
    // Comment: Constructs an instance with a given identifier (e.g., from main.cpp) and defaults to non-averaged mode.

    Hausdorff2d(std::string id, bool average_) : Observable(id) {  // Parameterized constructor
        name = "hausdorff2d";  // Consistent naming for output identification.
        average = average_;    // Allows customization of whether to average distances or report raw values.
    }
    // Comment: Alternative constructor providing flexibility to toggle averaging behavior, initialized with an id and average flag from the caller (e.g., main.cpp).

    void process();
    // Comment: Declares the process() method, overridden from Observable, to compute the 2D Hausdorff dimension (Sec. 3.4).
    // Likely calculates sphere sizes or distance distributions within a slice, implemented in hausdorff2d.cpp.
    // HPC Target [OpenMP #2]: If process() involves loops over vertices or slices, it can be parallelized for performance gains (4-16x speedup).

private:
    int max_epsilon;
    // Comment: Defines the maximum distance (epsilon) or radius for Hausdorff measurement; not initialized here, suggesting it's set dynamically in process().
    // Represents the upper bound for distance sampling or sphere growth analysis (Sec. 3.4).

    bool average;
    // Comment: Boolean flag determining whether distances are averaged (e.g., mean Hausdorff dimension) or reported as raw values (e.g., sphere sizes per epsilon).
    // Provides flexibility in how the Hausdorff dimension is presented (Sec. 3.4).

    void initialize() { }
    // Comment: Overrides Observable::initialize() with an empty implementation, indicating no persistent state setup is needed beyond the base class.
    // Suggests that Hausdorff2d computes its measurements on-the-fly using Universe data (Sec. 3.4).

    std::vector<int> distanceList2d(Vertex::Label origin);
    // Comment: Declares a helper method to compute a list of 2D distances from a given vertex (origin) within its spatial slice.
    // Likely uses Observable::sphere2d() or Observable::distance() to gather distances or sphere sizes (Sec. 3.4).
    // Returns a vector of integers, possibly representing sphere sizes at increasing radii or pairwise distances.
    // HPC Target [GPU #8]: If BFS-based (e.g., via sphere2d), this could be GPU-accelerated for significant speedup (10-50x).
};

// HPC Targets Summary:
// [OpenMP #2]: Parallelize process() if it loops over vertices or slices to compute measurements (Sec. 3.4).
//             Potential speedup: 4-16x by distributing work across threads, depending on vertex count or slice size.
// [OpenMP #3]: Parallelize BFS in distanceList2d() if it uses breadth-first search (e.g., sphere2d), offering 2-4x speedup per call (Sec. 3.4).
// [GPU #8]: GPU-accelerate BFS operations in distanceList2d() for large max_epsilon, potentially achieving 10-50x speedup (Sec. 3.4).
// [General #10]: Pre-allocate vectors (e.g., in process() or distanceList2d()) to avoid dynamic resizing, improving cache efficiency (20-50% potential, Sec. 3.1).