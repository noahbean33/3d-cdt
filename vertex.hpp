// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Ensures header is included only once during compilation
// Comment: Standard header guard to prevent multiple inclusions.

#include "pool.hpp"  // Includes Pool class template for memory management (Sec. 3.1.1 of paper)
// Comment: Pool provides O(1) operations for simplex storage and access.

class Tetra;  // Forward declaration of Tetra class
// Comment: Avoids circular dependency; Tetra is defined elsewhere (likely tetra.hpp).

class Vertex : public Pool<Vertex> {  // Vertex inherits from Pool for memory management
// Comment: Represents a vertex in the 3D CDT triangulation (Sec. 3.2.2). Inherits from Pool to manage instances efficiently (Sec. 3.1.1).

public:
    static const unsigned pool_size = 3000000;  // Maximum number of vertices in Pool
    // Comment: Defines capacity C for Pool (Sec. 3.1.1), set to 3 million vertices. Limits memory usage but caps system size.
    // HPC Target [General #10]: Fixed size could be pre-allocated for cache efficiency, but consider dynamic sizing for larger systems.

    int time;  // Slice number
    // Comment: Discrete time label τ for the spatial slice containing this vertex (Sec. 2.3). For S^1 x S^2 topology, vertices are ordered by time.

    int scnum;  // Spatial coordination number
    // Comment: Number of spatial (spacelike) neighbors in the same slice (Sec. 2.3.2). Used for flip move checks (Sec. 2.3.2, Fig. 4).

    int cnum;  // Coordination number
    // Comment: Total number of tetrahedra containing this vertex (Sec. 3.2.2). Updated continuously for O(1) delete move checks (Sec. 2.3.1).
    // HPC Note: Continuous updates incur overhead but avoid costly on-the-fly reconstruction (Sec. 3.2.2 footnote 11).

    Pool<Tetra>::Label tetra = -1;  // Some 31-simplex containing this vertex in its base
    // Comment: Label of an arbitrary (3,1)-tetrahedron with this vertex in its base (not apex) (Sec. 3.2.2, Fig. 7). Default -1 indicates unassigned.
    // Purpose: Starting point for neighborhood reconstruction (Sec. 3.2.2). Label type is an alias from Pool<Tetra>.

    bool neighborsVertex(Vertex::Label v);
    // Comment: Checks if vertex 'v' is a neighbor of this vertex (connected by an edge).
    // Usage: Likely used in move validation (e.g., flip move, Sec. 2.3.2) or observable calculations (Sec. 3.4).
    // HPC Target [OpenMP #3, GPU #8]: If used in BFS (e.g., sphere/distance), parallelization could speed up graph traversal.

private:
    // Comment: No private members currently; public access simplifies implementation but risks unintended modification.
    // Suggestion: Future versions might encapsulate reconstruction logic here.
};

// HPC Targets Summary:
// [General #10]: pool_size could be pre-allocated or dynamically sized for cache efficiency.
// [OpenMP #3, GPU #8]: neighborsVertex() could be part of BFS parallelization for observables (Sec. 3.4).