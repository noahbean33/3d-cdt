// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include "pool.hpp"  // Pool template for O(1) memory management (Sec. 3.1.1)

// Comment: Forward declarations of related classes.
class Vertex;
class Triangle;
class Tetra;

// Comment: HalfEdge represents a directed edge in the spatial triangulation of 3D CDT (Sec. 3.2).
class HalfEdge : public Pool<HalfEdge> {  // Inherits from Pool for efficient management
// Comment: Manages half-edges, which are edges of spatial triangles in (3,1)-tetrahedra (Sec. 3.2.2).

public:
    static const unsigned pool_size = 5000000;  // Maximum number of half-edges in Pool
    // Comment: Defines capacity C (Sec. 3.1.1), set to 5 million; limits system size.
    // HPC Target [General #10]: Fixed size aids cache but could be dynamic for scalability.

    std::array<Pool<Vertex>::Label, 2> vs = {-1, -1};  // Two vertex labels (start, end)
    // Comment: Represents directed edge endpoints; initialized to -1 (invalid).

    Pool<HalfEdge>::Label adj = -1;   // Adjacent half-edge (opposite direction)
    Pool<HalfEdge>::Label next = -1;  // Next half-edge in triangle (counterclockwise)
    Pool<HalfEdge>::Label prev = -1;  // Previous half-edge in triangle
    // Comment: Connectivity fields linking to other half-edges; -1 indicates unset (Sec. 3.2).

    Pool<Tetra>::Label tetra = -1;    // Containing (3,1)-tetrahedron
    // Comment: Links to tetra using this half-edge in its base triangle (Sec. 3.2.2).

    Pool<Triangle>::Label triangle = -1;  // Containing triangle
    // Comment: Links to spatial triangle formed by this edge (Sec. 3.2.2).

    void setVertices(Pool<Vertex>::Label ve, Pool<Vertex>::Label vf) {
        vs[0] = ve;  // Start vertex
        vs[1] = vf;  // End vertex
    }
    // Comment: Sets edge endpoints; called in updateHalfEdgeData() (Sec. 3.2).

    HalfEdge::Label getAdjacent() {
        return adj;  // Returns adjacent half-edge
    }
    // Comment: Accessor for opposite-directed edge; used in triangle neighbor linking (Sec. 3.2).

    void setAdjacent(HalfEdge::Label he) {
        he->adj = *this;  // Set reciprocal adjacency
        adj = he;
    }
    // Comment: Links this half-edge to its opposite; ensures bidirectional adjacency (Sec. 3.2).

private:
    // Comment: No additional private members; public access simplifies usage but risks unintended modification.
};

// HPC Targets Summary:
// [General #10]: pool_size could be pre-allocated or dynamic; std::array usage is cache-friendly.
// [OpenMP #3, GPU #8]: No direct parallelization here, but usage in BFS (e.g., Observable::sphere2dDual) could benefit.