// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include "pool.hpp"     // Pool template for O(1) memory management (Sec. 3.1.1)
#include "vertex.hpp"   // Vertex class for vertex references
#include "halfedge.hpp" // HalfEdge class for edge connectivity

// Comment: Triangle represents a 2-simplex in the 3D CDT triangulation (Sec. 3.2).

class Triangle : public Pool<Triangle> {  // Inherits from Pool for efficient management
// Comment: Manages triangles (spatial faces of (3,1)-tetrahedra, Sec. 3.2.2) with O(1) operations (Sec. 3.1.1).

public:
    static const unsigned pool_size = 1000000;  // Maximum number of triangles in Pool
    // Comment: Defines capacity C for Pool (Sec. 3.1.1), set to 1 million. Limits system size.
    // HPC Target [General #10]: Fixed size aids cache but could be dynamic for scalability.

    int time;  // Slice number
    // Comment: Time slice (S^1 direction, Sec. 2.3) of the spatial S^2 slice containing this triangle.

    void setVertices(Vertex::Label v0, Vertex::Label v1, Vertex::Label v2) {
        vs = {v0, v1, v2};  // Set triangle vertices
        // Comment: Assigns three vertex labels (Sec. 3.2.2).

        assert(v0->time == v1->time && v0->time == v2->time);  // Validate same time slice
        time = v0->time;  // Set triangle’s time
        // Comment: Ensures triangle lies within one spatial slice (S^2, Sec. 2.3).
    }
    // Comment: Called during geometry update (Sec. 3.2, updateTriangleData).

    void setHalfEdges(HalfEdge::Label h0, HalfEdge::Label h1, HalfEdge::Label h2) {
        hes = {h0, h1, h2};  // Set half-edge labels
        // Comment: Links triangle to its three bounding half-edges (Sec. 3.2).
    }
    // Comment: Establishes connectivity with edges (updateTriangleData).

    void setTriangleNeighbors(Triangle::Label tr0, Triangle::Label tr1, Triangle::Label tr2) {
        trnbr = {tr0, tr1, tr2};  // Set neighboring triangle labels
        // Comment: Links to adjacent triangles via shared edges (Sec. 3.2.2).
    }
    // Comment: Completes adjacency for measurements (Sec. 3.4).

    bool hasVertex(Vertex::Label v) {
        for (int i = 0; i < 3; i++) {  // Check if vertex is in triangle
            if (vs[i] == v) return true;
        }
        return false;
    }
    // Comment: Utility to test vertex inclusion, used in connectivity checks or moves (e.g., Sec. 2.3).

    std::array<Vertex::Label, 3> vs;     // Three vertex labels
    std::array<HalfEdge::Label, 3> hes;  // Three half-edge labels
    std::array<Triangle::Label, 3> trnbr; // Three neighboring triangle labels
    // Comment: Fixed-size arrays store connectivity data (Sec. 3.2.2, Fig. 7 analogue).
    // HPC Target [General #10]: Contiguous storage aids cache efficiency.

private:
    // Comment: No private members; public access simplifies use but risks modification.
};

// HPC Targets Summary:
// [General #10]: pool_size could be pre-allocated or dynamic; arrays already cache-friendly.
// [OpenMP #3, GPU #8]: hasVertex() could be part of parallel BFS for observables (Sec. 3.4).