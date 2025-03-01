// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <iostream>      // For printf in log()
#include <typeinfo>      // Unused; possibly for debugging (not in this code)
#include "pool.hpp"      // Pool template for O(1) memory management (Sec. 3.1.1)
#include "vertex.hpp"    // Vertex class for vertex references
#include "halfedge.hpp"  // HalfEdge class for edge connectivity

// Comment: Tetra represents a 3-simplex (tetrahedron) in 3D CDT (Sec. 2.3, 3.2.2).

class Tetra : public Pool<Tetra> {  // Inherits from Pool for efficient management
// Comment: Manages tetrahedra, the fundamental building blocks of 3D CDT (Sec. 2.3).

public:
    static const unsigned pool_size = 5000000;  // Maximum number of tetrahedra in Pool
    // Comment: Defines capacity C (Sec. 3.1.1), set to 5 million. Limits system size.
    // HPC Target [General #10]: Fixed size aids cache but could be dynamic for larger N_3.

    enum Type { THREEONE, ONETHREE, TWOTWO };  // Tetrahedron types (Sec. 2.3)
    // Comment: THREEONE: (3,1), ONETHREE: (1,3), TWOTWO: (2,2) based on vertex time distribution.
    int time;  // Slab number
    // Comment: Time slice (S^1 direction, Sec. 2.3) of the tetra’s "lower" vertices.

    Type type;  // Stores tetrahedron type
    // Comment: Indicates (3,1), (1,3), or (2,2) configuration (Sec. 2.3).

    inline const char* ToString(Tetra::Type t) {  // Converts type to string
        switch (t) {
            case THREEONE: return "31";
            case ONETHREE: return "13";
            case TWOTWO: return "22";
        }
    }
    // Comment: Utility for logging; returns shorthand for tetra type (e.g., "31" for (3,1)).

    void setVertices(Pool<Vertex>::Label v0, Pool<Vertex>::Label v1, Pool<Vertex>::Label v2, Pool<Vertex>::Label v3) {
        // Comment: Sets four vertex labels and determines tetra type (Sec. 3.2.2).
        if (v0->time == v1->time && v0->time == v2->time) type = THREEONE;  // (3,1): 3 vertices at same time
        if (v1->time == v2->time && v1->time == v3->time) type = ONETHREE;  // (1,3): 3 vertices at same time
        if (v0->time == v1->time && v2->time == v3->time) type = TWOTWO;    // (2,2): 2 pairs at same time
        assert(v0->time != v3->time);  // Ensure tetra spans time slices

        vs = {v0, v1, v2, v3};  // Assign vertices
        time = v0->time;  // Set slab time (typically lower time)
    }
    // Comment: Called during initialization or moves (e.g., Sec. 2.3.1 move26).

    void setHalfEdges(Pool<HalfEdge>::Label h0, Pool<HalfEdge>::Label h1, Pool<HalfEdge>::Label h2) {
        hes = {h0, h1, h2};  // Set half-edge labels for (3,1)-tetra base
    }
    // Comment: Links to spatial triangle edges (Sec. 3.2.2); only for (3,1)-tetras.

    HalfEdge::Label getHalfEdgeFrom(Vertex::Label v) {
        for (int i = 0; i < 3; i++) {  // Find half-edge starting at v
            if (hes[i]->vs[0] == v) return hes[i];
        }
        return false;  // Returns -1 if not found (assuming Label is int)
    }
    // Comment: Used in connectivity (e.g., updateHalfEdgeData, Sec. 3.2).

    HalfEdge::Label getHalfEdgeTo(Vertex::Label v) {
        for (int i = 0; i < 3; i++) {  // Find half-edge ending at v
            if (hes[i]->vs[1] == v) return hes[i];
        }
        return false;  // Returns -1 if not found
    }
    // Comment: Supports edge traversal (e.g., Sec. 3.2).

    void setTetras(Pool<Tetra>::Label t0, Pool<Tetra>::Label t1, Pool<Tetra>::Label t2, Pool<Tetra>::Label t3) {
        tnbr = {t0, t1, t2, t3};  // Set neighboring tetra labels
    }
    // Comment: Defines adjacency (Sec. 2.3, Fig. 3-5); updated in moves.

    bool is31() { return type == THREEONE; }  // Checks if (3,1)-tetra
    bool is13() { return type == ONETHREE; }  // Checks if (1,3)-tetra
    bool is22() { return type == TWOTWO; }    // Checks if (2,2)-tetra
    // Comment: Type queries for move conditions (e.g., Sec. 2.3.1).

    bool hasVertex(Pool<Vertex>::Label v) {
        for (int i = 0; i < 4; i++) {  // Check if v is a vertex
            if (vs[i] == v) return true;
        }
        return false;
    }
    // Comment: Used in connectivity and move validation (e.g., Sec. 2.3).
    // HPC Target [OpenMP #3, GPU #8]: Could be part of parallel BFS.

    bool neighborsTetra(Pool<Tetra>::Label t) {
        for (int i = 0; i < 4; i++) {  // Check if t is a neighbor
            if (tnbr[i] == t) return true;
        }
        return false;
    }
    // Comment: Validates adjacency (e.g., in check(), Sec. 1.3).

    Tetra::Label getTetraOpposite(Vertex::Label v) {
        assert(hasVertex(v));  // Ensure v is in tetra
        for (int i = 0; i < 4; i++) {
            if (vs[i] == v) return tnbr[i];  // Return neighbor opposite v
        }
        assert(false);  // Should never reach here
    }
    // Comment: Finds neighbor sharing face opposite v (Sec. 2.3 moves).

    Vertex::Label getVertexOpposite(Vertex::Label v) {
        auto tn = getTetraOpposite(v);  // Neighbor opposite v
        std::array<Vertex::Label, 3> face;  // Face shared with tn
        int i = 0;
        for (auto tv : vs) {  // Build face excluding v
            if (tv != v) face[i++] = tv;
        }
        for (auto tnv : tn->vs) {  // Find vertex not in face
            if (tnv != face[0] && tnv != face[1] && tnv != face[2]) return tnv;
        }
        assert(false);  // Should never reach here
    }
    // Comment: Gets vertex opposite v via neighbor (Sec. 2.3 moves).

    Vertex::Label getVertexOppositeTetra(Tetra::Label tn) {
        for (int i = 0; i < 4; i++) {  // Find vertex opposite neighbor tn
            if (tnbr[i] == tn) return vs[i];
        }
        assert(false);  // Neighbor not found
    }
    // Comment: Reverse lookup for moves (e.g., Sec. 2.3.2 flip).

    void exchangeTetraOpposite(Vertex::Label v, Tetra::Label tn) {
        for (int i = 0; i < 4; i++) {  // Replace neighbor opposite v
            if (vs[i] == v) tnbr[i] = tn;
        }
    }
    // Comment: Updates connectivity during moves (e.g., Sec. 2.3.1).

    void log() {  // Debug logging
        Pool<Tetra>::Label t = *this;  // Current tetra label
        printf("t: %d - %s\n", static_cast<int>(t), ToString(type));
        printf("\t");
        for (int i = 0; i < 4; i++) printf("v%d: %d ", i, vs[i]);  // Print vertices
        printf("\n\t");
        for (int i = 0; i < 4; i++) printf("t%d: %d ", i, tnbr[i]);  // Print neighbors
        printf("\n");
    }
    // Comment: Diagnostic tool for tetra state (Sec. 1.3 check).

    // TODO(JorenB): make private
    std::array<Pool<Tetra>::Label, 4> tnbr;  // Four neighboring tetra labels
    std::array<Pool<Vertex>::Label, 4> vs;   // Four vertex labels
    std::array<Pool<HalfEdge>::Label, 3> hes = {-1, -1, -1};  // Three half-edge labels (for (3,1))
    // Comment: Core data (Sec. 3.2.2, Fig. 7). Public access noted for refactoring.
    // HPC Target [General #10]: Contiguous arrays optimize cache.
};

// HPC Targets Summary:
// [General #10]: pool_size could be dynamic; arrays are cache-friendly.
// [OpenMP #3, GPU #8]: hasVertex(), neighborsTetra() could be parallelized in BFS (Sec. 3.4).