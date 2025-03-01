// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include <vector>          // For std::vector used in BFS
#include <unordered_map>   // For tracking visited tetrahedra
#include <algorithm>       // Potentially unused; might be for std::find or similar (not in this code)
#include "vertex.hpp"      // Vertex class definition
#include "tetra.hpp"       // Tetra class definition for tetrahedron access

// Comment: Implements Vertex member function for 3D CDT triangulation (Sec. 3.2.2).

bool Vertex::neighborsVertex(Vertex::Label v) {
    // Comment: Checks if vertex 'v' is a neighbor of this vertex (connected by an edge).
    // Purpose: Used in move validation (e.g., flip move, Sec. 2.3.2) or observable calculations (Sec. 3.4).
    // Method: Breadth-First Search (BFS) starting from 'tetra', exploring tetrahedra containing this vertex.

    Vertex::Label vc = *this;  // Current vertex label (implicit conversion from Vertex to Label)
    // Comment: 'this' is the current Vertex object; cast to Label (likely an int) for comparison.
    if (v == vc) return false;  // Early exit if checking self as neighbor
    // Comment: Prevents counting a vertex as its own neighbor.

    auto t = tetra;  // Starting (3,1)-tetrahedron containing this vertex (Sec. 3.2.2)
    // Comment: 'tetra' is a Pool<Tetra>::Label (Sec. 3.1.1), initialized in Vertex (default -1 if unassigned).

    std::unordered_map<int, bool> done;  // Tracks visited tetrahedra to avoid cycles
    done.reserve(v->cnum);  // Pre-allocate for expected number of tetrahedra (Sec. 3.2.2)
    // Comment: Reserves space for efficiency; 'v->cnum' approximates tetrahedra around 'vc'.
    // HPC Target [General #10]: Pre-allocation aids cache, but dynamic resizing could still occur.

    std::vector<Tetra::Label> current = {t};  // Current BFS frontier, initialized with starting tetra
    std::vector<Tetra::Label> next = {};      // Next BFS frontier
    // Comment: BFS uses two vectors to process tetrahedra level-by-level (Sec. 3.4 BFS description).

    do {
        for (auto tc : current) {  // Iterate over current frontier tetrahedra
            // Comment: 'tc' is a Tetra::Label from Pool<Tetra>.
            for (auto tcn : tc->tnbr) {  // Check neighboring tetrahedra (Sec. 3.2.2, Fig. 7)
                // Comment: 'tnbr' is Tetra’s array of 4 neighboring tetra labels.
                if (!tcn->hasVertex(vc)) continue;  // Skip if neighbor doesn’t contain this vertex
                // Comment: Ensures we stay within tetrahedra containing 'vc' (Sec. 3.2.2 reconstruction).

                if (done.find(tcn) == done.end()) {  // If tetra not yet visited
                    if (tcn->hasVertex(v)) return true;  // Found 'v' as neighbor; exit
                    // Comment: Direct neighbor check; returns true if 'v' is in this tetra.
                    done[tcn] = true;  // Mark tetra as visited
                    next.push_back(tcn);  // Add to next frontier for exploration
                    // Comment: Expands BFS to unvisited tetrahedra containing 'vc'.
                }
            }
        }
        current = next;  // Move to next frontier
        next.clear();    // Reset next frontier
        // Comment: Classic BFS level progression (Sec. 3.4).
    } while (current.size() > 0);  // Continue until no more tetrahedra to explore

    return false;  // 'v' not found as neighbor
    // Comment: Exhaustive search completed; 'v' is not adjacent to 'vc'.
}

// HPC Targets Summary:
// [OpenMP #3]: Parallelize BFS outer loop or use level-synchronous approach for observable speedup.
// [GPU #8]: Port BFS to CUDA kernel for massive parallelism in distance/sphere calculations.
// [General #10]: Optimize cache with pre-allocated contiguous arrays instead of vectors/maps.