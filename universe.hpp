// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <vector>      // For dynamic arrays (e.g., slabSizes)
#include <random>      // For std::default_random_engine
#include <string>      // For std::string (e.g., fID)
#include "vertex.hpp"  // Vertex class for vertices
#include "halfedge.hpp" // HalfEdge class (likely for connectivity)
#include "triangle.hpp" // Triangle class (sub-simplices)
#include "tetra.hpp"   // Tetra class for tetrahedra
#include "pool.hpp"    // Pool template for O(1) memory management (Sec. 3.1.1)
#include "bag.hpp"     // Bag template for random selection (Sec. 3.1.2)

// Comment: Universe manages the 3D CDT triangulation state (Sec. 3 of paper).

class Universe {
public:
    static int nSlices;  // Number of spatial slices (S^1 direction, Sec. 2.3)
    // Comment: Tracks time slices in S^1 x S^2 topology.

    static std::vector<int> slabSizes;  // Sizes of slabs (tetrahedra between slices)
    static std::vector<int> sliceSizes; // Sizes of spatial slices (S^2 topology)
    // Comment: Store geometric data for volume profile (Sec. 3.4).

    static std::string fID;      // File identifier for output
    static std::string OutFile;  // Output geometry file name
    // Comment: Used for I/O (Sec. 3.3).

    static int strictness;      // Manifold strictness level (Sec. 1.3)
    static int volfix_switch;   // Toggle for volume fixing (Sec. 2.4)
    // Comment: Configuration parameters from main.cpp.

    static Bag<Tetra, Tetra::pool_size> tetrasAll;  // All tetrahedra
    static Bag<Tetra, Tetra::pool_size> tetras31;  // All (3,1)-tetrahedra
    static Bag<Vertex, Vertex::pool_size> verticesAll;  // All vertices
    static Bag<Vertex, Vertex::pool_size> verticesSix;  // Vertices with six tetrahedra
    // Comment: Bags enable O(1) random selection (Sec. 3.1.2). 'tetras31' for (2,6)-move, 'verticesSix' for (6,2)-move (Sec. 2.3.1).
    // HPC Target [General #10]: Large memory footprint; optimize with contiguous storage.

    static bool initialize(std::string geometryFilename, std::string fID, int strictness, int volfix_switch);
    // Comment: Sets up initial triangulation (Sec. 3.1), optionally loading from file. Returns success status.
    // HPC Target [General #10]: Pre-allocate Pools/Bags for cache efficiency.

    static bool exportGeometry(std::string geometryFilename);
    // Comment: Saves triangulation to file (Sec. 3). Returns success status.

    static bool move26(Tetra::Label t);  // (2,6)-move (add move, Sec. 2.3.1)
    static bool move62(Vertex::Label v); // (6,2)-move (delete move, Sec. 2.3.1)
    // Comment: Add/delete moves modify volume (Fig. 3). Return success status.
    // HPC Target [GPU #7]: Batch move attempts on GPU for speedup.

    static bool move44(Tetra::Label t123, Tetra::Label t134);  // (4,4)-move (flip move, Sec. 2.3.2)
    // Comment: Flips spatial edge (Fig. 4). Takes two (3,1)-tetrahedra as input.

    static bool move23u(Tetra::Label t31, Tetra::Label t22);  // (2,3)-move upward (shift, Sec. 2.3.3)
    static bool move32u(Tetra::Label t31, Tetra::Label t22l, Tetra::Label t22r);  // (3,2)-move upward (ishift)
    static bool move23d(Tetra::Label t13, Tetra::Label t22);  // (2,3)-move downward
    static bool move32d(Tetra::Label t13, Tetra::Label t22l, Tetra::Label t22r);  // (3,2)-move downward
    // Comment: Shift/ishift moves adjust interpolation between slices (Fig. 5). 'u'/'d' denote direction.
    // HPC Target [GPU #7]: Parallelize move attempts across tetrahedra.

    static void updateVertexData();    // Updates Vertex data (e.g., cnum, Sec. 3.2.2)
    static void updateHalfEdgeData();  // Updates HalfEdge data
    static void updateTriangleData();  // Updates Triangle data
    // Comment: Maintains geometric consistency post-move (Sec. 3.2).

    static void updateGeometry();  // Full geometry update
    // Comment: Likely calls above updates; prepares for measurements (Sec. 3.2).

    static std::vector<Vertex::Label> vertices;            // All vertex labels
    static std::vector<Tetra::Label> tetras;              // All tetra labels
    static std::vector<HalfEdge::Label> halfEdges;        // All half-edge labels
    static std::vector<Triangle::Label> triangles;        // All triangle labels
    // Comment: Store full simplex lists for reconstruction (Sec. 3.2).
    // HPC Target [General #10]: Pre-allocate to avoid resizing.

    static std::vector<std::vector<Vertex::Label>> vertexNeighbors;       // Vertex adjacency lists
    static std::vector<std::array<Triangle::Label, 3>> triangleNeighbors; // Triangle adjacency (3 neighbors each)
    // Comment: Reconstructed connectivity for measurements (Sec. 3.2.2).
    // HPC Target [OpenMP #3, GPU #8]: Used in BFS; parallelize construction/use.

    static void check();  // Validates triangulation (e.g., manifold properties, Sec. 1.3)
    // Comment: Debugging/diagnostic tool.

private:
    Universe() {}  // Private constructor; enforces static-only usage
    // Comment: Singleton-like design; all members static (Sec. 3).

    static std::default_random_engine rng;  // Random number generator
    // Comment: Used for move selection (Sec. 2.1).
    // HPC Target [MPI #4]: Needs per-rank instances for distributed runs.
};

// HPC Targets Summary:
// [General #10]: Pre-allocate vectors/Bags for cache efficiency.
// [GPU #7]: Accelerate move functions (e.g., move26) on GPU.
// [OpenMP #3, GPU #8]: Parallelize BFS using vertexNeighbors/triangleNeighbors.
// [MPI #4]: Distribute RNG and simulation instances.