// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "vertex.hpp"    // Vertex class for vertex references
#include "halfedge.hpp"  // HalfEdge class definition
#include "tetra.hpp"     // Tetra class for tetrahedron references

// Comment: This is the implementation file for the HalfEdge class, representing directed edges in the spatial triangulation of 3D CDT (Sec. 3.2.2).
// Note: Contains no additional code beyond header inclusions; all HalfEdge methods (setVertices, getAdjacent, setAdjacent) are defined inline in halfedge.hpp.

// Purpose: Serves as a compilation unit for the HalfEdge class, though all functionality is implemented inline in the header (Sec. 3.2). This file ensures the class is compiled and linked properly with other components (e.g., universe.cpp, observable.cpp).

// Role in CDT: The HalfEdge class manages half-edges, which are directed edges of spatial triangles in (3,1)-tetrahedra (Sec. 3.2.2). It supports connectivity between vertices, tetrahedra, and triangles, facilitating geometry reconstruction (e.g., in Universe::updateHalfEdgeData) and traversal (e.g., Observable::sphere2dDual).

// HPC Targets Summary:
// [General #10]: No dynamic allocation here; cache efficiency relies on halfedge.hpp’s std::array usage for vs and Pool’s fixed-size elements array (Sec. 3.1.1).
// [OpenMP #3, GPU #8]: No direct parallelizable methods here, but HalfEdge usage in BFS (e.g., Observable::sphere2dDual via Triangle::trnbr and HalfEdge::adj) could benefit from parallelization (Sec. 3.4).