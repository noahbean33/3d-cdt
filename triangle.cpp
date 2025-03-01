// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "triangle.hpp"  // Includes Triangle class definition

// Comment: This is the implementation file for the Triangle class, representing spatial 2-simplices in 3D CDT (Sec. 3.2.2).
// Note: Contains no additional code beyond the header inclusion; all Triangle methods are defined inline in triangle.hpp.

// Purpose: Serves as a compilation unit for the Triangle class, though all functionality (setVertices, setHalfEdges, 
// setTriangleNeighbors, hasVertex) is implemented in the header (Sec. 3.2). This file ensures the class is compiled 
// and linked properly with other components (e.g., universe.cpp).

// HPC Targets Summary:
// [General #10]: No dynamic allocation here; cache efficiency relies on triangle.hpp’s std::array usage.
// [OpenMP #3, GPU #8]: No direct parallelizable methods here, but hasVertex() usage in BFS (e.g., universe.cpp’s 
// updateTriangleData or observables) could benefit from parallelization (Sec. 3.4).