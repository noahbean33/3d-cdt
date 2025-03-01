// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "tetra.hpp"  // Includes Tetra class definition

// Comment: This is the implementation file for the Tetra class, representing 3-simplices (tetrahedra) in 3D CDT (Sec. 2.3, 3.2.2).
// Note: Contains no additional code beyond the header inclusion; all Tetra methods are defined inline in tetra.hpp.

// Purpose: Serves as a compilation unit for the Tetra class, though all functionality (setVertices, setHalfEdges, 
// setTetras, is31, hasVertex, neighborsTetra, getTetraOpposite, etc.) is implemented inline in the header (Sec. 3.2.2). 
// This file ensures the class is compiled and linked properly with other components (e.g., universe.cpp).

// Role in CDT: The Tetra class is central to the triangulation, managing tetrahedron types ((3,1), (1,3), (2,2)), 
// vertices, neighbors, and half-edges (for (3,1)-tetras). It supports Monte Carlo moves (Sec. 2.3) and geometry 
// reconstruction (Sec. 3.2.2), with instances created and updated in universe.cpp (e.g., initialize, move26).

// HPC Targets Summary:
// [General #10]: No dynamic allocation here; cache efficiency relies on tetra.hpp’s std::array usage for vs, tnbr, hes.
// [OpenMP #3, GPU #8]: No direct parallelizable methods here, but hasVertex(), neighborsTetra() usage in BFS 
// (e.g., universe.cpp’s updateVertexData or observables) could benefit from parallelization (Sec. 3.4).