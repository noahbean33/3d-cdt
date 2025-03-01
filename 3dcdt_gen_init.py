#!/usr/bin/env python
# coding: utf-8
# parameters: g T outfile

import sys
import gc
import numpy as np
import pylab as pl  # Unused; likely for debugging
import math as m    # Unused; possibly vestigial
import time         # Unused; possibly for profiling
import random       # Unused; possibly vestigial
import itertools    # Unused; possibly vestigial
import re           # Unused; possibly vestigial
import fileinput    # Unused; possibly vestigial

# Comment: Script generates initial 3D CDT triangulation with genus g and T time slices, output to outfile (Sec. 3.1).

# Equations defining the triangulation topology (Sec. 1.3, 3.2):
# I.      X*t - N1TL + N2TL - N22 = 0
# II.     N2SL + N2TL = 2(N31 + N22)
# III.    N31 = 4/3 N1SL = 2N2SL
# IV.     N1SL = 3/2 N2SL
# V.      N0 = X*t + 1/2 N2SL
# Derived: 2N3 + N22 + N1TL - X*t = 0
#          N31 = -3/2 N22 - 1/2 N1TL + X*t = 0
# Comment: X is Euler characteristic; N0, N1SL, N1TL, N2SL, N2TL, N22, N31 are simplex counts.

def prepare_vertex_array_sphere(T):  # Vertex array for genus 0 (sphere, S^2)
    start = 0
    V = []
    for t in range(T):  # T slices
        v = np.arange(start, start + 5, 1)  # 5 vertices per slice
        start = max(v) + 1
        V.append(v)
    vertex = np.zeros(T * 5)  # Total vertices
    for t in range(T):
        for item in V[t]:
            vertex[item] = t  # Assign time slice
    return vertex  # Comment: Returns vertex times for S^1 x S^2 topology (Sec. 2.3).

def prepare_vertex_array(g, T):  # Vertex array for genus g > 0
    V = []
    start = 0
    v = GetInnerVertices(g, T, start)  # Inner vertices
    V.append(v)
    for s in range(1, g):  # Genus segments
        start = max(v.flatten()) + 1
        v = GetInnerVertices(g, T, start)
        V.append(v)
    Missing = GetMissingVertices(g, T)  # Additional vertices
    missing_max = max(np.concatenate(Missing).flatten()) + 1
    Corner = GetCornerVertices(g, T, missing_max)  # Corner vertices
    vertex = np.zeros(1 + max(Corner.flatten()))  # Total vertex array
    for structure in V:  # Inner vertices
        for t in range(T):
            for item in structure[t]:
                vertex[item] = t
    if g > 1:  # Multi-genus case
        for structure in Missing:
            for t in range(T):
                for item in structure[t]:
                    vertex[item] = t
        for structure in Corner:
            for t in range(T):
                for item in structure[t]:
                    vertex[item] = t
    if g == 1:  # Torus case
        for t in range(T):
            for item in Missing[t]:
                vertex[item] = t
            for item in Corner[t]:
                vertex[item] = t
    return np.asarray(vertex, dtype='int')  # Comment: Assigns time slices (Sec. 3.1).

def FindTrianglePairs(list_):  # Finds tetrahedron pairs sharing triangles
    pairs = []
    for i in range(len(list_)):
        p0, p1, p2, p3 = list_[i]
        lista0 = np.argwhere(list_.T == p0).T[1]
        lista1 = np.argwhere(list_.T == p1).T[1]
        lista2 = np.argwhere(list_.T == p2).T[1]
        ov01 = list(set(lista0).intersection(lista1))
        ov12 = list(set(lista1).intersection(lista2))
        triangle = list(set(ov01).intersection(ov12))
        pairs.append(triangle)
    return np.asarray(pairs)  # Comment: Identifies neighbor tetrahedra (Sec. 3.2).

def getTrianglelList(list3):  # Extracts unique triangles from tetrahedra
    list2 = []
    for simp in list3:
        t1 = list(np.sort([simp[0], simp[1], simp[2]]))
        t2 = list(np.sort([simp[0], simp[1], simp[3]]))
        t3 = list(np.sort([simp[0], simp[2], simp[3]]))
        t4 = list(np.sort([simp[1], simp[2], simp[3]]))
        for t in [t1, t2, t3, t4]:
            if t not in list2:
                list2.append(t)
    return np.asarray(list2)  # Comment: Builds triangle list (Sec. 3.2).

def getLinklList(list2):  // Extracts unique links (edges) from triangles
    list1 = []
    for tria in list2:
        for el in [[tria[0], tria[1]], [tria[0], tria[2]], [tria[1], tria[2]]]:
            if el not in list1:
                list1.append(el)
    return np.asarray(list1)  # Comment: Builds edge list (Sec. 3.2).

def PrepDat(vertexlist, simplexlist, list_):  # Prepares output data
    dat = [len(vertexlist)]  # Number of vertices
    dat.extend(vertexlist)   # Vertex times
    dat.append(len(vertexlist))  # Vertex count check
    dat.append(len(simplexlist))  # Number of tetrahedra
    for i in range(len(simplexlist)):
        dat.extend(list_[i])     # Tetrahedron vertices
        dat.extend(simplexlist[i])  # Neighbor indices
    dat.append(len(simplexlist))  # Tetra count check
    return np.asarray(dat)  # Comment: Formats for CDT file (Sec. 3.1).

# Block generation functions for tetrahedra (Sec. 3.2)
def GetBlocks_up_left(a, b, c, d, e, f):
    return np.asarray([a, b, c, e]), np.asarray([a, c, d, e]), np.asarray([a, d, e, f])

def GetBlocks_up_right(a, b, c, d, e, f):
    return np.asarray([a, b, c, d]), np.asarray([b, c, d, f]), np.asarray([b, d, e, f])

def GetBlocks_down_left(a, b, c, d, e, f):
    return np.asarray([a, b, c, e]), np.asarray([a, c, e, f]), np.asarray([a, d, e, f])

def GetBlocks_down_right(a, b, c, d, e, f):
    return np.asarray([a, b, c, f]), np.asarray([a, b, e, f]), np.asarray([a, d, e, f])

def GetTR(list3):  # Validates tetrahedron triangulation
    TR, bad_tetra = [], []
    for i, tetra in enumerate(list3):
        triangles = [
            np.asarray([tetra[0], tetra[1], tetra[2]]),
            np.asarray([tetra[0], tetra[1], tetra[3]]),
            np.asarray([tetra[0], tetra[2], tetra[3]]),
            np.asarray([tetra[1], tetra[2], tetra[3]])
        ]
        for tri in triangles:
            p1 = np.argwhere(list3.T == tri[0]).T[1]
            p2 = np.argwhere(list3.T == tri[1]).T[1]
            p3 = np.argwhere(list3.T == tri[2]).T[1]
            int12 = list(set(p1).intersection(p2))
            int13 = list(set(p1).intersection(p3))
            TR.append(len(list(set(int12).intersection(int13))))
            if len(list(set(int12).intersection(int13))) != 2:
                bad_tetra.append(np.asarray([i, tri, len(list(set(int12).intersection(int13)))]))
    not2 = [i for i, tr in enumerate(TR) if tr != 2]
    return np.asarray(not2), np.asarray(bad_tetra)  # Comment: Checks each face has 2 tetrahedra.

def GetInnerVertices(g, T, start):  # Inner vertices for genus g
    V = []
    prev_max = start
    for t in range(T):
        v_cur = np.arange(prev_max, prev_max + 12, 1)  # 12 vertices per slice
        prev_max = max(v_cur) + 1
        V.append(v_cur)
    V.append(np.arange(start, start + 12, 1))  # Wrap around
    return np.asarray(V)  # Comment: Core structure vertices (Sec. 3.1).

def GenCenter(g, T):  # Generates central tetrahedra
    start = 0
    Simps = [GenInnerStructure(GetInnerVertices(g, T, start))]
    for s in range(1, g):
        start = max(Simps[-1].flatten()) + 1
        Simps.append(GenInnerStructure(GetInnerVertices(g, T, start)))
    Simps = np.asarray(Simps).reshape(-1, 4)  # Flatten to tetra list
    return Simps  # Comment: Central structure for genus g (Sec. 3.1).

def GenInnerStructure(v):  # Inner tetrahedra for one segment
    Simps = []
    for t in range(T):
        vv, vv_2 = v[t], v[t + 1]
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 = vv
        a0_2, a1_2, a2_2, a3_2, a4_2, a5_2, a6_2, a7_2, a8_2, a9_2, a10_2, a11_2 = vv_2
        for s in [
            GetBlocks_up_right(a4, a0, a1, a4_2, a0_2, a1_2),
            GetBlocks_up_right(a5, a1, a2, a5_2, a1_2, a2_2),
            # ... (additional calls omitted for brevity)
            GetBlocks_down_right(a3, a10, a11, a3_2, a10_2, a11_2)
        ]:
            Simps.extend(s)
    return np.asarray(Simps)  # Comment: Generates tetrahedra blocks (Sec. 3.1).

# Additional functions (GetInnerOuters, GetVsingle, GetV0, GetVmid, GetVlast, GetCornerVertices, GetMissingVertices, GenMissing, GenSphere, FindPairs) omitted for brevity but follow similar structure: generate vertices or tetrahedra for specific topology parts.

g = int(sys.argv[1])  # Genus (0 for sphere, >0 for higher genus)
T = int(sys.argv[2])  # Number of time slices

if g == 0:
    list3 = GenSphere(T)  # Sphere topology (S^1 x S^2)
else:
    S_center = GenCenter(g, T)  # Central tetrahedra
    S_missing = GenMissing(g, T)  # Boundary tetrahedra
    list3 = np.vstack((S_center, S_missing))  # Combine

Simplex_list = FindPairs(list3)  # Neighboring tetrahedra pairs
num0 = max(list3.flatten()) + 1  # Total vertex count

if g > 0:
    vertex = prepare_vertex_array(g, T)  # Higher genus vertex times
else:
    vertex = prepare_vertex_array_sphere(T)  # Sphere vertex times

dat = PrepDat(vertex, Simplex_list, list3)  # Format output data

outfile = sys.argv[3]  # Output file path

with open(outfile, 'w') as file_handler:
    file_handler.write("0\n")  # Ordered flag (0 = unordered)
    for item in dat:
        file_handler.write("{}\n".format(int(item)))  # Write integers
# Comment: Generates triangulation file for Universe::initialize (Sec. 3.1).