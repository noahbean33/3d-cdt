// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include <unordered_map>    // For std::unordered_map (sliceEdges, done)
#include <algorithm>        // For std::sort, std::find
#include "minbu.hpp"        // Minbu class definition

// Comment: Implements Minbu to measure minimal bunching (Minbu) observable (Sec. 3.4).
void Minbu::process() {
    std::string tmp = "";  // Output string

    int slice;  // Target time slice
    for (int i = 0; i < Universe::nSlices; i++) {
        if (Universe::sliceSizes[i] == Simulation::target2Volume) {  // Find slice matching target2Volume
            slice = i;
            break;
        }
    }
    // Comment: Selects a slice with vertex count matching target2Volume (Sec. 2.4).

    std::unordered_map<int, HalfEdge::Label> sliceEdges;  // Half-edges in target slice
    for (auto he : Universe::halfEdges) {
        if (he->vs[0]->time == slice) sliceEdges[he] = he;  // Collect edges where start vertex is in slice
    }
    // Comment: Maps half-edge indices to labels for slice-specific edges (Sec. 3.2).

    std::unordered_map<int, bool> done;  // Tracks processed half-edges
    std::vector<std::array<Vertex::Label, 3>> minNecks;  // Minimal necks (triangles)

    for (auto& it : sliceEdges) {  // Iterate over slice edges
        auto he = it.second;
        if (done.find(he) != done.end()) continue;  // Skip processed edges

        auto start = he->next->adj->next;  // Start of front cycle
        std::vector<HalfEdge::Label> fronts;  // Front half-edges
        HalfEdge::Label cur = start;
        do {  // Traverse front cycle
            fronts.push_back(cur);
            cur = cur->adj->next;
        } while (cur->adj->next->adj != he);

        start = he->prev->adj->prev;  // Start of back cycle
        std::vector<HalfEdge::Label> backs;  // Back half-edges
        cur = start;
        do {  // Traverse back cycle
            backs.push_back(cur);
            cur = cur->adj->prev;
        } while (cur->adj->prev->adj != he);

        for (auto f : fronts) {  // Find minimal necks
            for (auto b : backs) {
                if (f->vs[1] == b->vs[0] && done.find(b) == done.end() && done.find(f) == done.end()) {
                    std::array<Vertex::Label, 3> nek = {f->vs[0], b->vs[1], f->vs[1]};  // Triangle vertices
                    std::sort(nek.begin(), nek.end());  // Sort for uniqueness
                    minNecks.push_back(nek);  // Add minimal neck
                    for (auto t : Universe::tetras31) {  // Verify no tetra contains this triangle
                        if (t->hasVertex(nek[0]) && t->hasVertex(nek[1]) && t->hasVertex(nek[2])) {
                            t->log();
                            std::cout << std::flush;
                            assert(false);  // Should not happen in minimal neck
                        }
                    }
                }
            }
        }
        done[he] = true;       // Mark edge processed
        done[he->adj] = true;  // Mark adjacent edge processed
    }
    // Comment: Identifies minimal necks (triangles separating regions) in the slice’s spatial graph.

    std::vector<std::array<Vertex::Label, 2>> neckLinks;  // Edges of minimal necks
    std::sort(minNecks.begin(), minNecks.end());  // Ensure uniqueness
    for (auto n : minNecks) {
        neckLinks.push_back({n[0], n[1]});
        neckLinks.push_back({n[1], n[2]});
        neckLinks.push_back({n[0], n[2]});
    }
    // Comment: Extracts edges from minimal neck triangles.

    printf("------\n");

    std::vector<int> histogram;  // Histogram of region sizes
    histogram.resize(Simulation::target2Volume / 2 + 1, 0);  // Size bins up to half target2Volume
    for (auto neck : minNecks) {
        Triangle::Label tr;  // Find a triangle containing neck vertex
        for (auto tri : Universe::triangles) {
            if (tri->hasVertex(neck[0]) || tri->hasVertex(neck[1]) || tri->hasVertex(neck[2])) {
                tr = tri;
                break;
            }
        }

        auto origin = tr;  // Start BFS from this triangle
        std::vector<Triangle::Label> tdone;  // Visited triangles
        std::vector<Triangle::Label> thisDepth;  // Current BFS depth
        std::vector<Triangle::Label> nextDepth;  // Next BFS depth

        tdone.push_back(origin);
        thisDepth.push_back(origin);

        int totalTr = 0;  // Count of triangles in region
        int currentDepth = 0;
        do {  // BFS to count triangles in region
            for (auto t : thisDepth) {
                for (auto he : t->hes) {
                    auto v1 = he->vs[0];
                    auto v2 = he->vs[1];
                    if ((v1 == neck[0] || v1 == neck[1] || v1 == neck[2]) && 
                        (v2 == neck[0] || v2 == neck[1] || v2 == neck[2])) {
                        continue;  // Skip neck edges to separate regions
                    }
                    auto neighbor = he->adj->triangle;  // Adjacent triangle via half-edge
                    if (std::find(tdone.begin(), tdone.end(), neighbor) == tdone.end()) {
                        nextDepth.push_back(neighbor);
                        tdone.push_back(neighbor);
                        totalTr++;  // Increment region size
                    }
                }
            }
            thisDepth = nextDepth;
            nextDepth.clear();
            currentDepth++;
        } while (thisDepth.size() > 0);

        if (totalTr + 1 < Simulation::target2Volume / 2) {  // Bin region size
            histogram.at(totalTr + 1) += 1;  // Smaller region
        } else {
            histogram.at(Simulation::target2Volume - totalTr - 1) += 1;  // Larger region’s complement
        }
    }
    // Comment: BFS counts triangles in regions separated by minimal necks; builds histogram (Sec. 3.4).

    for (auto h : histogram) {  // Format output
        tmp += std::to_string(h) + " ";
    }
    output = tmp;  // Set for Observable::write(); trailing space not removed
    // HPC Targets [OpenMP #2, GPU #8]: Parallelize neck loop; optimize BFS and histogram.
}