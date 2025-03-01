// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include <fstream>        // For file I/O (initialize, exportGeometry)
#include <algorithm>      // For std::find, std::fill
#include <unordered_map>  // For vertex/tetra mapping in exportGeometry
#include "universe.hpp"   // Universe class definition

// Static member initializations (Sec. 3)
int Universe::nSlices = 0;  // Number of time slices (S^1, Sec. 2.3)
std::vector<int> Universe::slabSizes;  // Tetrahedra between slices
std::vector<int> Universe::sliceSizes; // Vertices per spatial slice (S^2)
std::string Universe::fID;             // File identifier
std::string Universe::OutFile;         // Output file name
int Universe::strictness;              // Manifold strictness (Sec. 1.3)
int Universe::volfix_switch;           // Volume fixing toggle (Sec. 2.4)
std::default_random_engine Universe::rng(0);  // RNG for moves (Sec. 2.1)
// HPC Target [MPI #4]: Static RNG needs per-rank instances for distributed runs.
Bag<Tetra, Tetra::pool_size> Universe::tetrasAll(rng);  // All tetrahedra (Sec. 3.1.2)
Bag<Tetra, Tetra::pool_size> Universe::tetras31(rng);   // (3,1)-tetrahedra
Bag<Vertex, Vertex::pool_size> Universe::verticesAll(rng);  // All vertices
Bag<Vertex, Vertex::pool_size> Universe::verticesSix(rng);  // Vertices with cnum=6 (Sec. 2.3.1)

std::vector<Vertex::Label> Universe::vertices;         // Vertex labels
std::vector<Tetra::Label> Universe::tetras;           // Tetra labels
std::vector<HalfEdge::Label> Universe::halfEdges;     // Half-edge labels
std::vector<Triangle::Label> Universe::triangles;     // Triangle labels
std::vector<std::vector<Vertex::Label>> Universe::vertexNeighbors;  // Vertex adjacency
std::vector<std::array<Triangle::Label, 3>> Universe::triangleNeighbors;  // Triangle adjacency
// HPC Target [General #10]: Pre-allocate vectors for cache efficiency.

bool Universe::initialize(std::string geometryFilename, std::string fID_, int strictness_, int volfix_switch_) {
    fID = fID_;  // Set file ID
    std::ifstream infile(geometryFilename.c_str());  // Open input geometry file
    assert(!infile.fail());  // Ensure file opened successfully

    strictness = strictness_;  // Set manifold strictness
    volfix_switch = volfix_switch_;  // Set volume fixing switch

    bool ordered;  // Flag for tetrahedron ordering convention
    infile >> ordered;  // Read ordering flag

    int n0;  // Number of vertices
    infile >> n0;
    printf("n0: %d\n", n0);  // Log vertex count
    int line;  // Temp variable for file validation

    int maxTime = 0;  // Max time slice
    std::vector<Vertex::Label> vs(n0);  // Temporary vertex storage
    for (int i = 0; i < n0; i++) {
        infile >> line;  // Read vertex time
        auto v = Vertex::create();  // Create vertex (Sec. 3.1.1)
        verticesAll.add(v);  // Add to Bag
        v->time = line;  // Set time slice
        vs.at(i) = v;
        if (v->time > maxTime) maxTime = v->time;  // Track max time
    }
    infile >> line;
    if (line != n0) return false;  // Validate vertex count

    nSlices = maxTime + 1;  // Set number of slices
    slabSizes.resize(maxTime + 1);  // Resize slab sizes
    sliceSizes.resize(maxTime + 1); // Resize slice sizes
    std::fill(slabSizes.begin(), slabSizes.end(), 0);  // Initialize to 0
    std::fill(sliceSizes.begin(), sliceSizes.end(), 0);

    int n3;  // Number of tetrahedra
    infile >> n3;
    printf("n3: %d\n", n3);  // Log tetra count
    for (int i = 0; i < n3; i++) {
        auto t = Tetra::create();  // Create tetrahedron (Sec. 3.1.1)
        int tvs[4];  // Vertex labels
        for (int j = 0; j < 4; j++) infile >> tvs[j];  // Read vertices
        int tts[4];  // Neighbor tetra labels
        for (int j = 0; j < 4; j++) infile >> tts[j];  // Read neighbors

        t->setVertices(tvs[0], tvs[1], tvs[2], tvs[3]);  // Set tetra vertices
        if (t->is31()) {  // If (3,1)-tetra (Sec. 2.3)
            for (int j = 0; j < 3; j++) Pool<Vertex>::at(tvs[j])->tetra = t;  // Link to base vertices
        }
        t->setTetras(tts[0], tts[1], tts[2], tts[3]);  // Set tetra neighbors

        tetrasAll.add(t);  // Add to all tetra Bag
        if (t->is31()) tetras31.add(t);  // Add to (3,1) Bag
        slabSizes.at(t->vs[1]->time) += 1;  // Update slab size
        if (t->is31()) sliceSizes.at(t->vs[0]->time) += 1;  // Update slice size
    }
    infile >> line;
    if (line != n3) return false;  // Validate tetra count
    printf("read %s\n", geometryFilename.c_str());

    if (!ordered) {  // Reorder tetrahedra if not in convention (Sec. 3.2.2)
        for (auto t : tetrasAll) {
            auto tnbr = t->tnbr;
            Tetra::Label t012 = -1, t013 = -1, t023 = -1, t123 = -1;
            for (auto tn : tnbr) {  // Assign neighbors by opposite vertex
                if (!tn->hasVertex(t->vs[3])) { t012 = tn; continue; }
                if (!tn->hasVertex(t->vs[2])) { t013 = tn; continue; }
                if (!tn->hasVertex(t->vs[1])) { t023 = tn; continue; }
                if (!tn->hasVertex(t->vs[0])) { t123 = tn; continue; }
            }
            assert(t012 >= 0 && t013 >= 0 && t023 >= 0 && t123 >= 0);
            t->setTetras(t123, t023, t013, t012);  // Reorder to standard convention
        }
    }

    for (auto v : vs) {  // Compute initial coordination numbers (Sec. 3.2.2)
        int cnum = 0, scnum = 0;
        for (auto t : tetrasAll) {
            if (t->hasVertex(v)) cnum++;  // Total tetrahedra count
            if (!t->is31()) continue;
            if (t->hasVertex(v) && t->vs[3] != v) scnum++;  // Spatial neighbors in (3,1)
        }
        v->scnum = scnum;
        v->cnum = cnum;
    }
    return true;  // Initialization successful
    // HPC Target [General #10]: Pre-allocate vs, slabSizes for cache efficiency.
}

bool Universe::exportGeometry(std::string geometryFilename) {
    updateGeometry();  // Ensure latest state (Sec. 3.2)

    std::unordered_map<int, int> vertexMap;  // Map Vertex::Label to file index
    std::vector<Vertex::Label> intVMap(vertices.size());  // Reverse mapping
    int i = 0;
    for (auto v : vertices) {
        vertexMap[v] = i;
        intVMap[i] = v;
        i++;
    }

    std::unordered_map<int, int> tetraMap;  // Map Tetra::Label to file index
    std::vector<Tetra::Label> intTMap(tetrasAll.size());  // Reverse mapping
    i = 0;
    for (auto t : tetrasAll) {
        tetraMap[t] = i;
        intTMap[i] = t;
        i++;
    }

    std::string out = "1\n";  // Indicate ordered tetrahedra
    out += std::to_string(vertices.size()) + "\n";  // Vertex count
    for (auto v : intVMap) out += std::to_string(v->time) + "\n";  // Vertex times
    out += std::to_string(vertices.size()) + "\n";  // Vertex count check
    out += std::to_string(tetrasAll.size()) + "\n";  // Tetra count

    for (auto t : intTMap) {  // Write tetrahedra data
        for (auto v : t->vs) out += std::to_string(vertexMap[v]) + "\n";
        for (auto tn : t->tnbr) out += std::to_string(tetraMap[tn]) + "\n";
    }
    out += std::to_string(tetrasAll.size());  // Tetra count check

    std::ofstream file(geometryFilename, std::ios::out | std::ios::trunc);
    assert(file.is_open());
    file << out << "\n";
    file.close();
    std::cout << geometryFilename << "\n";
    return true;  // Export successful
    // HPC Note: String concatenation inefficient; minimal impact unless frequent.
}

bool Universe::move26(Tetra::Label t) {  // (2,6)-move (Sec. 2.3.1, Fig. 3)
    assert(t->is31());  // Input must be (3,1)-tetra
    int time = t->vs[0]->time;  // Base time slice
    auto tv = t->tnbr[3];  // Vertical (1,3) neighbor
    assert(tv->is13());

    auto vn = Vertex::create();  // New vertex at center
    verticesAll.add(vn);
    vn->time = time;
    vn->scnum = 3;  // 3 spatial neighbors
    vn->cnum = 6;   // 6 tetrahedra

    auto v0 = t->vs[0], v1 = t->vs[1], v2 = t->vs[2], vt = t->vs[3], vb = tv->vs[0];  // Existing vertices
    auto tn01 = Tetra::create(), tn12 = Tetra::create(), tn20 = Tetra::create();  // New (3,1)-tetras
    auto tvn01 = Tetra::create(), tvn12 = Tetra::create(), tvn20 = Tetra::create();  // New (1,3)-tetras

    tetrasAll.add(tn01); tetrasAll.add(tn12); tetrasAll.add(tn20);
    tetras31.add(tn01); tetras31.add(tn12); tetras31.add(tn20);
    tetrasAll.add(tvn01); tetrasAll.add(tvn12); tetrasAll.add(tvn20);

    auto to0 = t->getTetraOpposite(v0), to1 = t->getTetraOpposite(v1), to2 = t->getTetraOpposite(v2);  // Old neighbors
    auto tvo0 = tv->getTetraOpposite(v0), tvo1 = tv->getTetraOpposite(v1), tvo2 = tv->getTetraOpposite(v2);

    tn01->setVertices(v0, v1, vn, vt);  // Set new tetra vertices
    tn12->setVertices(v1, v2, vn, vt);
    tn20->setVertices(v2, v0, vn, vt);
    tvn01->setVertices(vb, v0, v1, vn);
    tvn12->setVertices(vb, v1, v2, vn);
    tvn20->setVertices(vb, v2, v0, vn);

    tn01->setTetras(tn12, tn20, to2, tvn01);  // Set new tetra neighbors
    tn12->setTetras(tn20, tn01, to0, tvn12);
    tn20->setTetras(tn01, tn12, to1, tvn20);
    tvn01->setTetras(tn01, tvn12, tvn20, tvo2);
    tvn12->setTetras(tn12, tvn20, tvn01, tvo0);
    tvn20->setTetras(tn20, tvn01, tvn12, tvo1);

    to0->exchangeTetraOpposite(t->getVertexOpposite(v0), tn12);  // Update external neighbors
    to1->exchangeTetraOpposite(t->getVertexOpposite(v1), tn20);
    to2->exchangeTetraOpposite(t->getVertexOpposite(v2), tn01);
    tvo0->exchangeTetraOpposite(tv->getVertexOpposite(v0), tvn12);
    tvo1->exchangeTetraOpposite(tv->getVertexOpposite(v1), tvn20);
    tvo2->exchangeTetraOpposite(tv->getVertexOpposite(v2), tvn01);

    slabSizes[time] += 2;  // Update geometry stats
    slabSizes[(time - 1 + nSlices) % nSlices] += 2;
    sliceSizes[time] += 2;

    tetrasAll.remove(t); tetras31.remove(t); tetrasAll.remove(tv);  // Remove old tetras
    Tetra::destroy(t); Tetra::destroy(tv);

    vn->tetra = tn01;  // Update vertex tetra links
    v0->tetra = tn01; v1->tetra = tn12; v2->tetra = tn20;

    v0->scnum++; v1->scnum++; v2->scnum++;  // Update coordination numbers
    v0->cnum += 2; v1->cnum += 2; v2->cnum += 2; vt->cnum += 2; vb->cnum += 2;

    return true;  // Move successful
    // HPC Target [GPU #7]: Batch multiple move26 calls on GPU.
}

bool Universe::move62(Vertex::Label v) {  // (6,2)-move (Sec. 2.3.1, Fig. 3)
    assert(v->cnum == 6);  // Input must have 6 tetrahedra
    int time = v->time;
    auto t01 = v->tetra;  // Starting (3,1)-tetra
    auto tv01 = t01->tnbr[3];  // Vertical (1,3)-tetra

    int vpos = -1;
    for (int i = 0; i < 3; i++) if (t01->vs[i] == v) { vpos = i; break; }  // Find v’s position
    assert(vpos >= 0);

    auto v0 = t01->vs[(vpos + 1) % 3], v1 = t01->vs[(vpos + 2) % 3], v2 = t01->getVertexOpposite(v0);  // Neighbor vertices
    auto t12 = t01->getTetraOpposite(v0), t20 = t01->getTetraOpposite(v1);  // Other (3,1)-tetras
    auto tv12 = tv01->getTetraOpposite(v0), tv20 = tv01->getTetraOpposite(v1);  // Other (1,3)-tetras

    assert(t01->is31() && t12->is31() && t20->is31() && tv01->is13() && tv12->is13() && tv20->is13());

    auto to01 = t01->getTetraOpposite(v), to12 = t12->getTetraOpposite(v), to20 = t20->getTetraOpposite(v);  // External neighbors
    auto tvo01 = tv01->getTetraOpposite(v), tvo12 = tv12->getTetraOpposite(v), tvo20 = tv20->getTetraOpposite(v);

    // Commented-out checks for tadpoles/self-energy (Sec. 1.3); replaced by strictness levels
    if (strictness == 1) {  // Disallow tadpoles
        if (v0->scnum < 3 || v1->scnum < 3 || v2->scnum < 3) return false;
    } else if (strictness >= 2) {  // Disallow self-energy
        if (v0->scnum < 4 || v1->scnum < 4 || v2->scnum < 4) return false;
    }

    auto tn = Tetra::create(), tvn = Tetra::create();  // New (3,1) and (1,3) tetras
    auto vt = t01->vs[3], vb = tv01->vs[0];  // Top/bottom vertices

    tetrasAll.add(tn); tetras31.add(tn); tetrasAll.add(tvn);

    tn->setVertices(v0, v1, v2, vt);  // Set new tetra vertices
    tvn->setVertices(vb, v0, v1, v2);

    tn->setTetras(to12, to20, to01, tvn);  // Set new tetra neighbors
    tvn->setTetras(tn, tvo12, tvo20, tvo01);

    v0->tetra = tn; v1->tetra = tn; v2->tetra = tn;  // Update vertex links

    v0->scnum--; v1->scnum--; v2->scnum--;  // Update coordination numbers
    v0->cnum -= 2; v1->cnum -= 2; v2->cnum -= 2; vt->cnum -= 2; vb->cnum -= 2;

    to01->exchangeTetraOpposite(t01->getVertexOpposite(v), tn);  // Update external neighbors
    to12->exchangeTetraOpposite(t12->getVertexOpposite(v), tn);
    to20->exchangeTetraOpposite(t20->getVertexOpposite(v), tn);
    tvo01->exchangeTetraOpposite(tv01->getVertexOpposite(v), tvn);
    tvo12->exchangeTetraOpposite(tv12->getVertexOpposite(v), tvn);
    tvo20->exchangeTetraOpposite(tv20->getVertexOpposite(v), tvn);

    tetrasAll.remove(t01); tetrasAll.remove(t12); tetrasAll.remove(t20);  // Remove old tetras
    tetras31.remove(t01); tetras31.remove(t12); tetras31.remove(t20);
    tetrasAll.remove(tv01); tetrasAll.remove(tv12); tetrasAll.remove(tv20);
    Tetra::destroy(t01); Tetra::destroy(t12); Tetra::destroy(t20);
    Tetra::destroy(tv01); Tetra::destroy(tv12); Tetra::destroy(tv20);

    verticesAll.remove(v); Vertex::destroy(v);  // Remove central vertex

    slabSizes[time] -= 2; slabSizes[(time - 1 + nSlices) % nSlices] -= 2; sliceSizes[time] -= 2;  // Update stats

    return true;
    // HPC Target [GPU #7]: Batch move62 calls on GPU.
}

bool Universe::move44(Tetra::Label t012, Tetra::Label t230) {  // (4,4)-move (Sec. 2.3.2, Fig. 4)
    Vertex::Label v0, v1, v2, v3;  // Identify vertices
    v1 = t012->getVertexOppositeTetra(t230);  // Opposite vertex in t012
    v3 = t230->getVertexOppositeTetra(t012);  // Opposite vertex in t230
    for (int i = 0; i < 3; i++) {  // Find shared vertices
        if (t012->vs[i] == v1) {
            v2 = t012->vs[(i + 1) % 3];
            v0 = t012->vs[(i + 2) % 3];
            break;
        }
    }

    auto tv012 = t012->tnbr[3], tv230 = t230->tnbr[3];  // Vertical (1,3) neighbors

    if (strictness >= 1 && v1 == v3) return false;  // Prevent degeneracy
    if (strictness >= 2 && (v0->scnum == 3 || v2->scnum == 3)) return false;  // Prevent self-energy
    if (strictness >= 3 && v1->neighborsVertex(v3)) return false;  // Prevent existing link (Sec. 2.3.2)

    auto vt = t012->vs[3], vb = tv012->vs[0];  // Top/bottom vertices
    auto ta01 = t012->getTetraOpposite(v2), ta12 = t012->getTetraOpposite(v0);  // External neighbors
    auto ta23 = t230->getTetraOpposite(v0), ta30 = t230->getTetraOpposite(v2);
    auto tva01 = tv012->getTetraOpposite(v2), tva12 = tv012->getTetraOpposite(v0);
    auto tva23 = tv230->getTetraOpposite(v0), tva30 = tv230->getTetraOpposite(v2);

    if (ta01 == t230 || ta23 == t012 || tva01 == tv230 || tva23 == tv012) return false;  // Prevent tadpoles

    auto tn013 = t230, tn123 = t012, tvn013 = tv230, tvn123 = tv012;  // Reuse tetras for flip
    tn013->setVertices(v0, v1, v3, vt); tn123->setVertices(v1, v2, v3, vt);  // New vertex sets
    tvn013->setVertices(vb, v0, v1, v3); tvn123->setVertices(vb, v1, v2, v3);

    tn013->setTetras(tn123, ta30, ta01, tvn013);  // New neighbor sets
    tn123->setTetras(ta23, tn013, ta12, tvn123);
    tvn013->setTetras(tn013, tvn123, tva30, tva01);
    tvn123->setTetras(tn123, tva23, tvn013, tva12);

    ta01->exchangeTetraOpposite(t012->getVertexOpposite(v2), tn013);  // Update external neighbors
    ta23->exchangeTetraOpposite(t230->getVertexOpposite(v0), tn123);
    tva01->exchangeTetraOpposite(tv012->getVertexOpposite(v2), tvn013);
    tva23->exchangeTetraOpposite(tv230->getVertexOpposite(v0), tvn123);

    v0->scnum--; v1->scnum++; v2->scnum--; v3->scnum++;  // Update spatial coordination
    v0->cnum -= 2; v1->cnum += 2; v2->cnum -= 2; v3->cnum += 2;  // Update total coordination
    v0->tetra = tn013; v2->tetra = tn123;  // Update tetra links

    if (strictness >= 2) assert(v0->scnum >= 3 && v2->scnum >= 3);  // Validate post-move

    return true;
    // HPC Target [GPU #7]: Batch move44 calls on GPU.
}

bool Universe::move23u(Tetra::Label t31, Tetra::Label t22) {  // (2,3)-move upward (Sec. 2.3.3, Fig. 5)
    Vertex::Label v0 = t31->getVertexOppositeTetra(t22), v1 = t22->getVertexOppositeTetra(t31);
    int v0pos = -1;
    for (int i = 0; i < 3; i++) if (t31->vs[i] == v0) { v0pos = i; break; }
    assert(v0pos >= 0);
    Vertex::Label v2 = t31->vs[(v0pos + 1) % 3], v4 = t31->vs[(v0pos + 2) % 3], v3 = t31->vs[3];

    Tetra::Label ta023 = t31->getTetraOpposite(v4), ta034 = t31->getTetraOpposite(v2);
    Tetra::Label ta123 = t22->getTetraOpposite(v4), ta124 = t22->getTetraOpposite(v3), ta134 = t22->getTetraOpposite(v2);

    if (ta023->hasVertex(v1) || ta123->hasVertex(v0) || ta034->hasVertex(v1) || ta134->hasVertex(v0) || v0->neighborsVertex(v1)) return false;  // Prevent invalid configs

    auto tn31 = Tetra::create(), tn22l = Tetra::create(), tn22r = Tetra::create();
    tetrasAll.add(tn31); tetras31.add(tn31); tetrasAll.add(tn22l); tetrasAll.add(tn22r);

    tn31->setVertices(v0, v2, v4, v1); tn22l->setVertices(v0, v2, v1, v3); tn22r->setVertices(v0, v4, v1, v3);
    tn31->setTetras(ta124, tn22r, tn22l, t31->tnbr[3]);
    tn22l->setTetras(ta123, tn22r, ta023, tn31);
    tn22r->setTetras(ta134, tn22l, ta034, tn31);

    int time = tn31->vs[0]->time;
    slabSizes[time] += 1;

    t31->tnbr[3]->exchangeTetraOpposite(t31->tnbr[3]->vs[0], tn31);
    ta023->exchangeTetraOpposite(t31->getVertexOpposite(v4), tn22l);
    ta034->exchangeTetraOpposite(t31->getVertexOpposite(v2), tn22r);
    ta123->exchangeTetraOpposite(t22->getVertexOpposite(v4), tn22l);
    ta124->exchangeTetraOpposite(t22->getVertexOpposite(v3), tn31);
    ta134->exchangeTetraOpposite(t22->getVertexOpposite(v2), tn22r);

    v0->cnum += 2; v1->cnum += 2;
    tetrasAll.remove(t31); tetras31.remove(t31); tetrasAll.remove(t22);
    Tetra::destroy(t31); Tetra::destroy(t22);

    tn31->vs[0]->tetra = tn31; tn31->vs[1]->tetra = tn31; tn31->vs[2]->tetra = tn31;

    return true;
    // HPC Target [GPU #7]: Batch move23u calls on GPU.
}

bool Universe::move32u(Tetra::Label t31, Tetra::Label t22l, Tetra::Label t22r) {  // (3,2)-move upward (Sec. 2.3.3)
    Vertex::Label v1 = t31->vs[3], v3 = t22l->getVertexOppositeTetra(t31), v4 = t31->getVertexOppositeTetra(t22l);
    int v4pos = -1;
    for (int i = 0; i < 3; i++) if (t31->vs[i] == v4) { v4pos = i; break; }
    assert(v4pos >= 0);
    Vertex::Label v0 = t31->vs[(v4pos + 1) % 3], v2 = t31->vs[(v4pos + 2) % 3];

    Tetra::Label ta023 = t22l->getTetraOpposite(v1), ta034 = t22r->getTetraOpposite(v1);
    Tetra::Label ta123 = t22l->getTetraOpposite(v0), ta124 = t31->getTetraOpposite(v0), ta134 = t22r->getTetraOpposite(v0);

    if (ta023->hasVertex(v4) || ta123->hasVertex(v4) || ta034->hasVertex(v2) || ta124->hasVertex(v3) || ta134->hasVertex(v2)) return false;

    auto tn31 = Tetra::create(), tn22 = Tetra::create();
    tetrasAll.add(tn31); tetras31.add(tn31); tetrasAll.add(tn22);

    tn31->setVertices(v0, v2, v4, v3); tn22->setVertices(v2, v4, v1, v3);
    tn31->setTetras(tn22, ta034, ta023, t31->tnbr[3]);
    tn22->setTetras(ta134, ta123, tn31, ta124);

    t31->tnbr[3]->exchangeTetraOpposite(t31->tnbr[3]->vs[0], tn31);
    ta023->exchangeTetraOpposite(t22l->getVertexOpposite(v1), tn31);
    ta034->exchangeTetraOpposite(t22r->getTetraOpposite(v1), tn31);
    ta123->exchangeTetraOpposite(t22l->getVertexOpposite(v0), tn22);
    ta124->exchangeTetraOpposite(t31->getVertexOpposite(v0), tn22);
    ta134->exchangeTetraOpposite(t22r->getTetraOpposite(v0), tn22);

    v0->cnum -= 2; v1->cnum -= 2;
    tetrasAll.remove(t31); tetras31.remove(t31); tetrasAll.remove(t22l); tetrasAll.remove(t22r);
    Tetra::destroy(t31); Tetra::destroy(t22l); Tetra::destroy(t22r);

    int time = tn31->vs[0]->time;
    slabSizes[time] -= 1;

    tn31->vs[0]->tetra = tn31; tn31->vs[1]->tetra = tn31; tn31->vs[2]->tetra = tn31;

    return true;
    // HPC Target [GPU #7]: Batch move32u calls on GPU.
}

bool Universe::move23d(Tetra::Label t13, Tetra::Label t22) {  // (2,3)-move downward
    Vertex::Label v0 = t13->getVertexOppositeTetra(t22), v1 = t22->getVertexOppositeTetra(t13);
    auto t31 = t13->tnbr[0];  // Access via (3,1) neighbor due to ordering
    int v0pos = -1;
    for (int i = 0; i < 3; i++) if (t31->vs[i] == v0) { v0pos = i; break; }
    assert(v0pos >= 0);
    Vertex::Label v2 = t31->vs[(v0pos + 1) % 3], v4 = t31->vs[(v0pos + 2) % 3], v3 = t13->vs[0];

    Tetra::Label ta023 = t13->getTetraOpposite(v4), ta034 = t13->getTetraOpposite(v2);
    Tetra::Label ta123 = t22->getTetraOpposite(v4), ta124 = t22->getTetraOpposite(v3), ta134 = t22->getTetraOpposite(v2);

    if (ta023->hasVertex(v1) || ta123->hasVertex(v0) || ta034->hasVertex(v1) || ta134->hasVertex(v0) || v0->neighborsVertex(v1)) return false;

    auto tn13 = Tetra::create(), tn22l = Tetra::create(), tn22r = Tetra::create();
    tetrasAll.add(tn13); tetrasAll.add(tn22l); tetrasAll.add(tn22r);

    tn13->setVertices(v1, v0, v2, v4); tn22l->setVertices(v1, v3, v0, v2); tn22r->setVertices(v1, v3, v0, v4);
    tn13->setTetras(t13->tnbr[0], ta124, tn22r, tn22l);
    tn22l->setTetras(ta023, tn13, ta123, tn22r);
    tn22r->setTetras(ta034, tn13, ta134, tn22l);

    int time = t31->vs[0]->time;
    slabSizes[time] += 1;

    t13->tnbr[0]->exchangeTetraOpposite(t13->tnbr[0]->vs[3], tn13);
    ta023->exchangeTetraOpposite(t13->getVertexOpposite(v4), tn22l);
    ta034->exchangeTetraOpposite(t13->getVertexOpposite(v2), tn22r);
    ta123->exchangeTetraOpposite(t22->getVertexOpposite(v4), tn22l);
    ta124->exchangeTetraOpposite(t22->getTetraOpposite(v3), tn13);
    ta134->exchangeTetraOpposite(t22->getTetraOpposite(v2), tn22r);

    v0->cnum += 2; v1->cnum += 2;
    tetrasAll.remove(t13); tetrasAll.remove(t22);
    Tetra::destroy(t13); Tetra::destroy(t22);

    return true;
    // HPC Target [GPU #7]: Batch move23d calls on GPU.
}

bool Universe::move32d(Tetra::Label t13, Tetra::Label t22l, Tetra::Label t22r) {  // (3,2)-move downward
    Vertex::Label v1 = t13->vs[0], v3 = t22l->getVertexOppositeTetra(t13), v4 = t13->getVertexOppositeTetra(t22l);
    auto t31 = t13->tnbr[0];
    int v4pos = -1;
    for (int i = 0; i < 3; i++) if (t31->vs[i] == v4) { v4pos = i; break; }
    assert(v4pos >= 0);
    Vertex::Label v0 = t31->vs[(v4pos + 1) % 3], v2 = t31->vs[(v4pos + 2) % 3];

    Tetra::Label ta023 = t22l->getTetraOpposite(v1), ta034 = t22r->getTetraOpposite(v1);
    Tetra::Label ta123 = t22l->getTetraOpposite(v0), ta124 = t13->getTetraOpposite(v0), ta134 = t22r->getTetraOpposite(v0);

    if (ta023->hasVertex(v4) || ta123->hasVertex(v4) || ta034->hasVertex(v2) || ta124->hasVertex(v3) || ta134->hasVertex(v2)) return false;

    auto tn13 = Tetra::create(), tn22 = Tetra::create();
    tetrasAll.add(tn13); tetrasAll.add(tn22);

    tn13->setVertices(v3, v0, v2, v4); tn22->setVertices(v1, v3, v2, v4);
    tn13->setTetras(t13->tnbr[0], tn22, ta034, ta023);
    tn22->setTetras(tn13, ta124, ta134, ta123);

    t13->tnbr[0]->exchangeTetraOpposite(t13->tnbr[0]->vs[3], tn13);
    ta023->exchangeTetraOpposite(t22l->getTetraOpposite(v1), tn13);
    ta034->exchangeTetraOpposite(t22r->getTetraOpposite(v1), tn13);
    ta123->exchangeTetraOpposite(t22l->getTetraOpposite(v0), tn22);
    ta124->exchangeTetraOpposite(t13->getTetraOpposite(v0), tn22);
    ta134->exchangeTetraOpposite(t22r->getTetraOpposite(v0), tn22);

    v0->cnum -= 2; v1->cnum -= 2;
    tetrasAll.remove(t13); tetrasAll.remove(t22l); tetrasAll.remove(t22r);
    Tetra::destroy(t13); Tetra::destroy(t22l); Tetra::destroy(t22r);

    int time = tn13->vs[3]->time;
    slabSizes[time] -= 1;

    return true;
    // HPC Target [GPU #7]: Batch move32d calls on GPU.
}

void Universe::updateVertexData() {  // Updates vertex connectivity (Sec. 3.2.2)
    vertices.clear();
    int max = 0;
    for (auto v : verticesAll) { vertices.push_back(v); if (v > max) max = v; }  // Collect vertex labels

    vertexNeighbors.clear();
    vertexNeighbors.resize(max + 1);  // Resize for max label

    for (auto v : verticesAll) {  // BFS to find neighbors
        std::vector<Vertex::Label> nbr;
        auto t = v->tetra;
        std::vector<Tetra::Label> current = {t}, next, done;

        do {
            for (auto tc : current) {
                for (auto tcn : tc->tnbr) {
                    if (!tcn->hasVertex(v)) continue;
                    if (std::find(done.begin(), done.end(), tcn) == done.end()) {
                        done.push_back(tcn);
                        next.push_back(tcn);
                    }
                }
            }
            current = next;
            next.clear();
        } while (current.size() > 0);

        for (auto td : done) {  // Extract unique neighbors
            for (auto vd : td->vs) {
                if (std::find(nbr.begin(), nbr.end(), vd) == nbr.end() && vd != v) nbr.push_back(vd);
            }
        }
        vertexNeighbors[v] = nbr;
    }
    // HPC Target [OpenMP #3, GPU #8]: Parallelize BFS for speedup.
}

void Universe::updateHalfEdgeData() {  // Updates half-edge connectivity (Sec. 3.2)
    for (int i = halfEdges.size() - 1; i >= 0; i--) HalfEdge::destroy(halfEdges[i]);  // Clear old half-edges
    assert(HalfEdge::size() == 0);
    halfEdges.clear();

    for (auto t : tetras31) {  // Create half-edges for (3,1)-tetras
        std::array<HalfEdge::Label, 3> these;
        for (int i = 0; i < 3; i++) {
            auto h = HalfEdge::create();
            h->setVertices(t->vs[i], t->vs[(i + 1) % 3]);
            h->tetra = t;
            these[i] = h;
            halfEdges.push_back(h);
        }
        t->setHalfEdges(these[0], these[1], these[2]);
        for (int i = 0; i < 3; i++) { these[i]->next = these[(i + 1) % 3]; these[i]->prev = these[(i - 1 + 3) % 3]; }
    }

    for (auto t : tetras31) {  // Link adjacent half-edges
        for (int i = 0; i < 3; i++) {
            auto v = t->vs[i], vt = t->vs[3];
            auto tc = t->getTetraOpposite(v);
            Tetra::Label tn;
            v = vt;
            while (tc->is22()) {  // Traverse (2,2)-tetras
                tn = tc->getTetraOpposite(v);
                auto vo = v;
                v = tc->vs[2] == v ? tc->vs[3] : tc->vs[2];
                tc = tn;
            }
            assert(tc->is31());
            auto hthis = t->hes[(i + 1) % 3], hthat = tc->getHalfEdgeTo(t->vs[(i + 1) % 3]);
            hthis->adj = hthat; hthat->adj = hthis;
        }
    }
}

void Universe::updateTriangleData() {  // Updates triangle connectivity (Sec. 3.2)
    for (int i = triangles.size() - 1; i >= 0; i--) Triangle::destroy(triangles[i]);  // Clear old triangles
    triangles.clear();

    for (auto t : tetras31) {  // Create triangles for (3,1)-tetras
        auto tr = Triangle::create();
        tr->setVertices(t->vs[0], t->vs[1], t->vs[2]);
        tr->setHalfEdges(t->hes[0], t->hes[1], t->hes[2]);
        t->hes[0]->triangle = tr; t->hes[1]->triangle = tr; t->hes[2]->triangle = tr;
        triangles.push_back(tr);
    }
    triangleNeighbors.resize(triangles.size());

    for (auto tr : triangles) {  // Link triangle neighbors
        tr->setTriangleNeighbors(tr->hes[0]->getAdjacent()->triangle, tr->hes[1]->getAdjacent()->triangle, tr->hes[2]->getAdjacent()->triangle);
        triangleNeighbors[tr] = tr->trnbr;
    }
}

void Universe::updateGeometry() {  // Full geometry update (Sec. 3.2)
    updateVertexData();
    updateHalfEdgeData();
    updateTriangleData();
}

void Universe::check() {  // Diagnostic check (Sec. 1.3)
    printf("====================================================\n");
    printf("c tetras: %d\n", tetrasAll.size());
    int count = 0;

    assert(tetrasAll.size() == Tetra::size());
    for (auto t : tetrasAll) {
        count++;
        for (int i = 0; i < 4; i++) {  // Validate vertices
            assert(verticesAll.contains(t->vs[i]));
            for (int j = i + 1; j < 4; j++) assert(t->vs[i] != t->vs[j]);
        }
        for (int i = 0; i < 4; i++) {  // Validate neighbors
            assert(tetrasAll.contains(t->tnbr[i]));
            assert(t->tnbr[i]->neighborsTetra(t));
            assert(t->tnbr[i] != t && t->tnbr[i] >= 0);
            int sv = 0;
            for (int j = 0; j < 4; j++) if (t->hasVertex(t->tnbr[i]->vs[j])) sv++;
            assert(sv >= 3);

            if (t->is31()) {
                if (i < 3) assert(t->tnbr[i]->is31() || t->tnbr[i]->is22());
                else assert(t->tnbr[i]->is13());
            } else if (t->is13()) {
                if (i == 0) assert(t->tnbr[i]->is31());
                else assert(t->tnbr[i]->is13() || t->tnbr[i]->is22());
            }
        }
        for (int i = 0; i < 4; i++) {  // Validate opposite consistency
            assert(t->getTetraOpposite(t->vs[i]) == t->tnbr[i]);
            assert(t->tnbr[i]->getTetraOpposite(t->getVertexOpposite(t->vs[i])) == t);
        }
    }

    for (auto v : verticesAll) {  // Validate vertex constraints
        assert(tetrasAll.contains(v->tetra));
        if (strictness == 1) assert(v->scnum >= 2);  // Tadpole restriction
        if (strictness == 2) assert(v->scnum >= 3);  // Self-energy restriction
    }

    for (auto tr : triangles) {  // Validate triangle neighbors
        for (auto trn : tr->trnbr) {
            bool found = false;
            for (auto trnn : trn->trnbr) if (trnn == tr) found = true;
            assert(found);
        }
    }
    printf("====================================================\n");
}