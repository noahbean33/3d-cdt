// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "simulation.hpp"  // Simulation class definition
#include "observable.hpp"  // Observable class definition

// Static member initializations (Sec. 3)
double Simulation::k0 = 0;        // Gravitational coupling (Sec. 2.3)
double Simulation::k3 = 0;        // Cosmological coupling (Sec. 2.3)
std::default_random_engine Simulation::rng(0);  // RNG with default seed 0
// TODO(JorenB): seed properly - Currently fixed; updated in start()
// HPC Target [MPI #4]: Static RNG needs per-rank instances.
int Simulation::targetVolume = 0;    // Target total volume (Sec. 2.4)
int Simulation::target2Volume = 0;   // Secondary target volume
double Simulation::epsilon = 0.00004;// Volume fixing strength (Sec. 2.4)
bool Simulation::measuring = false;  // Measurement phase flag (Sec. 3.3)
std::vector<Observable*> Simulation::observables3d;  // 3D observables
std::vector<Observable*> Simulation::observables2d;  // 2D observables
std::array<int, 3> Simulation::moveFreqs = {0, 0, 0}; // Move frequencies

void Simulation::start(double k0_, double k3_, int sweeps, int thermalSweeps, int ksteps, int targetVolume_, int target2Volume_, int seed, std::string OutFile, int v1, int v2, int v3) {
    Simulation::moveFreqs = {v1, v2, v3};  // Set move frequencies
    // Comment: v1, v2, v3 likely weight add/delete, flip, shift moves (Sec. 2.3).

    targetVolume = targetVolume_;  // Set target volumes
    target2Volume = target2Volume_;
    k3 = k3_;  // Set couplings
    k0 = k0_;

    for (auto o : observables3d) o->clear();  // Reset 3D observables
    for (auto o : observables2d) o->clear();  // Reset 2D observables
    // Comment: Clears prior measurement data (Sec. 3.4).

    rng.seed(seed);  // Set RNG seed
    // Comment: Overrides default seed 0 (Sec. 2.1).

    measuring = true;  // Initially true; toggled later
    //////////////////////////////////////////////////////////////////////
    // ********************** START THERMAL SWEEPS ******************** //
    //////////////////////////////////////////////////////////////////////
    printf("k0: %g, k3: %g, epsilon: %g \t thermal: %d \t sweeps: %d Target: %d\t Target2d: %d\t \n", k0, k3, epsilon, thermalSweeps, sweeps, targetVolume, target2Volume);

    for (int i = 1; i <= thermalSweeps; i++) {  // Thermalization phase (Sec. 3.3.2)
        int total2v = 0;  // Total 2D volume (sum of slice sizes)
        for (auto ss : Universe::sliceSizes) total2v += ss;

        double n31 = Universe::tetras31.size();  // Number of (3,1)-tetras
        int n3 = Universe::tetrasAll.size();     // Total tetra count

        printf("Thermal: i: %d\t  Tetra::size: %d tetras31:  %g k3: %g \n", i, n3, n31, k3);

        performSweep(ksteps * 1000);  // Perform sweep (Sec. 3.3.2)
        // HPC Target [OpenMP #1]: Parallelize this loop.

        tune();  // Adjust k3 (Sec. 3.3.1)

        if (i % (thermalSweeps / 10) == 0)  // Periodic export
            Universe::exportGeometry(OutFile);

        prepare();  // Update geometry (Sec. 3.2)
        for (auto o : observables3d) o->measure();  // Measure 3D observables
        // HPC Target [OpenMP #2]: Parallelize measurement loop.
    }

    //////////////////////////////////////////////////////////////////////
    // ********************** END THERMAL SWEEPS ******************** //
    //////////////////////////////////////////////////////////////////////
    printf("======\n");
    // ********************** START MEASURE SWEEPS ****************** //
    printf("k0: %g, k3: %g, epsilon: %g\n", k0, k3, epsilon);

    for (int i = 1; i <= sweeps; i++) {  // Measurement phase (Sec. 3.3.3)
        int total2v = 0;
        for (auto ss : Universe::sliceSizes) total2v += ss;
        int avg2v = total2v / Universe::nSlices;  // Average slice size

        printf("SWEEPS: i: %d\t Target: %d\t Target2d: %d\t CURRENT: %d avgslice: %d\n", i, targetVolume, target2Volume, Tetra::size(), avg2v);

        performSweep(ksteps * 1000);  // Perform sweep
        // HPC Target [OpenMP #1]: Parallelize this loop.

        if (i % (sweeps / 10) == 0)  // Periodic export
            Universe::exportGeometry(OutFile);

        if (observables3d.size() > 0) {  // Measure 3D observables
            int vol_switch = Universe::volfix_switch;
            if (targetVolume > 0) {  // Volume fixing loop
                int cmp;
                do {
                    attemptMove();
                    cmp = vol_switch == 0 ? Universe::tetras31.size() : Universe::tetrasAll.size();
                } while (cmp != targetVolume);  // Adjust until target hit
            }

            prepare();  // Update geometry
            // HPC Target [OpenMP #3]: Parallelize if BFS involved.
            for (auto o : observables3d) o->measure();  // Measure
            // HPC Target [OpenMP #2]: Parallelize this loop.
        }

        if (target2Volume > 0) {  // Measure 2D observables
            bool hit = false;
            do {
                attemptMove();
                for (auto s : Universe::sliceSizes) {
                    if (s == target2Volume) { hit = true; break; }
                }
            } while (!hit);  // Adjust until a slice matches target2Volume

            prepare();
            for (auto o : observables2d) o->measure();
            // HPC Target [OpenMP #2]: Parallelize this loop.
        }
    }
    // HPC Targets: [MPI #4] (distribute runs), [GPU #7] (accelerate moves), [General #12] (tune sweep size)
}

int Simulation::attemptMove() {  // Attempts a random move (Sec. 2.1)
    std::array<int, 3> cumFreqs = {moveFreqs[0], moveFreqs[0] + moveFreqs[1], moveFreqs[0] + moveFreqs[1] + moveFreqs[2]};
    int freqTotal = moveFreqs[0] + moveFreqs[1] + moveFreqs[2];  // Total frequency

    std::uniform_int_distribution<> moveGen(0, freqTotal - 1);  // Move selector
    std::uniform_int_distribution<> binGen(0, 1);               // Binary selector

    int move = moveGen(rng);  // Random move index

    if (move < cumFreqs[0]) {  // Add/Delete moves (Sec. 2.3.1)
        if (binGen(rng) == 0) {
            if (moveAdd()) return 1; else return -1;  // 1: success, -1: fail
        } else {
            if (moveDelete()) return 2; else return -2;
        }
    } else if (cumFreqs[0] <= move && move < cumFreqs[1]) {  // Flip move (Sec. 2.3.2)
        if (moveFlip()) return 3; else return -3;
    } else if (cumFreqs[1] <= move) {  // Shift/Ishift moves (Sec. 2.3.3)
        if (binGen(rng) == 0) {
            if (binGen(rng) == 0) {
                if (moveShift()) return 4; else return -4;  // Upward shift
            } else {
                if (moveShiftD()) return 4; else return -4; // Downward shift
            }
        } else {
            if (binGen(rng) == 0) {
                if (moveShiftI()) return 5; else return -5;  // Upward ishift
            } else {
                if (moveShiftID()) return 5; else return -5; // Downward ishift
            }
        }
    }
    return 0;  // No move attempted (shouldn’t occur with valid freqTotal)
    // HPC Target [GPU #7]: Batch move attempts on GPU.
}

std::vector<int> Simulation::performSweep(int n) {  // Executes n move attempts (Sec. 3.3.2, 3.3.3)
    std::vector<int> moves(6, 0);      // Successful moves (0 unused, 1-5 for types)
    std::vector<int> failed_moves(6, 0); // Failed moves

    for (int i = 0; i < n; i++) {  // Attempt n moves
        int move_num = attemptMove();
        int move = abs(move_num);  // Move type (1-5)
        moves[move]++;             // Count attempt
        if (move_num < 0) failed_moves[move]++;  // Count failure
    }
    // HPC Target [OpenMP #1]: Parallelize this loop.

    int m1 = moves[1] + moves[2];  // Add/Delete total
    int m2 = moves[3];             // Flip total
    int m3 = moves[4] + moves[5];  // Shift/Ishift total
    int f1 = failed_moves[1] + failed_moves[2];  // Failed Add/Delete
    int f2 = failed_moves[3];                    // Failed Flip
    int f3 = failed_moves[4] + failed_moves[5];  // Failed Shift/Ishift

    // Avoid division by zero (bug workaround?)
    if (m1 == f1) m1++, m1++, f1++;
    if (m2 == f2) m2++, m2++, f2++;
    if (m3 == f3) m3++, m3++, f3++;

    double r1 = static_cast<double>(m1) / f1;  // Success ratios
    double r2 = static_cast<double>(m2) / f2;
    double r3 = static_cast<double>(m3) / f3;

    printf("%g\t%g\t%g\t\n", r1, r2, r3);  // Log ratios
    return moves;  // Return move counts
    // HPC Target [General #12]: Adjust n dynamically based on r1, r2, r3.
}

bool Simulation::moveAdd() {  // Attempts (2,6)-move (Sec. 2.3.1)
    double n31 = Universe::tetras31.size();
    int n3 = Universe::tetrasAll.size();
    int vol_switch = Universe::volfix_switch;

    double edS = exp(1 * k0 - 4 * k3);  // Action change (Sec. 2.3.1, eq. 25)
    double rg = n31 / (n31 + 2.0);     // Selection probability ratio
    double ar = edS * rg;              // Acceptance ratio

    if (vol_switch == 0) {  // Fix (3,1) count
        if (targetVolume > 0) ar *= exp(4 * epsilon * (targetVolume - n31 - 1));  // Volume fixing term
    } else {  // Fix total count
        if (targetVolume > 0) ar *= exp(8 * epsilon * (targetVolume - n3 - 2));
    }

    if (ar < 1.0) {  // Metropolis acceptance (Sec. 2.2)
        std::uniform_real_distribution<> uniform(0.0, 1.0);
        if (uniform(rng) > ar) return false;
    }

    Tetra::Label t = Universe::tetras31.pick();  // Random (3,1)-tetra
    Universe::move26(t);  // Perform move
    return true;
    // HPC Target [GPU #7]: Batch on GPU.
}

bool Simulation::moveDelete() {  // Attempts (6,2)-move (Sec. 2.3.1)
    double n31 = Universe::tetras31.size();
    int n3 = Universe::tetrasAll.size();
    int vol_switch = Universe::volfix_switch;

    double edS = exp(-1 * k0 + 4 * k3);  // Action change (Sec. 2.3.1, eq. 27)
    double rg = n31 / (n31 - 2.0);      // Selection probability ratio
    double ar = edS * rg;

    if (vol_switch == 0) {
        if (targetVolume > 0) ar *= exp(-4 * epsilon * (targetVolume - n31 - 1));
    } else {
        if (targetVolume > 0) ar *= exp(-8 * epsilon * (targetVolume - n3 - 2));
    }

    if (ar < 1.0) {
        std::uniform_real_distribution<> uniform(0.0, 1.0);
        if (uniform(rng) > ar) return false;
    }

    Vertex::Label v = Universe::verticesAll.pick();  // Random vertex
    if (v->cnum != 6 || v->scnum != 3) return false; // Check move condition (Sec. 2.3.1)

    Universe::move62(v);
    return true;
    // HPC Target [GPU #7]: Batch on GPU.
}

bool Simulation::moveFlip() {  // Attempts (4,4)-move (Sec. 2.3.2)
    Tetra::Label t012 = Universe::tetras31.pick();  // Random (3,1)-tetra
    std::uniform_int_distribution<> neighborGen(0, 2);
    Tetra::Label t230 = t012->tnbr[neighborGen(rng)];  // Random spatial neighbor

    if (!t230->is31()) return false;  // Must be (3,1)
    if (!t012->tnbr[3]->neighborsTetra(t230->tnbr[3])) return false;  // Check vertical neighbors

    return Universe::move44(t012, t230);  // No Metropolis step (ar=1, Sec. 2.3.2)
    // HPC Target [GPU #7]: Batch on GPU.
}

bool Simulation::moveShift() {  // Attempts (2,3)-move upward (Sec. 2.3.3)
    double edS = exp(-1.0 * k3);  // Action change (Sec. 2.3.3, eq. 28)
    double rg = 1.0;              // Selection probability ratio
    double ar = edS * rg;
    int n3 = Universe::tetrasAll.size();
    int vol_switch = Universe::volfix_switch;

    if (vol_switch == 1 && targetVolume > 0) ar *= exp(epsilon * (2 * targetVolume - 2 * n3 - 1));

    if (ar < 1.0) {
        std::uniform_real_distribution<> uniform(0.0, 1.0);
        if (uniform(rng) > ar) return false;
    }

    Tetra::Label t = Universe::tetras31.pick();
    std::uniform_int_distribution<> neighborGen(0, 2);
    Tetra::Label tn = t->tnbr[neighborGen(rng)];

    if (!tn->is22()) return false;  // Must be (2,2)-tetra

    return Universe::move23u(t, tn);
    // HPC Target [GPU #7]: Batch on GPU.
}

bool Simulation::moveShiftD() {  // Attempts (2,3)-move downward
    double edS = exp(-1.0 * k3);
    double rg = 1.0;
    double ar = edS * rg;
    int n3 = Universe::tetrasAll.size();
    int vol_switch = Universe::volfix_switch;

    if (vol_switch == 1 && targetVolume > 0) ar *= exp(epsilon * (2 * targetVolume - 2 * n3 - 1));

    if (ar < 1.0) {
        std::uniform_real_distribution<> uniform(0.0, 1.0);
        if (uniform(rng) > ar) return false;
    }

    Tetra::Label tv = Universe::tetras31.pick();
    auto t = tv->tnbr[3];  // Vertical (1,3)-tetra
    std::uniform_int_distribution<> neighborGen(1, 3);
    Tetra::Label tn = t->tnbr[neighborGen(rng)];

    if (!tn->is22()) return false;

    return Universe::move23d(t, tn);
    // HPC Target [GPU #7]: Batch on GPU.
}

bool Simulation::moveShiftI() {  // Attempts (3,2)-move upward (Sec. 2.3.3)
    double edS = exp(1 * k3);  // Action change (Sec. 2.3.3, eq. 29)
    double rg = 1.0;
    double ar = edS * rg;
    int n3 = Universe::tetrasAll.size();
    int vol_switch = Universe::volfix_switch;

    if (vol_switch == 1 && targetVolume > 0) ar *= exp(-epsilon * (2 * targetVolume - 2 * n3 - 1));

    if (ar < 1.0) {
        std::uniform_real_distribution<> uniform(0.0, 1.0);
        if (uniform(rng) > ar) return false;
    }

    Tetra::Label t = Universe::tetras31.pick();
    std::uniform_int_distribution<> neighborGen(0, 2);
    int neighbor = neighborGen(rng);
    Tetra::Label t22l = t->tnbr[neighbor], t22r = t->tnbr[(neighbor + 2) % 3];

    if (!t22l->is22() || !t22r->is22() || !t22l->neighborsTetra(t22r)) return false;

    int sv = 0;  // Shared vertices
    for (int i = 0; i < 4; i++) if (t22r->hasVertex(t22l->vs[i])) sv++;
    if (sv != 3) return false;  // Must share a face

    return Universe::move32u(t, t22l, t22r);
    // HPC Target [GPU #7]: Batch on GPU.
}

bool Simulation::moveShiftID() {  // Attempts (3,2)-move downward
    double edS = exp(1 * k3);
    double rg = 1.0;
    double ar = edS * rg;
    int n3 = Universe::tetrasAll.size();
    int vol_switch = Universe::volfix_switch;

    if (vol_switch == 1 && targetVolume > 0) ar *= exp(-epsilon * (2 * targetVolume - 2 * n3 - 1));

    if (ar < 1.0) {
        std::uniform_real_distribution<> uniform(0.0, 1.0);
        if (uniform(rng) > ar) return false;
    }

    Tetra::Label t = Universe::tetras31.pick()->tnbr[3];  // Vertical (1,3)-tetra
    std::uniform_int_distribution<> neighborGen(0, 2);
    int neighbor = neighborGen(rng);
    Tetra::Label t22l = t->tnbr[1 + neighbor], t22r = t->tnbr[1 + (neighbor + 2) % 3];

    if (!t22l->is22() || !t22r->is22() || !t22l->neighborsTetra(t22r)) return false;

    int sv = 0;
    for (int i = 0; i < 4; i++) if (t22r->hasVertex(t22l->vs[i])) sv++;
    if (sv != 3) return false;

    return Universe::move32d(t, t22l, t22r);
    // HPC Target [GPU #7]: Batch on GPU.
}

void Simulation::prepare() {  // Prepares for measurements (Sec. 3.3.3)
    Universe::updateGeometry();  // Updates connectivity (Sec. 3.2)
    // HPC Target [OpenMP #3]: Parallelize if BFS involved.
}

void Simulation::tune() {  // Tunes k3 to pseudocritical value (Sec. 3.3.1)
    double delta_k3 = 0.000001;  // Step size
    double ratio = 100;          // Scaling factor (unused in code?)

    int border_far = targetVolume * 0.5;      // Large deviation threshold
    int border_close = targetVolume * 0.05;   // Close threshold
    int border_vclose = targetVolume * 0.002; // Very close threshold
    int border_vvclose = targetVolume * 0.0001; // Extremely close threshold

    int vol_switch = Universe::volfix_switch;
    int fixvolume = vol_switch == 0 ? Universe::tetras31.size() : Universe::tetrasAll.size();

    int diff = targetVolume - fixvolume;  // Volume deviation
    if (diff > border_far) k3 -= delta_k3 * 1000;         // Large decrease
    else if (diff < -border_far) k3 += delta_k3 * 1000;   // Large increase
    else if (diff > border_close) k3 -= delta_k3 * 1000;  // Moderate decrease
    else if (diff < -border_close) k3 += delta_k3 * 1000; // Moderate increase
    else if (diff > border_vclose) k3 -= delta_k3 * 100;  // Small decrease
    else if (diff < -border_vclose) k3 += delta_k3 * 100; // Small increase
    else if (diff > border_vvclose) k3 -= delta_k3 * 20;  // Tiny decrease
    else if (diff < -border_vvclose) k3 += delta_k3 * 20; // Tiny increase
    // Comment: Adaptive tuning based on volume deviation (Sec. 3.3.1).
}