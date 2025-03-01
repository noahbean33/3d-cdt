// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <random>       // For std::default_random_engine
#include <vector>       // For std::vector (observables, sweep results)
#include <string>       // For std::string (outFile)
#include "universe.hpp" // Universe class for triangulation management

// Comment: Forward declaration of Observable class (defined elsewhere, e.g., observable.hpp).
class Observable;

// Comment: Simulation manages the Monte Carlo simulation process for 3D CDT (Sec. 2, 3.3).
class Simulation {
public:
    static double lambda;  // Cosmological coupling constant (Sec. 2.3)
    // Comment: Likely a legacy or misnamed variable; paper uses k3 for cosmological coupling.

    static void start(double k0, double k3, int sweeps, int thermalSweeps, int ksteps, int targetVolume_, int target2Volume, int seed, std::string outFile, int v1, int v2, int v3);
    // Comment: Initiates the simulation with parameters (Sec. 3.3):
    // - k0: Gravitational coupling (1/G_N, Sec. 2.3)
    // - k3: Cosmological coupling (Sec. 2.3)
    // - sweeps: Measurement sweeps (Sec. 3.3.3)
    // - thermalSweeps: Thermalization sweeps (Sec. 3.3.2)
    // - ksteps: Tuning steps (Sec. 3.3.1)
    // - targetVolume_, target2Volume: Target volumes (Sec. 2.4)
    // - seed: RNG seed (Sec. 2.1)
    // - outFile: Output file name
    // - v1, v2, v3: Custom parameters (possibly move weights or debug flags)
    // HPC Targets: [OpenMP #1], [MPI #4], [GPU #7], [General #12]

    static void addObservable3d(Observable& o) {
        observables3d.push_back(&o);  // Adds a 3D observable
    }
    // Comment: Registers a 3D observable (e.g., VolumeProfile, Sec. 3.4).

    static void addObservable2d(Observable& o) {
        observables2d.push_back(&o);  // Adds a 2D observable
    }
    // Comment: Registers a 2D observable (e.g., Ricci2d, Sec. 3.4); unexpected in 3D context (see main.cpp note).

    static std::array<int, 3> moveFreqs;  // Frequencies of move attempts
    // Comment: Likely tracks add/delete, flip, shift move counts (Sec. 2.3).

    static int attemptMove();  // Attempts a random Monte Carlo move
    // Comment: Selects and tries a move (Sec. 2.3); returns move type or success indicator.
    // HPC Target: [GPU #7]

    static int targetVolume;   // Target total volume (Sec. 2.4)
    static int target2Volume;  // Secondary target volume (possibly for tuning)

    static int thermal;  // Thermalization status or counter
    // Comment: Tracks thermalization phase (Sec. 3.3.2).

    static double k3_s;  // Stored k3 value (starting or tuned?)
    static double k0;    // Gravitational coupling (Sec. 2.3)
    static double k3;    // Cosmological coupling (Sec. 2.3)
    // Comment: Simulation parameters set by start().

private:
    static std::default_random_engine rng;  // RNG for move selection
    // Comment: Randomness source (Sec. 2.1).
    // HPC Target: [MPI #4]

    static double epsilon;  // Volume fixing strength (Sec. 2.4)
    // Comment: Controls fluctuation size in S_fix = ε(N - N̄)^2.

    static bool measuring;  // Flag for measurement phase
    // Comment: Distinguishes thermalization vs. measurement (Sec. 3.3).

    static std::vector<Observable*> observables3d;  // 3D observables list
    static std::vector<Observable*> observables2d;  // 2D observables list
    // Comment: Stores pointers to observables (Sec. 3.4).
    // HPC Target: [OpenMP #2]

    static std::vector<int> performSweep(int n);  // Executes n move attempts
    // Comment: Performs a sweep (Sec. 3.3.2, 3.3.3); returns move stats.
    // HPC Target: [OpenMP #1]

    // Move attempt methods (Sec. 2.3)
    static bool moveAdd();     // (2,6)-move (Sec. 2.3.1)
    static bool moveDelete();  // (6,2)-move (Sec. 2.3.1)
    static bool moveFlip();    // (4,4)-move (Sec. 2.3.2)
    static bool moveShift();   // (2,3)-move upward (Sec. 2.3.3)
    static bool moveShiftD();  // (2,3)-move downward
    static bool moveShiftI();  // (3,2)-move upward (ishift)
    static bool moveShiftID(); // (3,2)-move downward (ishift)
    // Comment: Wrappers for Universe moves; return success status.
    // HPC Target: [GPU #7]

    static void prepare();  // Prepares for measurements
    // Comment: Likely calls Universe::updateGeometry (Sec. 3.2, 3.3.3).
    // HPC Target: [OpenMP #3]

    static void tune();  // Tunes coupling constants
    // Comment: Adjusts k3 to pseudocritical value (Sec. 3.3.1).
};

// HPC Targets Summary:
// [OpenMP #1]: Parallelize performSweep() for sweeps (Sec. 3.3).
// [OpenMP #2]: Parallelize observable measurements (Sec. 3.4).
// [OpenMP #3]: Parallelize prepare() if BFS involved (Sec. 3.2).
// [MPI #4]: Distribute start() across nodes; fix static rng.
// [GPU #7]: Accelerate move attempts in attemptMove(), move*().
// [General #12]: Dynamically tune sweep size in performSweep().