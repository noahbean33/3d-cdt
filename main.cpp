// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#include "config.hpp"        // Configuration file reader
#include "universe.hpp"      // Manages the CDT triangulation state
#include "simulation.hpp"    // Controls Monte Carlo simulation workflow

#include "observables/volume_profile.hpp"  // Observable for volume distribution
#include "observables/ricci2d.hpp"         // 2D Ricci curvature observable (note: unexpected in 3D context)

std::default_random_engine rng(1);  // Global random number generator with fixed seed 1
// Comment: Provides randomness for Monte Carlo moves (Sec. 2.1 of paper). Fixed seed ensures reproducibility.
// HPC Target [MPI #4]: Static RNG needs per-rank instances for distributed runs to avoid identical sequences.

int main(int argc, const char * argv[]) {
    std::string fname;
    if (argc > 1) {
        fname = std::string(argv[1]);  // Read config file name from command line argument
        printf("%s\n", fname.c_str()); // Print config file name
    }
    // Comment: Allows user to specify config file via command line (e.g., ./main config.txt).
    // HPC Target [General #10]: Minimal I/O impact here, but consider buffered output for large runs.

    ConfigReader cfr;              // Object to parse configuration file
    cfr.read(fname);               // Load parameters from file

    // Read simulation parameters from config file (Sec. 3.3 of paper)
    double k0 = cfr.getDouble("k0");              // Gravitational coupling (Sec. 2.3)
    double k3_s = cfr.getDouble("k3");            // Cosmological coupling (Sec. 2.3)
    int genus = cfr.getInt("genus");              // Topology genus (Sec. 1.3)
    int targetVolume = cfr.getInt("targetvolume"); // Target total volume (Sec. 2.4)
    int target2Volume = cfr.getInt("target2volume"); // Target volume for secondary fixing (if applicable)
    int volfixSwitch = cfr.getInt("volfixswitch");   // Volume fixing toggle (Sec. 2.4)
    int seed = cfr.getInt("seed");                // Random seed for simulation
    std::string outputDir = cfr.getString("outputdir"); // Directory for output files
    std::string fID = cfr.getString("fileid");    // Unique file identifier
    int thermalSweeps = cfr.getInt("thermalsweeps");  // Number of thermalization sweeps (Sec. 3.3.2)
    int sweeps = cfr.getInt("measuresweeps");    // Number of measurement sweeps (Sec. 3.3.3)
    int kSteps = cfr.getInt("ksteps");            // Steps for coupling tuning (Sec. 3.3.1)
    int strictness = cfr.getInt("strictness");    // Manifold strictness parameter (Sec. 1.3)
    int v1 = cfr.getInt("v1");                    // Custom parameter (possibly vertex-related)
    int v2 = cfr.getInt("v2");                    // Custom parameter
    int v3 = cfr.getInt("v3");                    // Custom parameter
    std::string inFile = cfr.getString("infile"); // Input triangulation file
    std::string outFile = cfr.getString("outfile"); // Output file for results
    // Comment: Parameters configure simulation per Sec. 3 of paper. 'v1-v3' unclear without code context; possibly move-specific.

    printf("fID: %s\n", fID.c_str());    // Print file ID
    printf("seed: %d\n", seed);          // Print seed
    printf("strictness: %d\n", strictness); // Print strictness level
    // Comment: Diagnostic output for key parameters.

    Observable::data_dir = outputDir;    // Set output directory for observables (Sec. 3.4)
    // Comment: Static member suggests shared state across observables; potential concurrency issue.

    Universe::initialize(inFile, fID, strictness, volfixSwitch);
    // Comment: Sets up initial triangulation state (Sec. 3). Reads from 'inFile' if provided, enforces manifold properties.
    // HPC Target [General #10]: Initialization could pre-allocate memory (e.g., Pool capacity) for cache efficiency.

    printf("\n\n#######################\n");
    printf("* * * Initialized * * *\n");  // Status message
    printf("#######################\n\n");

    VolumeProfile vp3(fID);              // Observable for 3D volume profile (Sec. 3.4)
    Simulation::addObservable3d(vp3);    // Register with 3D simulation
    // Comment: Tracks volume distribution across time slices (S^1 x S^2 topology in Sec. 2.3).

    Ricci2d ricci2d(fID);                // Observable for 2D Ricci curvature
    Simulation::addObservable2d(ricci2d); // Register with 2D simulation (unexpected in 3D context)
    // Comment: Appears to be a mismatch; paper focuses on 3D CDT (Sec. 2.3), but 'Ricci2d' suggests 2D observable.
    // Potential Bug: Should this be a 3D Ricci observable? Verify against codebase intent.

    Simulation::start(k0, k3_s, sweeps, thermalSweeps, kSteps, targetVolume, target2Volume, seed, outFile, v1, v2, v3);
    // Comment: Runs full simulation: tuning (Sec. 3.3.1), thermalization (Sec. 3.3.2), and measurement (Sec. 3.3.3).
    // HPC Targets:
    // [OpenMP #1]: Parallelize sweeps within 'start()' (e.g., lines 97-105 in HPC plan).
    // [MPI #4]: Distribute independent runs across nodes here.
    // [GPU #7]: Accelerate move attempts in sweeps.
    // [General #12]: Tune sweep size dynamically based on acceptance rates.

    printf("\n\n####################\n");
    printf("* * * Finished * * *\n");  // Status message
    printf("####################\n\n");

    printf("t31: %d\n", Universe::tetras31.size()); // Print number of (3,1)-tetrahedra
    // Comment: Diagnostic output for triangulation state (Sec. 2.3). Useful for verifying volume fixing.

    return 0;  // Exit successfully
}