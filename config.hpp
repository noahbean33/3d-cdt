// Copyright 2021 Joren Brunekreef, Daniel Nemeth and Andrzej Görlich
#pragma once  // Prevents multiple inclusions
// Comment: Standard header guard.

#include <fstream>       // For std::ifstream to read config file
#include <string>        // For std::string (keys, values, filename)
#include <cassert>       // For assert() to enforce required keys
#include <unordered_map> // For std::unordered_map to store key-value pairs

// Comment: ConfigReader reads simulation parameters from a file for 3D CDT (Sec. 3).
class ConfigReader {
public:
    void read(std::string fname) {  // Reads config file into dictionary
        std::ifstream infile(fname);  // Open file
        assert(infile.is_open());     // Ensure file opened successfully
        std::string key, value;

        while (infile >> key >> value) {  // Read key-value pairs
            dict[key] = value;            // Store in map
        }

        // Enforce required keys (Sec. 3)
        assert(dict.find("k0") != dict.end());           // Gravitational coupling (Sec. 2.3)
        assert(dict.find("k3") != dict.end());           // Cosmological coupling (Sec. 2.3)
        assert(dict.find("genus") != dict.end());        // Topology genus (Sec. 1.3)
        assert(dict.find("targetvolume") != dict.end()); // Target total volume (Sec. 2.4)
        assert(dict.find("target2volume") != dict.end()); // Secondary target volume
        assert(dict.find("volfixswitch") != dict.end());  // Volume fixing toggle (Sec. 2.4)
        assert(dict.find("seed") != dict.end());         // RNG seed (Sec. 2.1)
        assert(dict.find("outputdir") != dict.end());    // Output directory
        assert(dict.find("fileid") != dict.end());       // File identifier
        assert(dict.find("thermalsweeps") != dict.end()); // Thermalization sweeps (Sec. 3.3.2)
        assert(dict.find("measuresweeps") != dict.end()); // Measurement sweeps (Sec. 3.3.3)
        assert(dict.find("ksteps") != dict.end());       // Tuning steps (Sec. 3.3.1)
        assert(dict.find("strictness") != dict.end());   // Manifold strictness (Sec. 1.3)
        assert(dict.find("v1") != dict.end());           // Custom parameter (e.g., move frequency)
        assert(dict.find("v2") != dict.end());           // Custom parameter
        assert(dict.find("v3") != dict.end());           // Custom parameter
        assert(dict.find("infile") != dict.end());       // Input geometry file
        assert(dict.find("outfile") != dict.end());      // Output geometry file
        // Comment: Called by main.cpp to load parameters; asserts ensure all keys present.
    }

    int getInt(std::string key) {  // Retrieves integer value
        return std::stoi(dict[key]);
    }
    // Comment: Converts string value to int (e.g., "seed", "thermalsweeps").

    double getDouble(std::string key) {  // Retrieves double value
        return std::stod(dict[key]);
    }
    // Comment: Converts string value to double (e.g., "k0", "k3").

    std::string getString(std::string key) {  // Retrieves string value
        return dict[key];
    }
    // Comment: Returns string as-is (e.g., "outputdir", "fileid").

private:
    std::unordered_map<std::string, std::string> dict;  // Stores key-value pairs
    // Comment: Hash map for efficient lookup (Sec. 3).
};

// HPC Targets Summary:
// [General #10]: Pre-allocating dict with reserve() could improve cache efficiency for large configs.
// No direct parallelization targets (OpenMP #3, GPU #8) as this is I/O setup, executed once.