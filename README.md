# Monte Carlo Simulation of (2+1)-Dimensional Causal Dynamical Triangulations

This codebase samples the partition sum of (2+1)-dimensional Causal Dynamical Triangulations (CDT), a lattice model of 3D quantum gravity with \( S^1 \times S^2 \) topology. It implements Monte Carlo methods to simulate random geometries, supporting custom observables for computing quantum expectation values. Simulation parameters are set via a config file. This fork enhances the original codebase with high-performance computing (HPC) and machine learning (ML) upgrades.

## Usage
### Build with GNU Make:
```bash
make
```
### Run an example simulation from the example directory:
```bash
./run.sh
```

### Config file template (all parameters required):
```
k0              2.0
k3              0.5
targetVolume    100000
slices          50
seed            1
fileID          3d-test-100k-1
measurements    50
importGeom      true
```
- **k0**: Inverse bare Newton constant coupling (proportional to \(1/G_N\), Sec. 2.3).
- **k3**: Bare cosmological constant coupling.
- **targetVolume**: Target number of tetrahedra (\(N_3\)).
- **slices**: Number of discrete time slices.
- **seed**: Random number generator seed (fixed output for identical seeds).
- **fileID**: Identifier for observable output files.
- **measurements**: Number of measurements per observable.
- **importGeom**: Import a pre-existing geometry from `geom/` if available.

## Observables
Standard observables (e.g., volume profile, shell volumes for Hausdorff dimension) are in `observables/`. Add them in `main.cpp` as shown in the example. Custom observables can leverage `Universe` (access to `Vertex` and `Tetra` objects, including neighbor lists) and `Observable` methods (e.g., metric spheres, dual distances).

## Optimization Plan (2025)

### Timeline

## Citation
If you use or extend this codebase, please cite the original work:
```bibtex
@misc{brunekreef2023simulating,
  title = {Simulating {{CDT}} quantum gravity},
  author = {Brunekreef, Joren and G{"o}rlich, Andrzej and Loll, Renate},n  year = {2023},
  month = oct,
  number = {arXiv:2310.16744},
  eprint = {2310.16744},
  primaryclass = {gr-qc, physics:hep-lat, physics:hep-th},
  publisher = {{arXiv}},
  doi = {10.48550/arXiv.2310.16744},
  urldate = {2023-11-13},
  archiveprefix = {arxiv},
  keywords = {General Relativity and Quantum Cosmology, High Energy Physics - Lattice, High Energy Physics - Theory}
}
```
Please also acknowledge this fork: `https://github.com/noahbean33/3d-cdt`.
