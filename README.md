# PFFM

- **`apps/main_modeI.cpp`** contains C++ code implementing a finite-element method (FEM) model of a double cantilever beam, which uses the phase-field fracture model for interface fracture propagation.

- **`apps/main_modeII.cpp`** contains C++ code implementing a finite-element method (FEM) model of an end loaded split test, which uses the phase-field fracture model for interface fracture propagation.

- The diffuse phase-field damage interacts with the surrounding bulk material, which artificially increases the apparent interface fracture toughness (assuming the bulk material has a higher fracture toughness). This effect can be mitigated by using an effective fracture toughness value.

- **`matlab/BVPSolverForGceff.m`** is a MATLAB script that pre-computes the effective fracture toughness for pure mode-I PFFM simulations using various methods.

- **`matlab/BVPSolverForGceff_modeII.m`** *(Forthcoming)*: The script for pure mode-II effective toughness is currently being finalized for submission and will be released to this repository upon publication.

- **`legacy/`**: Contains the original, monolithic version of the solver scripts (Version 1). These are preserved for reference but are deprecated in favor of the optimised C++ architecture in apps/.

The associated academic journal article for the mode-I work is available [here](https://doi.org/10.1016/j.engfracmech.2025.111546).

## Requirements

To compile the C++ code, you must include the [Eigen 3.4.0 header-only C++ library](https://gitlab.com/libeigen/eigen/-/releases/3.4.0), which provides linear algebra functionalities.
More information about Eigen can be found on the [Eigen website](https://eigen.tuxfamily.org/index.php?title=Main_Page).

*(Optional)* SuiteSparse (CHOLMOD) is recommended for faster execution but is not required; the build system will automatically fall back to Eigen's internal solvers if it is not found. See the [SuiteSparse](./SuiteSparse_setup.txt) Guide for details.

## Build & Run

The C++ code has been successfully compiled on:
- **Ubuntu 20.04.6 LTS** using **GCC 9.4.0**

The project includes an automated script (`run.sh`) to handle CMake configuration, compilation, and execution.

### 1. Run Mode I (Double Cantilever Beam)

```bash
./run.sh I
```

### 2. Run Mode II (End Loaded Split)

```bash
./run.sh II
```

### Manual Build
If you prefer to build manually without the helper script:

```bash
mkdir build && cd build
cmake ..
make
./solver_modeI
```

**Note:**
For best performance, the CMake configuration automatically applies the following compile flags:
```
-O3 -fno-math-errno -DNDEBUG -march=native
```
These flags enable high optimisation, disable unnecessary math error checking, and allow Eigen to fully exploit your hardware’s capabilities (SIMD/AVX).

## Usage

Model parameters (geometry, material properties, solver settings, etc.) can be modified directly in the source code near the start of the main function in the `apps/` files. After changing these parameters, running the `run.sh` script will automatically recompile the code and save a snapshot of the parameters alongside the results for reproducibility.

## License

This project is licensed under the [MIT License](./LICENSE).  
See the [LICENSE](./LICENSE) file for details.

## Citation

C. A. Fear, S. Wang, C. M. Harvey (2025), [Effective fracture toughness in phase-field models for interface fracture](https://doi.org/10.1016/j.engfracmech.2025.111546), Engineering Fracture Mechanics (328), 111546.