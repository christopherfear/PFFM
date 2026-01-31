cmake_minimum_required(VERSION 3.10)
project(PFFM_Solver)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# 1. Settings
set(CMAKE_CXX_STANDARD 14)
# Performance Flags:
# -O3: Max optimisation
# -fno-math-errno: Faster math (skips setting error codes)
# -DNDEBUG: Disables Eigen range checks (CRITICAL for speed)
# -march=native: Uses CPU's specific instructions (AVX, etc.)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fno-math-errno -DNDEBUG -march=native")

# 2. Include Paths
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories($ENV{HOME}/PFFM/eigen-3.4.0) 
include_directories(/opt/OpenBLAS/include)
include_directories(/usr/local/cuda/include)

# 3. Handle Cholmod Optionality
# We look for the library explicitly. 
# (You can add more PATHS if you install it somewhere else later)
find_library(CHOLMOD_LIB NAMES cholmod PATHS /opt/SuiteSparse/lib /usr/local/lib /usr/lib)

if(CHOLMOD_LIB)
    message(STATUS ">> Cholmod FOUND at ${CHOLMOD_LIB}. High-performance solver enabled.")
    
    # Add the include path for SuiteSparse
    include_directories(/opt/SuiteSparse/include/suitesparse)
    include_directories(/usr/include/suitesparse)

    # Define the macro that the C++ code will check
    add_definitions(-DUSE_CHOLMOD)

    # List of libraries to link
    # (We assume if Cholmod is there, the others are too)
    link_directories(/opt/SuiteSparse/lib)
    set(SPARSE_LIBS cholmod amd colamd camd ccolamd suitesparseconfig openblas)

else()
    message(WARNING ">> Cholmod NOT FOUND. Defaulting to Eigen SimplicialLLT (Slower).")

    # FALLBACK: Pure Eigen
    # Eigen is header-only and compiles its own math kernels.
    set(SPARSE_LIBS "")
endif()

# 4. Create Core Library
file(GLOB SOURCE_FILES "src/*.cpp")
add_library(pffm_core ${SOURCE_FILES})

# 5. Create Executables
add_executable(solver_modeI apps/main_modeI.cpp)
target_link_libraries(solver_modeI pffm_core ${SPARSE_LIBS})

add_executable(solver_modeII apps/main_modeII.cpp)
target_link_libraries(solver_modeII pffm_core ${SPARSE_LIBS})