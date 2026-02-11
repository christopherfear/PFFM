cmake_minimum_required(VERSION 3.10)
project(PFFM_Solver)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE) #true for run false for debug

# 1. Settings
set(CMAKE_CXX_STANDARD 14)
# Optimisation flags

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall") #debugging
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fno-math-errno -DNDEBUG -march=native") #optimised running

# 2. Include Paths
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories($ENV{HOME}/PFFM/eigen-3.4.0) 
include_directories(/opt/OpenBLAS/include)
include_directories(/usr/local/cuda/include)

# 3. Handle Cholmod
find_library(CHOLMOD_LIB NAMES cholmod PATHS /opt/SuiteSparse/lib /usr/local/lib /usr/lib)
if(CHOLMOD_LIB)
    message(STATUS ">> Cholmod FOUND. High-performance solver enabled.")
    include_directories(/opt/SuiteSparse/include/suitesparse)
    include_directories(/usr/include/suitesparse)
    add_definitions(-DUSE_CHOLMOD)
    link_directories(/opt/SuiteSparse/lib)
    set(SPARSE_LIBS cholmod amd colamd camd ccolamd suitesparseconfig openblas)
else()
    message(WARNING ">> Cholmod NOT FOUND. Using Eigen fallback.")
    set(SPARSE_LIBS "")
endif()

# 4. Create Core Library
file(GLOB SOURCE_FILES "src/*.cpp")
add_library(pffm_core ${SOURCE_FILES})

# 5. Create Executables
add_executable(pffm_solver apps/main.cpp) 
target_link_libraries(pffm_solver pffm_core ${SPARSE_LIBS})
