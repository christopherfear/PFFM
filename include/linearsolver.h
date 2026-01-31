#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <Eigen/Sparse>
#include <iostream>

/**
 * Provides a robust linear solver 'waterfall'.
 * Tries Cholmod (if available), then SimplicialLLT, then LDLT, 
 * and finally SparseLU as a fallback for non-SPD systems.
*/


// Only include Cholmod headers if CMake found the library
#ifdef USE_CHOLMOD
    #include <cholmod.h>
    #include <Eigen/CholmodSupport>
#endif

class LinearSolver {
public:
    using SpMat = Eigen::SparseMatrix<double>;
    using Vec   = Eigen::VectorXd;

    LinearSolver();
    void analyzePattern(const SpMat& K);
    bool solve(const SpMat& K, const Vec& rhs, Vec& solution);

private:
    #ifdef USE_CHOLMOD
        Eigen::CholmodSupernodalLLT<SpMat> solverCM;
    #endif

    Eigen::SimplicialLLT<SpMat>  solverLLT;
    Eigen::SimplicialLDLT<SpMat> solverLDLT;
    
    // The "Nuclear Option"
    Eigen::SparseLU<SpMat> solverLU; 
};

#endif