#include "linearsolver.h"

LinearSolver::LinearSolver() {
    // Constructor
}

void LinearSolver::analyzePattern(const SpMat& K) {
    #ifdef USE_CHOLMOD
        solverCM.analyzePattern(K);
    #endif
    solverLLT.analyzePattern(K);
    solverLDLT.analyzePattern(K);
    solverLU.analyzePattern(K); // Analyze pattern for LU too
}

bool LinearSolver::solve(const SpMat& K, const Vec& rhs, Vec& solution) {
    // 1. CHOLMOD (high performance if available)
    #ifdef USE_CHOLMOD
        solverCM.factorize(K);
        if(solverCM.info() == Eigen::Success) {
            solution = solverCM.solve(rhs);
            if(solverCM.info() == Eigen::Success) return true;
        }
    #endif

    // 2. SimplicialLLT (fastest for SPD)
    solverLLT.factorize(K);
    if(solverLLT.info() == Eigen::Success) {
        solution = solverLLT.solve(rhs);
        if(solverLLT.info() == Eigen::Success) return true;
    }

    // 3. SimplicialLDLT (handles asymmetric matrices)
    solverLDLT.factorize(K);
    if(solverLDLT.info() == Eigen::Success) {
        solution = solverLDLT.solve(rhs);
        if(solverLDLT.info() == Eigen::Success) {
            std::cout << "[Solver] Fallback to LDLT succeeded." << std::endl;
            return true;
        }
    }

    // 4. SparseLU
    // LU is robust to non-symmetry and indefiniteness
    std::cout << "[Solver] Warning: Symmetric solvers failed. Trying SparseLU..." << std::endl;
    solverLU.factorize(K);
    if(solverLU.info() == Eigen::Success) {
        solution = solverLU.solve(rhs);
        if(solverLU.info() == Eigen::Success) {
            std::cout << "[Solver] SparseLU succeeded." << std::endl;
            return true;
        }
    }

    std::cerr << "[Solver] CRITICAL FAILURE: All solvers (Cholmod, LLT, LDLT, LU) failed." << std::endl;
    return false;
}