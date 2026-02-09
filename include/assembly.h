#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <Eigen/Sparse>
#include <vector>
#include <unordered_map>
#include <iostream>

/**
 * Manages the transition from local element physics to the global system.
 * Handles monolithic assembly, history field irreversibility, and 
 * Dirichlet boundary condition enforcement.
*/

// H state variable update
void UpdateH(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const Eigen::VectorXd& u, std::vector<double>& H, double E_int, double E_bulk, double nu_int, double nu_bulk, double h2, double hi, double ls, bool isPlaneStress, bool isQuadratic, int quadratureDegree, std::vector<double>& exx, std::vector<double>& eyy, std::vector<double>& exy, std::vector<double>& tr_e);

// main assembly (K is the global stifness matrix of size (3*numNodes, 3*numNodes), ordered u_x,u_y,s)
void AssembleAll(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const Eigen::VectorXd& u, const Eigen::VectorXd& s, const std::vector<double>& H, double E_int, double E_bulk, double nu_int, double nu_bulk, double Gc_eff, double Gc_bulk, double t, double h2, double hi, double ls, Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, std::vector<Eigen::Triplet<double>>& triplets, bool isPlaneStress, bool isQuadratic, int quadratureDegree);

// Boundary Conditions
void ApplyDirichletBC(const Eigen::SparseMatrix<double> &K, const Eigen::VectorXd &F, const std::vector<bool> &isConstrained, const std::unordered_map<int, double*> &constrainedMap, const std::vector<int> &freeMap, Eigen::SparseMatrix<double> &KR, Eigen::VectorXd &FR);

// Helpers
void AssembleLinearSolidMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m);
void AssembleLinearPFFMMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m);
void AssembleLinearPFFMVector(Eigen::VectorXd& F, const Eigen::VectorXd& f, int i, int j, int l, int m);
void AssembleQuadraticSolidMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m, int n, int o, int p, int q, int r);
void AssembleQuadraticPFFMMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m, int n, int o, int p, int q, int r);
void AssembleQuadraticPFFMVector(Eigen::VectorXd& F, const Eigen::VectorXd& f, int i, int j, int l, int m, int n, int o, int p, int q, int r);

#endif