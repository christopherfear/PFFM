#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

/**
 * Phase-field fracture element library
 * Includes implementations for both Linear (4-node) and Quadratic (9-node) 
 * solid a
nd PFFM elements using spectral decomposition.
*/

void Get2DQuadrature(int degree, Eigen::VectorXd& xis, Eigen::VectorXd& etas, Eigen::VectorXd& ws);

// solid mechanics (includes damage effect)
void LinearSolidElement(
    double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, 
    double E, double nu, double t, 
    const Eigen::Matrix<double, 8, 1>& u, const Eigen::Matrix<double, 4, 1>& s, // <--- FIXED SIZE
    bool plane_stress, 
    const Eigen::VectorXd& xi, const Eigen::VectorXd& eta, const Eigen::VectorXd& weights, 
    Eigen::Matrix<double, 8, 8>& k);

void QuadraticSolidElement(
    double x1, double y1, double x2, double y2, double x3, double y3, 
    double x4, double y4, double x5, double y5, double x6, double y6, 
    double x7, double y7, double x8, double y8, double x9, double y9, 
    double E, double nu, double t, 
    const Eigen::Matrix<double, 18, 1>& u, const Eigen::Matrix<double, 9, 1>& s,
    bool plane_stress, 
    const Eigen::VectorXd& xi, const Eigen::VectorXd& eta, const Eigen::VectorXd& weights, 
    Eigen::Matrix<double, 18, 18>& k);

// PFFM
void LinearPFFMElement(
    double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, 
    double Gc, double ls, 
    const Eigen::Matrix<double, 4, 1>& s, const double* H,
    const Eigen::VectorXd& xi, const Eigen::VectorXd& eta, const Eigen::VectorXd& weights,
    Eigen::MatrixXd& k, Eigen::VectorXd& f);

void QuadraticPFFMElement(
    double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, 
    double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, 
    double Gc, double ls, 
    const Eigen::Matrix<double, 9, 1>& s, const double* H,
    const Eigen::VectorXd& xi, const Eigen::VectorXd& eta, const Eigen::VectorXd& weights,
    Eigen::MatrixXd& k, Eigen::VectorXd& f);

//Elastic strain energy 
void GetLinearElementTensileElasticStrainEnergy(
    double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, 
    double E, double nu, const Eigen::VectorXd& u, bool plane_stress, 
    const Eigen::VectorXd& xi_points, const Eigen::VectorXd& eta_points, const Eigen::VectorXd& weights,
    double* H_new, 
    double* exx_i, double* eyy_i, double* exy_i, double* tr_e_i
);
void GetQuadraticElementTensileElasticStrainEnergy(
    double x1, double y1, double x2, double y2, double x3, double y3, 
    double x4, double y4, double x5, double y5, double x6, double y6, 
    double x7, double y7, double x8, double y8, double x9, double y9, 
    double E, double nu, const Eigen::VectorXd& u, bool plane_stress, 
    // INPUTS: Pre-calculated quadrature (const ref = fast)
    const Eigen::VectorXd& xi, const Eigen::VectorXd& eta, const Eigen::VectorXd& weights,
    // OUTPUTS: Pre-allocated buffers (pass by ref = no malloc)
    double* H_out, 
    double* exx_i, double* eyy_i, 
    double* exy_i, double* tr_e_i
);

// NOTE
// Node ordering for Quadratic Elements:
//  4---7---3
//  |   |   |
//  8---9---6
//  |   |   |
//  1---5---2


#endif