#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <vector>
#include <string>
#include <Eigen/Dense>

/**
 * Facilitates the retrieval of apparent fracture toughness for validation 
 * of effective properties and tracks crack propagation trends.
*/

// crack length 
double CalculateCrackLength(const std::vector<std::vector<double>>& nodes, const std::vector<int>& interface_nodes, const Eigen::VectorXd& d, double L, double s_threshold, double Delta);

// Mode I Method
double GetGcFromdPidA(const std::vector<double>& D, const std::vector<double>& P, const std::vector<double>& a, const std::vector<double>& U, double a_min);

// Mode II Method
double GetGcfromCC(const std::vector<double>& D1, const std::vector<double>& F1, const std::vector<double>& F2, const std::vector<double>& a, double a_min, double t);

// Output
void SaveHField(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const std::vector<std::vector<double>>& H, bool isQuadratic, int quadratureDegree, std::string filename, int prec, const std::vector<std::vector<double>>& exx, const std::vector<std::vector<double>>& eyy, const std::vector<std::vector<double>>& exy, const std::vector<std::vector<double>>& tr_e);

#endif