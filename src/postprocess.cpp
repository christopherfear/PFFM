#include "postprocess.h"
#include "elements.h"
#include <iostream>
#include <fstream>
#include <iomanip>

// Linear interpolation to find exact location of the damage threshold
double CalculateCrackLength(const std::vector<std::vector<double>>& nodes, const std::vector<int>& interface_nodes, const Eigen::VectorXd& d, double L, double s_threshold, double Delta)
{
	double intact_length = L;
	int numNodes = nodes.size();
	std::vector<double> processed_y; 
	std::vector<double> processed_x;
	for(size_t i = 0; i < interface_nodes.size(); ++i) {
		int n2 = interface_nodes[i];
		double s2 = d(2*numNodes + n2);
		double x2 = nodes[n2][0];
		double y2 = nodes[n2][1];
		int row_index = -1;
		for(size_t j = 0; j < processed_y.size(); ++j) {
			if(std::fabs(y2 - processed_y[j]) < Delta/10) { row_index = j; break; }
		}
		if(row_index != -1 && x2 >= processed_x[row_index]) continue;
		if(s2 <= s_threshold && i > 0) { //i>0 safety check to stop attempted interpolation off the edge of the domain
			int n1 = interface_nodes[i-1];
			double x1 = nodes[n1][0];
			double s1 = d(2*numNodes + n1);
			double current_intact = x1 + (s_threshold - s1)/(s2 - s1)*(x2 - x1);
			if(current_intact < intact_length) intact_length = current_intact;
			if(row_index != -1) processed_x[row_index] = x2;
			else { processed_y.push_back(y2); processed_x.push_back(x2); }
		}
	}
	return L - intact_length;
}

// Apparent fracture toughness extraction
double GetGcFromdPidA(const std::vector<double>& D, const std::vector<double>& P, const std::vector<double>& a, const std::vector<double>& U, double a_min)
{
	size_t size = 0;
	for(size_t i = 0; i < a.size(); ++i)
	{
		if(a[i] > a_min)
		{
			size++;
		}
	}
	if(size < 2) return 0;
	
	Eigen::VectorXd selected_a(size); // crack lengths to include in -dPi/dA fit calculation
	Eigen::VectorXd selected_U(size); // strain energy values to include in -dPi/dA fit calculation
	Eigen::VectorXd selected_W(size); // work values to include in -dPi/dA fit calculation
	int offset = a.size() - size;
	for(size_t i = offset; i < a.size(); ++i)
	{
		selected_a(i - offset) = a[i];
		selected_U(i - offset) = U[i];
		if(i == offset && i > 0)
		{
			selected_W(i - offset) = 0.5*(P[i] + P[i - 1])*(2.*D[i] - 2.*D[i-1]); // add dW with trapezium approximation
		}
		else
		{
			selected_W(i - offset) = selected_W(i - offset - 1) + 0.5*(P[i] + P[i - 1])*(2.*D[i] - 2.*D[i-1]); // add dW with trapezium approximation
		}
	}
	
	// Perform linear regression on:
	// 1) U vs a using: (X^T*X)*[dUda, intercept] = X^T*y1
	// 2) W vs a using: (X^T*X)*[dWda, intercept] = X^T*y2
	Eigen::MatrixXd X(size, 2); // Matrix for [a, 1]
	Eigen::VectorXd y1(size); // Vector for U
	Eigen::VectorXd y2(size); // Vector for W
	for (size_t i = 0; i < size; ++i)
	{
		X(i, 0) = selected_a(i); // a
		X(i, 1) = 1.0; // Intercept term (constant 1 for linear regression)
		y1(i) = selected_U(i); // U
		y2(i) = selected_W(i); // W
	}

	Eigen::VectorXd coeffs1 = (X.transpose()*X).ldlt().solve(X.transpose()*y1);
	Eigen::VectorXd coeffs2 = (X.transpose()*X).ldlt().solve(X.transpose()*y2);
	double dUda = coeffs1(0); // dU/da
	double dWda = coeffs2(0); // dW/da

	
	double Gc = -dUda + dWda;
	return Gc;
}

double GetGcfromCC(const std::vector<double>& D, const std::vector<double>& F1, const std::vector<double>& F2, const std::vector<double>& a, double a_min, double t)
{
	size_t size = 0; // Unsigned integer for counting elements
	for(size_t i = 0; i < a.size(); ++i)
	{
		if(a[i] > a_min)
		{
			size++;
		}
	}
	if(size < 2) return 0; // insufficient points beyond a_min for linear fitting
	
	Eigen::VectorXd acubed(size);
	Eigen::VectorXd C(size);

	int offset = a.size() - size;
	for(size_t i = offset; i < a.size(); ++i)
	{
		acubed(i - offset) = a[i]*a[i]*a[i];
		C(i - offset) = D[i]/(F1[i] + F2[i]);
	}
	
	// set up the matrix X (for linear regression) and vector y
	Eigen::MatrixXd X(size, 2); // matrix for [a^3, 1]
	Eigen::VectorXd y(size); // vector for C

	for (size_t i = 0; i < size; ++i)
	{
		X(i, 0) = acubed(i); // a^3
		X(i, 1) = 1.0; // intercept term (constant 1 for linear regression)
		y(i) = C(i); // C
	}

	// perform linear regression: X*[n, A] = y
	// solve for the coefficients using the normal equation: (X^T*X)*[n, A] = X^T*y
	Eigen::VectorXd coeffs = (X.transpose()*X).ldlt().solve(X.transpose()*y);

	double n = coeffs(0); // slope based on all points up to current after a_min
	double A = coeffs(1); // intercept based on all points up to current after a_min
	
	double dC_da_init = 3.0*n*a[offset]*a[offset]; // dC/da at initation based on best fit to all points up to current after a_min
	double P_init = F1[offset] + F2[offset]; // P at initiation
	double Gc_init = 0.5*P_init*P_init*dC_da_init/t; // Gc at initiation
	return Gc_init;
}

// Field extraction for visualisation
void SaveHField(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const std::vector<double>& H, bool isQuadratic, int quadratureDegree, std::string filename, int prec, const std::vector<double>& exx, const std::vector<double>& eyy, const std::vector<double>& exy, const std::vector<double>& tr_e)
{
	Eigen::VectorXd xi_points;
	Eigen::VectorXd eta_points;
	Eigen::VectorXd weights;
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	
	int numQP = weights.size();

	// save H field to CSV
	std::ofstream outfile(filename); // overwrite
	outfile << "x,y,H,exx,eyy,exy,tr_e" << std::endl;
	outfile << std::fixed << std::setprecision(prec);
	
	for(int i = 0; i < elements.size(); ++i)
	{
		if(isQuadratic)
		{
			int n1 = elements[i][0];
			int n2 = elements[i][1];
			int n3 = elements[i][2];
			int n4 = elements[i][3];
			int n5 = elements[i][4];
			int n6 = elements[i][5];
			int n7 = elements[i][6];
			int n8 = elements[i][7];
			int n9 = elements[i][8];
			double x1 = nodes[n1][0];
			double y1 = nodes[n1][1];
			double x2 = nodes[n2][0];
			double y2 = nodes[n2][1];
			double x3 = nodes[n3][0];
			double y3 = nodes[n3][1];
			double x4 = nodes[n4][0];
			double y4 = nodes[n4][1];
			double x5 = nodes[n5][0];
			double y5 = nodes[n5][1];
			double x6 = nodes[n6][0];
			double y6 = nodes[n6][1];
			double x7 = nodes[n7][0];
			double y7 = nodes[n7][1];
			double x8 = nodes[n8][0];
			double y8 = nodes[n8][1];
			double x9 = nodes[n9][0];
			double y9 = nodes[n9][1];

			for(int p = 0; p < weights.size(); ++p) // loop integration points
			{
				double xi = xi_points(p);
				double eta = eta_points(p);
				
				double N1 = (eta*xi*(eta - 1.)*(xi - 1.))/4.;
				double N2 = -(eta*(xi*xi - 1.)*(eta - 1.))/2.;
				double N3 = (eta*xi*(eta - 1.)*(xi + 1.))/4.;
				double N4 = -(xi*(eta*eta - 1.)*(xi + 1.))/2.;
				double N5 = (eta*xi*(eta + 1.)*(xi + 1.))/4.;
				double N6 = -(eta*(xi*xi - 1.)*(eta + 1.))/2.;
				double N7 = (eta*xi*(eta + 1.)*(xi - 1.))/4.;
				double N8 = -(xi*(eta*eta - 1.)*(xi - 1.))/2.;
				double N9 = (eta*eta - 1.)*(xi*xi - 1.);
				double xp = N1*x1 + N2*x2 + N3*x3 + N4*x4 + N5*x5 + N6*x6 + N7*x7 + N8*x8 + N9*x9;
				double yp = N1*y1 + N2*y2 + N3*y3 + N4*y4 + N5*y5 + N6*y6 + N7*y7 + N8*y8 + N9*y9;

				
                int idx = i * numQP + p; // FLATTENED INDEXING

                double Hval = H[idx]; 
                double exxval = exx[idx];
                double eyyval = eyy[idx];
                double exyval = exy[idx];
                double tr_eval = tr_e[idx];
                
                outfile << xp
                    << "," << yp
                    << "," << Hval
                    << "," << exxval
                    << "," << eyyval
                    << "," << exyval
                    << "," << tr_eval
                    << std::endl;
			}
		}
		else // linear elements (4 nodes)
		{
			int n1 = elements[i][0];
			int n2 = elements[i][1];
			int n3 = elements[i][2];
			int n4 = elements[i][3];
			double x1 = nodes[n1][0];
			double y1 = nodes[n1][1];
			double x2 = nodes[n2][0];
			double y2 = nodes[n2][1];
			double x3 = nodes[n3][0];
			double y3 = nodes[n3][1];
			double x4 = nodes[n4][0];
			double y4 = nodes[n4][1];
			
			for(int p = 0; p < weights.size(); ++p) // loop integration points
			{
				double xi = xi_points(p);
				double eta = eta_points(p);
				
				double N1 = (1. - eta)*(1. - xi)/4.;
				double N2 = (1. - eta)*(1. + xi)/4.;
				double N3 = (1. + eta)*(1. + xi)/4.;
				double N4 = (1. + eta)*(1. - xi)/4.;
				double xp = N1*x1 + N2*x2 + N3*x3 + N4*x4;
				double yp = N1*y1 + N2*y2 + N3*y3 + N4*y4;

		
                int idx = i * numQP + p;	// FLATTENED INDEXING

                double Hval = H[idx]; 
                double exxval = exx[idx];
                double eyyval = eyy[idx];
                double exyval = exy[idx];
                double tr_eval = tr_e[idx];
                
                outfile << xp
                    << "," << yp
                    << "," << Hval
                    << "," << exxval
                    << "," << eyyval
                    << "," << exyval
                    << "," << tr_eval
                    << std::endl;
			}
		}
	}
	outfile.close();
}
