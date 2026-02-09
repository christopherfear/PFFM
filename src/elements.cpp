#include "elements.h"
#include <iostream>

// Numerical integration 
void Get2DQuadrature(int degree, Eigen::VectorXd& xis, Eigen::VectorXd& etas, Eigen::VectorXd& ws)
{
	if(degree == 1) // 1 point, fully integrates a 1st-degree (linear) 2D polynomial
	{
		xis.resize(1); xis << 0.0;
		etas.resize(1); etas << 0.0;
		ws.resize(1); ws << 4.0;
	}
	else if(degree == 3) // 4 points, fully integrates a 3rd-degree (cubic) 2D polynomial
	{
		double p = 1./sqrt(3);
		xis.resize(4); xis << -p, p, -p, p;
		etas.resize(4); etas << -p, -p, p, p;
		ws.resize(4); ws << 1.0, 1.0, 1.0, 1.0;
	}
	else if(degree == 5) // 9 points, fully integrates a 5th-degree (quintic) 2D polynomial
	{
		double p = sqrt(0.6);
		xis.resize(9); xis << -p, 0, p, -p, 0, p, -p, 0, p;
		etas.resize(9); etas << -p, -p, -p, 0, 0, 0, p, p, p;
		ws.resize(9); ws << 25, 40, 25, 40, 64, 40, 25, 40, 25;
		ws *= 1./81;
	}
	else if(degree == 7) // 16 points, fully integrates a 7th-degree (septic) 2D polynomial
	{
		xis.resize(16); xis << -0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053, -0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053, -0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053, -0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053;
		etas.resize(16); etas << -0.861136311594053, -0.861136311594053, -0.861136311594053, -0.861136311594053, -0.339981043584856, -0.339981043584856, -0.339981043584856, -0.339981043584856, 0.339981043584856, 0.339981043584856, 0.339981043584856, 0.339981043584856, 0.861136311594053, 0.861136311594053, 0.861136311594053, 0.861136311594053;
		ws.resize(16); ws << 0.121002993285602, 0.226851851851852, 0.226851851851852, 0.121002993285602, 0.226851851851852, 0.425293303010694, 0.425293303010694, 0.226851851851852, 0.226851851851852, 0.425293303010694, 0.425293303010694, 0.226851851851852, 0.121002993285602, 0.226851851851852, 0.226851851851852, 0.121002993285602;
	}
	else if(degree == 9) // 25 points, fully integrates a 9th-degree (nonic) 2D polynomial
	{
		xis.resize(25); xis << 0, 0, 0, 0, 0, 0.538469310105683, 0.538469310105683, 0.538469310105683, 0.538469310105683, 0.538469310105683, 0.906179845938664, 0.906179845938664, 0.906179845938664, 0.906179845938664, 0.906179845938664, -0.538469310105683, -0.538469310105683, -0.538469310105683, -0.538469310105683, -0.538469310105683, -0.906179845938664, -0.906179845938664, -0.906179845938664, -0.906179845938664, -0.906179845938664;
		etas.resize(25); etas << 0, 0.538469310105683, 0.906179845938664, -0.538469310105683, -0.906179845938664, 0, 0.538469310105683, 0.906179845938664, -0.538469310105683, -0.906179845938664, 0, 0.538469310105683, 0.906179845938664, -0.538469310105683, -0.906179845938664, 0, 0.538469310105683, 0.906179845938664, -0.538469310105683, -0.906179845938664, 0, 0.538469310105683, 0.906179845938664, -0.538469310105683, -0.906179845938664;
		ws.resize(25); ws << 0.323634567901235, 0.272286532550751, 0.134785072387521, 0.272286532550751, 0.134785072387521, 0.272286532550751, 0.229085404223991, 0.1134, 0.229085404223991, 0.1134, 0.134785072387521, 0.1134, 0.0561343488624286, 0.1134, 0.0561343488624286, 0.272286532550751, 0.229085404223991, 0.1134, 0.229085404223991, 0.1134, 0.134785072387521, 0.1134, 0.0561343488624286, 0.1134, 0.0561343488624286;
	}
	else
	{
		std::cout << "Degree not supported!" << std::endl;
	}
}


// Displacement  field elements
void LinearSolidElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double E, double nu, double t,  const Eigen::Matrix<double, 8, 1>& u, const Eigen::Matrix<double, 4, 1>& s, bool plane_stress, const Eigen::VectorXd& xi_points, const Eigen::VectorXd& eta_points, const Eigen::VectorXd& weights, 
                        Eigen::Matrix<double, 8, 8>& k)
{
	// constitutive properties matrix
	Eigen::Matrix3d D;
	Eigen::Matrix3d Dvol;
	if(plane_stress == true)
	{
		D <<
			1., nu, 0.,
			nu, 1., 0.,
			0., 0., (1. - nu)/2.;
		D *= E/(1 - nu*nu);
		Dvol <<
			1., 1., 0.,
			1., 1., 0.,
			0., 0., 0.;
		Dvol *= E/2./(1. - nu);
	}
	else
	{
		D <<
			1. - nu, nu, 0.,
			nu, 1. - nu, 0.,
			0., 0., (1. - 2.*nu)/2.;
		D *= E/(1. + nu)/(1. - 2.*nu);
		Dvol <<
			1., 1., 0.,
			1., 1., 0.,
			0., 0., 0.;
		Dvol *= E/2./(1. + nu)/(1. - 2.*nu);
	}
	Eigen::Matrix3d Ddev = D - Dvol;
	
    k.setZero();
	for(int p = 0; p < weights.size(); ++p) // loop integration points
	{
		double xi = xi_points(p);
		double eta = eta_points(p);
		double weight = weights(p);
		
		// shape functions and shape functions vector
		double N1 = (1. - eta)*(1. - xi)/4.;
		double N2 = (1. - eta)*(1. + xi)/4.;
		double N3 = (1. + eta)*(1. + xi)/4.;
		double N4 = (1. + eta)*(1. - xi)/4.;
		Eigen::Matrix<double, 4, 1> N;
		N << N1, N2, N3, N4;

		// Jacobian matrix
		double dN1_dxi = -(1. - eta)/4;
		double dN1_deta = -(1. - xi)/4;
		double dN2_dxi = (1. - eta)/4;
		double dN2_deta = -(1. + xi)/4;
		double dN3_dxi = (1. + eta)/4;
		double dN3_deta = (1. + xi)/4;
		double dN4_dxi = -(1. + eta)/4;
		double dN4_deta = (1. - xi)/4;
		double dx_dxi = dN1_dxi*x1 + dN2_dxi*x2 + dN3_dxi*x3 + dN4_dxi*x4;
		double dx_deta = dN1_deta*x1 + dN2_deta*x2 + dN3_deta*x3 + dN4_deta*x4;
		double dy_dxi = dN1_dxi*y1 + dN2_dxi*y2 + dN3_dxi*y3 + dN4_dxi*y4;
		double dy_deta = dN1_deta*y1 + dN2_deta*y2 + dN3_deta*y3 + dN4_deta*y4;
		Eigen::Matrix2d J;
		J << dx_dxi, dy_dxi, dx_deta, dy_deta;
		Eigen::Matrix2d J_inv = J.inverse();
		double J_det = J.determinant();
		
		// strain-displacement matrix
		Eigen::Matrix<double, 2, 4> dN_local;
		dN_local << dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi,
			dN1_deta, dN2_deta, dN3_deta, dN4_deta;
		Eigen::Matrix<double, 2, 4> dN_global = J_inv*dN_local;
		double dN1_dx = dN_global(0, 0);
		double dN2_dx = dN_global(0, 1);
		double dN3_dx = dN_global(0, 2);
		double dN4_dx = dN_global(0, 3);
		double dN1_dy = dN_global(1, 0);
		double dN2_dy = dN_global(1, 1);
		double dN3_dy = dN_global(1, 2);
		double dN4_dy = dN_global(1, 3);
		Eigen::Matrix<double, 3, 8> B;
		B <<
			dN1_dx, 0, dN2_dx, 0, dN3_dx, 0, dN4_dx, 0,
			0, dN1_dy, 0, dN2_dy, 0, dN3_dy, 0, dN4_dy,
			dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx, dN4_dy, dN4_dx;
		Eigen::Matrix<double, 8, 3> B_T = B.transpose();
		
		
        Eigen::Vector3d sigma = D * (B * u.segment<8>(0)); // .segment<8>(0) forces the compiler to treat 'u' as a fixed-size vector of length 8
		double sigma_vol = (sigma(0) + sigma(1))/2.; // volumetric stress (positive in tension)
		bool tensile = (sigma_vol > -std::numeric_limits<double>::epsilon()); // sigma_vol >= 0 with some tolerance for floating point error
		
		// Calculate s scalar term efficiently using fixed sizes
        // .segment<4>(0) treats 's' as fixed size 4 to match N
        double s_term = (N.transpose() * s.segment<4>(0) * s.segment<4>(0).transpose() * N)(0,0);

		//integrate the damage effected stiffness with quadrature
        if(tensile)
        {
            k += t*B_T*(s_term)*D*B*weight*J_det;
        }
        else // compressive
        {
            k += t*B_T*(s_term*Ddev + Dvol)*B*weight*J_det;
        }
    }
}

void QuadraticSolidElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, 
                           double E, double nu, double t, 
                           const Eigen::Matrix<double, 18, 1>& u, const Eigen::Matrix<double, 9, 1>& s, 
                           bool plane_stress, 
                           const Eigen::VectorXd& xi_points, const Eigen::VectorXd& eta_points, const Eigen::VectorXd& weights, 
                           Eigen::Matrix<double, 18, 18>& k)
{
	// constitutive properties matrix
	Eigen::Matrix3d D;
	Eigen::Matrix3d Dvol;
	if(plane_stress == true)
	{
		D <<
			1., nu, 0.,
			nu, 1., 0.,
			0., 0., (1. - nu)/2.;
		D *= E/(1 - nu*nu);
		Dvol <<
			1., 1., 0.,
			1., 1., 0.,
			0., 0., 0.;
		Dvol *= E/2./(1. - nu);
	}
	else
	{
		D <<
			1. - nu, nu, 0.,
			nu, 1. - nu, 0.,
			0., 0., (1. - 2.*nu)/2.;
		D *= E/(1. + nu)/(1. - 2.*nu);
		Dvol <<
			1., 1., 0.,
			1., 1., 0.,
			0., 0., 0.;
		Dvol *= E/2./(1. + nu)/(1. - 2.*nu);
	}
	Eigen::Matrix3d Ddev = D - Dvol;
	
    k.setZero();
	for(int p = 0; p < weights.size(); ++p) // loop integration points
	{
		double xi = xi_points(p);
		double eta = eta_points(p);
		double weight = weights(p);
		
		// shape functions and shape functions vector
		double N1 = (eta*xi*(eta - 1.)*(xi - 1.))/4.;
		double N2 = -(eta*(xi*xi - 1.)*(eta - 1.))/2.;
		double N3 = (eta*xi*(eta - 1.)*(xi + 1.))/4.;
		double N4 = -(xi*(eta*eta - 1.)*(xi + 1.))/2.;
		double N5 = (eta*xi*(eta + 1.)*(xi + 1.))/4.;
		double N6 = -(eta*(xi*xi - 1.)*(eta + 1.))/2.;
		double N7 = (eta*xi*(eta + 1.)*(xi - 1.))/4.;
		double N8 = -(xi*(eta*eta - 1.)*(xi - 1.))/2.;
		double N9 = (eta*eta - 1.)*(xi*xi - 1.);
		Eigen::Matrix<double, 9, 1> N;
		N << N1, N2, N3, N4, N5, N6, N7, N8, N9;

		// Jacobian matrix
		double dN1_dxi = (eta*(2.*xi - 1.)*(eta - 1.))/4.;
		double dN2_dxi = -eta*xi*(eta - 1.);
		double dN3_dxi = (eta*(2.*xi + 1.)*(eta - 1.))/4.;
		double dN4_dxi = -((eta*eta - 1.)*(2.*xi + 1.))/2.;
		double dN5_dxi = (eta*(2.*xi + 1.)*(eta + 1.))/4.;
		double dN6_dxi = -eta*xi*(eta + 1.);
		double dN7_dxi = (eta*(2.*xi - 1.)*(eta + 1.))/4.;
		double dN8_dxi = -((eta*eta - 1.)*(2.*xi - 1.))/2.;
		double dN9_dxi = 2.*xi*(eta*eta - 1.);
		double dN1_deta = (xi*(2.*eta - 1.)*(xi - 1.))/4.;
		double dN2_deta = -((2.*eta - 1.)*(xi*xi - 1.))/2.;
		double dN3_deta = (xi*(2.*eta - 1.)*(xi + 1.))/4.;
		double dN4_deta = -eta*xi*(xi + 1.);
		double dN5_deta = (xi*(2.*eta + 1.)*(xi + 1.))/4.;
		double dN6_deta = -((2.*eta + 1.)*(xi*xi - 1.))/2.;
		double dN7_deta = (xi*(2.*eta + 1.)*(xi - 1.))/4.;
		double dN8_deta = -eta*xi*(xi - 1.);
		double dN9_deta = 2.*eta*(xi*xi - 1.);
		double dx_dxi = dN1_dxi*x1 + dN2_dxi*x2 + dN3_dxi*x3 + dN4_dxi*x4 + dN5_dxi*x5 + dN6_dxi*x6 + dN7_dxi*x7 + dN8_dxi*x8 + dN9_dxi*x9;
		double dx_deta = dN1_deta*x1 + dN2_deta*x2 + dN3_deta*x3 + dN4_deta*x4 + dN5_deta*x5 + dN6_deta*x6 + dN7_deta*x7 + dN8_deta*x8 + dN9_deta*x9;
		double dy_dxi = dN1_dxi*y1 + dN2_dxi*y2 + dN3_dxi*y3 + dN4_dxi*y4 + dN5_dxi*y5 + dN6_dxi*y6 + dN7_dxi*y7 + dN8_dxi*y8 + dN9_dxi*y9;
		double dy_deta = dN1_deta*y1 + dN2_deta*y2 + dN3_deta*y3 + dN4_deta*y4 + dN5_deta*y5 + dN6_deta*y6 + dN7_deta*y7 + dN8_deta*y8 + dN9_deta*y9;
		Eigen::Matrix2d J;
		J << dx_dxi, dy_dxi, dx_deta, dy_deta;
		Eigen::Matrix2d J_inv = J.inverse();
		double J_det = J.determinant();
		
		// strain-displacement matrix
		Eigen::Matrix<double, 2, 9> dN_local;
		dN_local << dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN5_dxi, dN6_dxi, dN7_dxi, dN8_dxi, dN9_dxi, 
			dN1_deta, dN2_deta, dN3_deta, dN4_deta, dN5_deta, dN6_deta, dN7_deta, dN8_deta, dN9_deta;
		Eigen::Matrix<double, 2, 9> dN_global = J_inv*dN_local;
		double dN1_dx = dN_global(0, 0);
		double dN2_dx = dN_global(0, 1);
		double dN3_dx = dN_global(0, 2);
		double dN4_dx = dN_global(0, 3);
		double dN5_dx = dN_global(0, 4);
		double dN6_dx = dN_global(0, 5);
		double dN7_dx = dN_global(0, 6);
		double dN8_dx = dN_global(0, 7);
		double dN9_dx = dN_global(0, 8);
		double dN1_dy = dN_global(1, 0);
		double dN2_dy = dN_global(1, 1);
		double dN3_dy = dN_global(1, 2);
		double dN4_dy = dN_global(1, 3);
		double dN5_dy = dN_global(1, 4);
		double dN6_dy = dN_global(1, 5);
		double dN7_dy = dN_global(1, 6);
		double dN8_dy = dN_global(1, 7);
		double dN9_dy = dN_global(1, 8);
		Eigen::Matrix<double, 3, 18> B;
		B <<
			dN1_dx, 0, dN2_dx, 0, dN3_dx, 0, dN4_dx, 0, dN5_dx, 0, dN6_dx, 0, dN7_dx, 0, dN8_dx, 0, dN9_dx, 0,
			0, dN1_dy, 0, dN2_dy, 0, dN3_dy, 0, dN4_dy, 0, dN5_dy, 0, dN6_dy, 0, dN7_dy, 0, dN8_dy, 0, dN9_dy,
			dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx, dN4_dy, dN4_dx, dN5_dy, dN5_dx, dN6_dy, dN6_dx, dN7_dy, dN7_dx, dN8_dy, dN8_dx, dN9_dy, dN9_dx;
		Eigen::Matrix<double, 18, 3> B_T = B.transpose();
		
		Eigen::Vector3d sigma = D * (B * u.segment<18>(0)); // .segment<18>(0) forces the compiler to treat 'u' as a fixed-size vector of length 18
		double sigma_vol = (sigma(0) + sigma(1))/2.; // volumetric stress (positive in tension)
		bool tensile = (sigma_vol > -std::numeric_limits<double>::epsilon()); // sigma_vol >= 0 with some tolerance for floating point error
		
		// .segment<9>(0) for s (9 DOFs)
        double s_term = (N.transpose() * s.segment<9>(0) * s.segment<9>(0).transpose() * N)(0,0);
		
		// integrate damage-affected stiffness matrix with quadrature
        if(tensile)
        {
            k += t*B_T*(s_term)*D*B*weight*J_det;
        }
        else // compressive
        {
            k += t*B_T*(s_term*Ddev + Dvol)*B*weight*J_det;
        }
    }
}

// Phase-field fracture elements 
void LinearPFFMElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, 
                       double Gc, double ls, 
                       const Eigen::Matrix<double, 4, 1>& s,
                       const double* H,
                       const Eigen::VectorXd& xi_points, const Eigen::VectorXd& eta_points, const Eigen::VectorXd& weights, 
                       Eigen::MatrixXd& k, Eigen::VectorXd& f)
{
    if(k.rows() != 4 || k.cols() != 4) k.resize(4, 4);
    k.setZero();
    
    if(f.size() != 4) f.resize(4);
    f.setZero();
	for(int p = 0; p < weights.size(); ++p) // loop integration points
	{
		double xi = xi_points(p);
		double eta = eta_points(p);
		double weight = weights(p);
		
		// shape functions and shape functions vector
		double N1 = (1. - eta)*(1. - xi)/4.;
		double N2 = (1. - eta)*(1. + xi)/4.;
		double N3 = (1. + eta)*(1. + xi)/4.;
		double N4 = (1. + eta)*(1. - xi)/4.;
		Eigen::Matrix<double, 4, 1> N;
		N << N1, N2, N3, N4;

		// Jacobian matrix
		double dN1_dxi = -(1. - eta)/4;
		double dN1_deta = -(1. - xi)/4;
		double dN2_dxi = (1. - eta)/4;
		double dN2_deta = -(1. + xi)/4;
		double dN3_dxi = (1. + eta)/4;
		double dN3_deta = (1. + xi)/4;
		double dN4_dxi = -(1. + eta)/4;
		double dN4_deta = (1. - xi)/4;
		double dx_dxi = dN1_dxi*x1 + dN2_dxi*x2 + dN3_dxi*x3 + dN4_dxi*x4;
		double dx_deta = dN1_deta*x1 + dN2_deta*x2 + dN3_deta*x3 + dN4_deta*x4;
		double dy_dxi = dN1_dxi*y1 + dN2_dxi*y2 + dN3_dxi*y3 + dN4_dxi*y4;
		double dy_deta = dN1_deta*y1 + dN2_deta*y2 + dN3_deta*y3 + dN4_deta*y4;
		Eigen::Matrix2d J;
		J << dx_dxi, dy_dxi, dx_deta, dy_deta;
		Eigen::Matrix2d J_inv = J.inverse();
		double J_det = J.determinant();
		
		// shape function derivatives in global coordinates
		Eigen::Matrix<double, 2, 4> dN_local;
		dN_local << dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, 
			dN1_deta, dN2_deta, dN3_deta, dN4_deta;
		Eigen::Matrix<double, 2, 4> dN_global = J_inv*dN_local; // rows = dNi/dx, dNi/dy; columns = dN1/d(x,y), dN2/d(x,y), ...
		
		for(int i = 0; i < 4; ++i)
		{
			for(int j = 0; j < 4; ++j)
			{
				k(i, j) += Gc*(dN_global(0, i)*dN_global(0, j) + dN_global(1, i)*dN_global(1, j))*weight*J_det;
				k(i, j) += 2./ls*H[p]*N(i)*N(j)*weight*J_det;
				k(i, j) += Gc/(ls*ls)*N(i)*N(j)*weight*J_det;
			}
			
			f(i) += Gc/(ls*ls)*N(i)*weight*J_det;
		}
	}
}

void QuadraticPFFMElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, 
                          double Gc, double ls, 
                          const Eigen::Matrix<double, 9, 1>& s,
                          const double* H,
                          const Eigen::VectorXd& xi_points, const Eigen::VectorXd& eta_points, const Eigen::VectorXd& weights, 
                          Eigen::MatrixXd& k, Eigen::VectorXd& f)
{
	//Resize only if required
	if(k.rows() != 9 || k.cols() != 9) k.resize(9, 9);
    k.setZero();
    if(f.size() != 9) f.resize(9);
    f.setZero();

	for(int p = 0; p < weights.size(); ++p) // loop integration points
	{
		double xi = xi_points(p);
		double eta = eta_points(p);
		double weight = weights(p);
		
		// shape functions and shape functions vector
		double N1 = (eta*xi*(eta - 1.)*(xi - 1.))/4.;
		double N2 = -(eta*(xi*xi - 1.)*(eta - 1.))/2.;
		double N3 = (eta*xi*(eta - 1.)*(xi + 1.))/4.;
		double N4 = -(xi*(eta*eta - 1.)*(xi + 1.))/2.;
		double N5 = (eta*xi*(eta + 1.)*(xi + 1.))/4.;
		double N6 = -(eta*(xi*xi - 1.)*(eta + 1.))/2.;
		double N7 = (eta*xi*(eta + 1.)*(xi - 1.))/4.;
		double N8 = -(xi*(eta*eta - 1.)*(xi - 1.))/2.;
		double N9 = (eta*eta - 1.)*(xi*xi - 1.);
		Eigen::Matrix<double, 9, 1> N;
		N << N1, N2, N3, N4, N5, N6, N7, N8, N9;

		// Jacobian matrix
		double dN1_dxi = (eta*(2.*xi - 1.)*(eta - 1.))/4.;
		double dN2_dxi = -eta*xi*(eta - 1.);
		double dN3_dxi = (eta*(2.*xi + 1.)*(eta - 1.))/4.;
		double dN4_dxi = -((eta*eta - 1.)*(2.*xi + 1.))/2.;
		double dN5_dxi = (eta*(2.*xi + 1.)*(eta + 1.))/4.;
		double dN6_dxi = -eta*xi*(eta + 1.);
		double dN7_dxi = (eta*(2.*xi - 1.)*(eta + 1.))/4.;
		double dN8_dxi = -((eta*eta - 1.)*(2.*xi - 1.))/2.;
		double dN9_dxi = 2.*xi*(eta*eta - 1.);
		double dN1_deta = (xi*(2.*eta - 1.)*(xi - 1.))/4.;
		double dN2_deta = -((2.*eta - 1.)*(xi*xi - 1.))/2.;
		double dN3_deta = (xi*(2.*eta - 1.)*(xi + 1.))/4.;
		double dN4_deta = -eta*xi*(xi + 1.);
		double dN5_deta = (xi*(2.*eta + 1.)*(xi + 1.))/4.;
		double dN6_deta = -((2.*eta + 1.)*(xi*xi - 1.))/2.;
		double dN7_deta = (xi*(2.*eta + 1.)*(xi - 1.))/4.;
		double dN8_deta = -eta*xi*(xi - 1.);
		double dN9_deta = 2.*eta*(xi*xi - 1.);
		double dx_dxi = dN1_dxi*x1 + dN2_dxi*x2 + dN3_dxi*x3 + dN4_dxi*x4 + dN5_dxi*x5 + dN6_dxi*x6 + dN7_dxi*x7 + dN8_dxi*x8 + dN9_dxi*x9;
		double dx_deta = dN1_deta*x1 + dN2_deta*x2 + dN3_deta*x3 + dN4_deta*x4 + dN5_deta*x5 + dN6_deta*x6 + dN7_deta*x7 + dN8_deta*x8 + dN9_deta*x9;
		double dy_dxi = dN1_dxi*y1 + dN2_dxi*y2 + dN3_dxi*y3 + dN4_dxi*y4 + dN5_dxi*y5 + dN6_dxi*y6 + dN7_dxi*y7 + dN8_dxi*y8 + dN9_dxi*y9;
		double dy_deta = dN1_deta*y1 + dN2_deta*y2 + dN3_deta*y3 + dN4_deta*y4 + dN5_deta*y5 + dN6_deta*y6 + dN7_deta*y7 + dN8_deta*y8 + dN9_deta*y9;
		Eigen::Matrix2d J;
		J << dx_dxi, dy_dxi, dx_deta, dy_deta;
		Eigen::Matrix2d J_inv = J.inverse();
		double J_det = J.determinant();
		
		// shape function derivatives in global coordinates
		Eigen::Matrix<double, 2, 9> dN_local;
		dN_local << dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN5_dxi, dN6_dxi, dN7_dxi, dN8_dxi, dN9_dxi, 
			dN1_deta, dN2_deta, dN3_deta, dN4_deta, dN5_deta, dN6_deta, dN7_deta, dN8_deta, dN9_deta;
		Eigen::Matrix<double, 2, 9> dN_global = J_inv*dN_local;
		
		for(int i = 0; i < 9; ++i)
		{
			for(int j = 0; j < 9; ++j)
			{
				k(i, j) += Gc*(dN_global(0, i)*dN_global(0, j) + dN_global(1, i)*dN_global(1, j))*weight*J_det;
				k(i, j) += 2./ls*H[p]*N(i)*N(j)*weight*J_det;
				k(i, j) += Gc/(ls*ls)*N(i)*N(j)*weight*J_det;
			}
			
			f(i) += Gc/(ls*ls)*N(i)*weight*J_det;
		}
	}
}

// Energy and state utilities
void GetLinearElementTensileElasticStrainEnergy(
    double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, 
    double E, double nu, const Eigen::VectorXd& u, bool plane_stress, 
    const Eigen::VectorXd& xi_points, const Eigen::VectorXd& eta_points, const Eigen::VectorXd& weights,
    double* H_new, 
	double* exx_i, double* eyy_i, double* exy_i, double* tr_e_i)
{
    // Constitutive Matrix D
    Eigen::MatrixXd D(3, 3);
    double K;
    if(plane_stress) {
        D << 1., nu, 0., nu, 1., 0., 0., 0., (1. - nu)/2.;
        D *= E/(1. - nu*nu);
        K = E/2./(1. - nu);
    } else {
        D << 1. - nu, nu, 0., nu, 1. - nu, 0., 0., 0., (1. - 2.*nu)/2.;
        D *= E/(1. + nu)/(1. - 2.*nu);
        K = E/2./(1. + nu)/(1. - 2.*nu);
    }

    // Ensure vectors are the right size (Cheap if size hasn't changed)
    int num_points = weights.size();
	for(int p = 0; p < weights.size(); ++p) // loop integration points
	{
		double xi = xi_points(p);
		double eta = eta_points(p);
		double weight = weights(p);
		
		// shape functions and shape functions vector
		double N1 = (1. - eta)*(1. - xi)/4.;
		double N2 = (1. - eta)*(1. + xi)/4.;
		double N3 = (1. + eta)*(1. + xi)/4.;
		double N4 = (1. + eta)*(1. - xi)/4.;
		Eigen::VectorXd N(4);
		N << N1, N2, N3, N4;

		// Jacobian matrix
		double dN1_dxi = -(1. - eta)/4;
		double dN1_deta = -(1. - xi)/4;
		double dN2_dxi = (1. - eta)/4;
		double dN2_deta = -(1. + xi)/4;
		double dN3_dxi = (1. + eta)/4;
		double dN3_deta = (1. + xi)/4;
		double dN4_dxi = -(1. + eta)/4;
		double dN4_deta = (1. - xi)/4;
		double dx_dxi = dN1_dxi*x1 + dN2_dxi*x2 + dN3_dxi*x3 + dN4_dxi*x4;
		double dx_deta = dN1_deta*x1 + dN2_deta*x2 + dN3_deta*x3 + dN4_deta*x4;
		double dy_dxi = dN1_dxi*y1 + dN2_dxi*y2 + dN3_dxi*y3 + dN4_dxi*y4;
		double dy_deta = dN1_deta*y1 + dN2_deta*y2 + dN3_deta*y3 + dN4_deta*y4;
		Eigen::Matrix2d J;
		J << dx_dxi, dy_dxi, dx_deta, dy_deta;
		Eigen::Matrix2d J_inv = J.inverse();
		
		// strain-displacement matrix
		Eigen::MatrixXd dN_local(2, 4);
		dN_local << dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi,
			dN1_deta, dN2_deta, dN3_deta, dN4_deta;
		Eigen::MatrixXd dN_global = J_inv*dN_local;
		double dN1_dx = dN_global(0, 0);
		double dN2_dx = dN_global(0, 1);
		double dN3_dx = dN_global(0, 2);
		double dN4_dx = dN_global(0, 3);
		double dN1_dy = dN_global(1, 0);
		double dN2_dy = dN_global(1, 1);
		double dN3_dy = dN_global(1, 2);
		double dN4_dy = dN_global(1, 3);
		Eigen::MatrixXd B(3, 8);
		B <<
			dN1_dx, 0, dN2_dx, 0, dN3_dx, 0, dN4_dx, 0,
			0, dN1_dy, 0, dN2_dy, 0, dN3_dy, 0, dN4_dy,
			dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx, dN4_dy, dN4_dx;
		Eigen::MatrixXd B_T = B.transpose();
		
		Eigen::Vector3d epsilon = B*u;
        Eigen::Vector3d sigma = D*epsilon;
        double epsilon_vol = epsilon(0) + epsilon(1);
        double sigma_vol = (sigma(0) + sigma(1))/2.;
        bool tensile = (sigma_vol > -1e-16); 

        double U = 0.5*sigma.transpose()*epsilon;
        double U_vol = 0.5*K*epsilon_vol*epsilon_vol;
        double U_dev = U - U_vol;

        // Write to output vectors
        H_new[p] = tensile ? U : U_dev;
        exx_i[p] = epsilon(0);
        eyy_i[p] = epsilon(1);
        exy_i[p] = epsilon(2);
        tr_e_i[p] = epsilon_vol;
    }
}

void GetQuadraticElementTensileElasticStrainEnergy(
    double x1, double y1, double x2, double y2, double x3, double y3, 
    double x4, double y4, double x5, double y5, double x6, double y6, 
    double x7, double y7, double x8, double y8, double x9, double y9, 
    double E, double nu, const Eigen::VectorXd& u, bool plane_stress, 
    const Eigen::VectorXd& xi_points, const Eigen::VectorXd& eta_points, const Eigen::VectorXd& weights,
	double* H_out, 
	double* exx_i, double* eyy_i, 
	double* exy_i, double* tr_e_i)
{
    Eigen::MatrixXd D(3, 3);
    double K;
    if(plane_stress) {
        D << 1., nu, 0., nu, 1., 0., 0., 0., (1. - nu)/2.;
        D *= E/(1. - nu*nu);
        K = E/2./(1. - nu);
    } else {
        D << 1. - nu, nu, 0., nu, 1. - nu, 0., 0., 0., (1. - 2.*nu)/2.;
        D *= E/(1. + nu)/(1. - 2.*nu);
        K = E/2./(1. + nu)/(1. - 2.*nu);
    }

    int nPoints = weights.size();

    for(int p = 0; p < nPoints; ++p) 
    {
        double xi = xi_points(p);
        double eta = eta_points(p);
		
		// shape functions and shape functions vector
		double N1 = (eta*xi*(eta - 1.)*(xi - 1.))/4.;
		double N2 = -(eta*(xi*xi - 1.)*(eta - 1.))/2.;
		double N3 = (eta*xi*(eta - 1.)*(xi + 1.))/4.;
		double N4 = -(xi*(eta*eta - 1.)*(xi + 1.))/2.;
		double N5 = (eta*xi*(eta + 1.)*(xi + 1.))/4.;
		double N6 = -(eta*(xi*xi - 1.)*(eta + 1.))/2.;
		double N7 = (eta*xi*(eta + 1.)*(xi - 1.))/4.;
		double N8 = -(xi*(eta*eta - 1.)*(xi - 1.))/2.;
		double N9 = (eta*eta - 1.)*(xi*xi - 1.);
		Eigen::VectorXd N(9);
		N << N1, N2, N3, N4, N5, N6, N7, N8, N9;

		// Jacobian matrix
		double dN1_dxi = (eta*(2.*xi - 1.)*(eta - 1.))/4.;
		double dN2_dxi = -eta*xi*(eta - 1.);
		double dN3_dxi = (eta*(2.*xi + 1.)*(eta - 1.))/4.;
		double dN4_dxi = -((eta*eta - 1.)*(2.*xi + 1.))/2.;
		double dN5_dxi = (eta*(2.*xi + 1.)*(eta + 1.))/4.;
		double dN6_dxi = -eta*xi*(eta + 1.);
		double dN7_dxi = (eta*(2.*xi - 1.)*(eta + 1.))/4.;
		double dN8_dxi = -((eta*eta - 1.)*(2.*xi - 1.))/2.;
		double dN9_dxi = 2.*xi*(eta*eta - 1.);
		double dN1_deta = (xi*(2.*eta - 1.)*(xi - 1.))/4.;
		double dN2_deta = -((2.*eta - 1.)*(xi*xi - 1.))/2.;
		double dN3_deta = (xi*(2.*eta - 1.)*(xi + 1.))/4.;
		double dN4_deta = -eta*xi*(xi + 1.);
		double dN5_deta = (xi*(2.*eta + 1.)*(xi + 1.))/4.;
		double dN6_deta = -((2.*eta + 1.)*(xi*xi - 1.))/2.;
		double dN7_deta = (xi*(2.*eta + 1.)*(xi - 1.))/4.;
		double dN8_deta = -eta*xi*(xi - 1.);
		double dN9_deta = 2.*eta*(xi*xi - 1.);
		double dx_dxi = dN1_dxi*x1 + dN2_dxi*x2 + dN3_dxi*x3 + dN4_dxi*x4 + dN5_dxi*x5 + dN6_dxi*x6 + dN7_dxi*x7 + dN8_dxi*x8 + dN9_dxi*x9;
		double dx_deta = dN1_deta*x1 + dN2_deta*x2 + dN3_deta*x3 + dN4_deta*x4 + dN5_deta*x5 + dN6_deta*x6 + dN7_deta*x7 + dN8_deta*x8 + dN9_deta*x9;
		double dy_dxi = dN1_dxi*y1 + dN2_dxi*y2 + dN3_dxi*y3 + dN4_dxi*y4 + dN5_dxi*y5 + dN6_dxi*y6 + dN7_dxi*y7 + dN8_dxi*y8 + dN9_dxi*y9;
		double dy_deta = dN1_deta*y1 + dN2_deta*y2 + dN3_deta*y3 + dN4_deta*y4 + dN5_deta*y5 + dN6_deta*y6 + dN7_deta*y7 + dN8_deta*y8 + dN9_deta*y9;
		Eigen::Matrix2d J;
		J << dx_dxi, dy_dxi, dx_deta, dy_deta;
		Eigen::Matrix2d J_inv = J.inverse();
		
		// strain-displacement matrix
		Eigen::MatrixXd dN_local(2, 9);
		dN_local << dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN5_dxi, dN6_dxi, dN7_dxi, dN8_dxi, dN9_dxi, 
			dN1_deta, dN2_deta, dN3_deta, dN4_deta, dN5_deta, dN6_deta, dN7_deta, dN8_deta, dN9_deta;
		Eigen::MatrixXd dN_global = J_inv*dN_local;
		double dN1_dx = dN_global(0, 0);
		double dN2_dx = dN_global(0, 1);
		double dN3_dx = dN_global(0, 2);
		double dN4_dx = dN_global(0, 3);
		double dN5_dx = dN_global(0, 4);
		double dN6_dx = dN_global(0, 5);
		double dN7_dx = dN_global(0, 6);
		double dN8_dx = dN_global(0, 7);
		double dN9_dx = dN_global(0, 8);
		double dN1_dy = dN_global(1, 0);
		double dN2_dy = dN_global(1, 1);
		double dN3_dy = dN_global(1, 2);
		double dN4_dy = dN_global(1, 3);
		double dN5_dy = dN_global(1, 4);
		double dN6_dy = dN_global(1, 5);
		double dN7_dy = dN_global(1, 6);
		double dN8_dy = dN_global(1, 7);
		double dN9_dy = dN_global(1, 8);
		Eigen::MatrixXd B(3, 18);
		B <<
			dN1_dx, 0, dN2_dx, 0, dN3_dx, 0, dN4_dx, 0, dN5_dx, 0, dN6_dx, 0, dN7_dx, 0, dN8_dx, 0, dN9_dx, 0,
			0, dN1_dy, 0, dN2_dy, 0, dN3_dy, 0, dN4_dy, 0, dN5_dy, 0, dN6_dy, 0, dN7_dy, 0, dN8_dy, 0, dN9_dy,
			dN1_dy, dN1_dx, dN2_dy, dN2_dx, dN3_dy, dN3_dx, dN4_dy, dN4_dx, dN5_dy, dN5_dx, dN6_dy, dN6_dx, dN7_dy, dN7_dx, dN8_dy, dN8_dx, dN9_dy, dN9_dx;
		Eigen::MatrixXd B_T = B.transpose();
		
		Eigen::Vector3d epsilon = B*u;
        Eigen::Vector3d sigma = D*epsilon;
        
        double epsilon_vol = epsilon(0) + epsilon(1);
        double sigma_vol = (sigma(0) + sigma(1))/2.;
        bool tensile = (sigma_vol > -1e-16);

        double U = 0.5*sigma.transpose()*epsilon;
        double U_vol = 0.5*K*epsilon_vol*epsilon_vol;
        double U_dev = U - U_vol;

        // Write directly to output buffer
        H_out[p] = tensile ? U : U_dev;
        
        exx_i[p] = epsilon(0);
        eyy_i[p] = epsilon(1);
        exy_i[p] = epsilon(2);
        tr_e_i[p] = epsilon_vol;
    }
}