#include "assembly.h"
#include "elements.h" 

// Helper Functions

void AssembleLinearSolidMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m)
{
	int nodes[4] = {i, j, l, m};
    
	for (int row = 0; row < 4; ++row)
	{
		for (int col = 0; col < 4; ++col)
		{
			int row_idx = 2*nodes[row];
			int col_idx = 2*nodes[col];

			triplets.emplace_back(row_idx, col_idx, k(2*row, 2*col));
			triplets.emplace_back(row_idx, col_idx + 1, k(2*row, 2*col + 1));
			triplets.emplace_back(row_idx + 1, col_idx, k(2*row + 1, 2*col));
			triplets.emplace_back(row_idx + 1, col_idx + 1, k(2*row + 1, 2*col + 1));
		}
	}
}

void AssembleLinearPFFMMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m)
{
	int nodes[4] = {i, j, l, m};
    
	for (int row = 0; row < 4; ++row)
	{
		for (int col = 0; col < 4; ++col)
		{
			int row_idx = nodes[row];
			int col_idx = nodes[col];

			triplets.emplace_back(row_idx, col_idx, k(row, col));
		}
	}
}

void AssembleLinearPFFMVector(Eigen::VectorXd& F, const Eigen::VectorXd& f, int i, int j, int l, int m)
{
	int nodes[4] = {i, j, l, m};
    
	for (int row = 0; row < 4; ++row)
	{
		int row_idx = nodes[row];
		F(row_idx) += f(row);
	}
}

void AssembleQuadraticSolidMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m, int n, int o, int p, int q, int r)
{
	int nodes[9] = {i, j, l, m, n, o, p, q, r};
    
	for (int row = 0; row < 9; ++row)
	{
		for (int col = 0; col < 9; ++col)
		{
			int row_idx = 2*nodes[row];
			int col_idx = 2*nodes[col];

			triplets.emplace_back(row_idx, col_idx, k(2*row, 2*col));
			triplets.emplace_back(row_idx, col_idx + 1, k(2*row, 2*col + 1));
			triplets.emplace_back(row_idx + 1, col_idx, k(2*row + 1, 2*col));
			triplets.emplace_back(row_idx + 1, col_idx + 1, k(2*row + 1, 2*col + 1));
		}
	}
}

void AssembleQuadraticPFFMMatrix(std::vector<Eigen::Triplet<double>>& triplets, const Eigen::MatrixXd& k, int i, int j, int l, int m, int n, int o, int p, int q, int r)
{
	int nodes[9] = {i, j, l, m, n, o, p, q, r};
    
	for (int row = 0; row < 9; ++row)
	{
		for (int col = 0; col < 9; ++col)
		{
			int row_idx = nodes[row];
			int col_idx = nodes[col];

			triplets.emplace_back(row_idx, col_idx, k(row, col));
		}
	}
}

void AssembleQuadraticPFFMVector(Eigen::VectorXd& F, const Eigen::VectorXd& f, int i, int j, int l, int m, int n, int o, int p, int q, int r)
{
	int nodes[9] = {i, j, l, m, n, o, p, q, r};
    
	for (int row = 0; row < 9; ++row)
	{
		int row_idx = nodes[row];
		F(row_idx) += f(row);
	}
}


//Global assmebly
void AssembleAll(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const Eigen::VectorXd& u, const Eigen::VectorXd& s, const std::vector<std::vector<double>>& H, double E_int, double E_bulk, double nu_int, double nu_bulk, double Gc_eff, double Gc_bulk, double t, double h2, double hi, double ls, Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F, bool isPlaneStress, bool isQuadratic, int quadratureDegree)
{
	// use triplets to efficiently construct sparse matrices
	std::vector<Eigen::Triplet<double>> triplets;
	if(isQuadratic)
	{
		triplets.reserve(324*elements.size() + 81*elements.size()); // (9*9*4 + 9*9)*elements.size()
	}
	else
	{
		triplets.reserve(64*elements.size() + 16*elements.size()); // (4*4*4 + 4*4)*elements.size()
	}
	
	F.setZero();
	Eigen::VectorXd xi, eta, weights;
    Get2DQuadrature(quadratureDegree, xi, eta, weights);

    // 2. STACK BUFFERS (Fastest Filling)
    Eigen::Matrix<double, 8, 1> u_el_lin;
    Eigen::Matrix<double, 4, 1> s_el_lin;
    Eigen::Matrix<double, 8, 8> k_solid_lin;
    
    // PFFM matrices (Dynamic is fine for these small 4x4s, or hoist them as MatrixXd)
    Eigen::MatrixXd k_pffm_lin(4, 4); 
    Eigen::VectorXd f_pffm_lin(4);

    Eigen::Matrix<double, 18, 1> u_el_quad;
    Eigen::Matrix<double, 9, 1> s_el_quad;
    Eigen::Matrix<double, 18, 18> k_solid_quad;
    
    Eigen::MatrixXd k_pffm_quad(9, 9);
    Eigen::VectorXd f_pffm_quad(9);
	
	std::cout << "Assembling element " << std::flush;
	for(int i = 0; i < elements.size(); ++i)
    {
		const std::vector<double>& H_el = H[i];
        if(isQuadratic)
        {
            // Node indices
            int n1 = elements[i][0]; int n2 = elements[i][1]; int n3 = elements[i][2];
            int n4 = elements[i][3]; int n5 = elements[i][4]; int n6 = elements[i][5];
            int n7 = elements[i][6]; int n8 = elements[i][7]; int n9 = elements[i][8];
            
            // Coordinates
            double x1 = nodes[n1][0]; double y1 = nodes[n1][1];
            double x2 = nodes[n2][0]; double y2 = nodes[n2][1];
            double x3 = nodes[n3][0]; double y3 = nodes[n3][1];
            double x4 = nodes[n4][0]; double y4 = nodes[n4][1];
            double x5 = nodes[n5][0]; double y5 = nodes[n5][1];
            double x6 = nodes[n6][0]; double y6 = nodes[n6][1];
            double x7 = nodes[n7][0]; double y7 = nodes[n7][1];
            double x8 = nodes[n8][0]; double y8 = nodes[n8][1];
            double x9 = nodes[n9][0]; double y9 = nodes[n9][1];
            
            double Gc = (y9 > h2 && y9 < h2 + hi) ? Gc_eff : Gc_bulk;
            double E = (y9 > h2 && y9 < h2 + hi) ? E_int : E_bulk;
            double nu = (y9 > h2 && y9 < h2 + hi) ? nu_int : nu_bulk;
            
            // Fill pre-allocated vectors (Fast copy)
            s_el_quad << s(n1), s(n2), s(n3), s(n4), s(n5), s(n6), s(n7), s(n8), s(n9);
            u_el_quad << u(2*n1), u(2*n1 + 1), u(2*n2), u(2*n2 + 1), u(2*n3), u(2*n3 + 1), 
                         u(2*n4), u(2*n4 + 1), u(2*n5), u(2*n5 + 1), u(2*n6), u(2*n6 + 1), 
                         u(2*n7), u(2*n7 + 1), u(2*n8), u(2*n8 + 1), u(2*n9), u(2*n9 + 1);
            
            // 1. SOLID: Pass buffer by reference
           	QuadraticSolidElement(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8, x9, y9, 
                                  E, nu, t, u_el_quad, s_el_quad, isPlaneStress, xi, eta, weights, 
                                  k_solid_quad);

            AssembleQuadraticSolidMatrix(triplets, k_solid_quad, n1, n2, n3, n4, n5, n6, n7, n8, n9);
            
            // 2. PFFM: Pass buffers by reference
            QuadraticPFFMElement(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8, x9, y9, 
                                 Gc, ls, s_el_quad, H_el, xi, eta, weights, 
                                 k_pffm_quad, f_pffm_quad);

            int offset = 2*nodes.size();
            AssembleQuadraticPFFMMatrix(triplets, k_pffm_quad, offset + n1, offset + n2, offset + n3, offset + n4, offset + n5, offset + n6, offset + n7, offset + n8, offset + n9);
            AssembleQuadraticPFFMVector(F, f_pffm_quad, offset + n1, offset + n2, offset + n3, offset + n4, offset + n5, offset + n6, offset + n7, offset + n8, offset + n9);
        }
        else // Linear Elements
        {
            int n1 = elements[i][0]; int n2 = elements[i][1]; 
            int n3 = elements[i][2]; int n4 = elements[i][3];
            
            double x1 = nodes[n1][0]; double y1 = nodes[n1][1];
            double x2 = nodes[n2][0]; double y2 = nodes[n2][1];
            double x3 = nodes[n3][0]; double y3 = nodes[n3][1];
            double x4 = nodes[n4][0]; double y4 = nodes[n4][1];
            
            double yc = (y1 + y2 + y3 + y4)/4;
            double Gc = (yc > h2 && yc < h2 + hi) ? Gc_eff : Gc_bulk;
            double E = (yc > h2 && yc < h2 + hi) ? E_int : E_bulk;
            double nu = (yc > h2 && yc < h2 + hi) ? nu_int : nu_bulk;
            
            // Fill pre-allocated vectors
            s_el_lin << s(n1), s(n2), s(n3), s(n4);
            u_el_lin << u(2*n1), u(2*n1 + 1), u(2*n2), u(2*n2 + 1), 
                        u(2*n3), u(2*n3 + 1), u(2*n4), u(2*n4 + 1);
            
			// 1. SOLID
            LinearSolidElement(x1, y1, x2, y2, x3, y3, x4, y4, 
                               E, nu, t, u_el_lin, s_el_lin, isPlaneStress, xi, eta, weights, 
                               k_solid_lin);

            AssembleLinearSolidMatrix(triplets, k_solid_lin, n1, n2, n3, n4);
            
            // 2. PFFM
            LinearPFFMElement(x1, y1, x2, y2, x3, y3, x4, y4, 
                              Gc, ls, s_el_lin, H_el, xi, eta, weights, 
                              k_pffm_lin, f_pffm_lin);

            int offset = 2*nodes.size();
            AssembleLinearPFFMMatrix(triplets, k_pffm_lin, offset + n1, offset + n2, offset + n3, offset + n4);
            AssembleLinearPFFMVector(F, f_pffm_lin, offset + n1, offset + n2, offset + n3, offset + n4);
        }
        
        if((i + 1) % 10000 == 0) std::cout << i + 1 << ", " << std::flush;
        else if(i == elements.size() - 1) std::cout << i + 1 << std::endl;
    }
    
    K.setFromTriplets(triplets.begin(), triplets.end());
}

void UpdateH(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const Eigen::VectorXd& u, std::vector<std::vector<double>>& H, double E_int, double E_bulk, double nu_int, double nu_bulk, double h2, double hi, double ls, bool isPlaneStress, bool isQuadratic, int quadratureDegree, std::vector<std::vector<double>>&exx, std::vector<std::vector<double>>&eyy, std::vector<std::vector<double>>&exy, std::vector<std::vector<double>>&tr_e)
{
    // 1. Calculate Quadrature ONCE (Performance Fix)
    Eigen::VectorXd xi_points, eta_points, weights;
    Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);

    // 2. Allocate reusable vectors ONCE (Performance Fix)
    int num_points = weights.size();
    std::vector<double> H_new;
    std::vector<double> exx_i;
    std::vector<double> eyy_i;
    std::vector<double> exy_i;
    std::vector<double> tr_e_i;

    // Reserve memory to prevent re-allocation
    H_new.reserve(num_points);
    exx_i.reserve(num_points);
    eyy_i.reserve(num_points);
    exy_i.reserve(num_points);
    tr_e_i.reserve(num_points);

    std::cout << "Updating H in element " << std::flush;
    
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
            double x1 = nodes[n1][0]; double y1 = nodes[n1][1];
            double x2 = nodes[n2][0]; double y2 = nodes[n2][1];
            double x3 = nodes[n3][0]; double y3 = nodes[n3][1];
            double x4 = nodes[n4][0]; double y4 = nodes[n4][1];
            double x5 = nodes[n5][0]; double y5 = nodes[n5][1];
            double x6 = nodes[n6][0]; double y6 = nodes[n6][1];
            double x7 = nodes[n7][0]; double y7 = nodes[n7][1];
            double x8 = nodes[n8][0]; double y8 = nodes[n8][1];
            double x9 = nodes[n9][0]; double y9 = nodes[n9][1];
            
            double E = (y9 > h2 && y9 < h2 + hi) ? E_int : E_bulk;
            double nu = (y9 > h2 && y9 < h2 + hi) ? nu_int : nu_bulk;
            
            Eigen::VectorXd u_el(18); 
            u_el << u(2*n1), u(2*n1 + 1), u(2*n2), u(2*n2 + 1), u(2*n3), u(2*n3 + 1), u(2*n4), u(2*n4 + 1), u(2*n5), u(2*n5 + 1), u(2*n6), u(2*n6 + 1), u(2*n7), u(2*n7 + 1), u(2*n8), u(2*n8 + 1), u(2*n9), u(2*n9 + 1);
            
            // Pass the pre-allocated vectors by reference
            GetQuadraticElementTensileElasticStrainEnergy(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8, x9, y9, E, nu, u_el, isPlaneStress, xi_points, eta_points, weights, H_new, exx_i, eyy_i, exy_i, tr_e_i);
            
            // Update global H with irreversibility (max)
            std::vector<double>& H_old = H[i];
            for(int k = 0; k < num_points; ++k) {
                H_old[k] = std::max(H_old[k], H_new[k]);
            }
            
            // Store strains
            exx[i] = exx_i;
            eyy[i] = eyy_i;
            exy[i] = exy_i;
            tr_e[i] = tr_e_i;
        }
        else // linear elements (4 nodes)
        {
            int n1 = elements[i][0];
            int n2 = elements[i][1];
            int n3 = elements[i][2];
            int n4 = elements[i][3];
            double x1 = nodes[n1][0]; double y1 = nodes[n1][1];
            double x2 = nodes[n2][0]; double y2 = nodes[n2][1];
            double x3 = nodes[n3][0]; double y3 = nodes[n3][1];
            double x4 = nodes[n4][0]; double y4 = nodes[n4][1];
            
            double yc = (y1 + y2 + y3 + y4)/4; 
            double E = (yc > h2 && yc < h2 + hi) ? E_int : E_bulk;
            double nu = (yc > h2 && yc < h2 + hi) ? nu_int : nu_bulk;
            
            Eigen::VectorXd u_el(8); 
            u_el << u(2*n1), u(2*n1 + 1), u(2*n2), u(2*n2 + 1), u(2*n3), u(2*n3 + 1), u(2*n4), u(2*n4 + 1);
            
            // Pass the pre-allocated vectors by reference
            GetLinearElementTensileElasticStrainEnergy(x1, y1, x2, y2, x3, y3, x4, y4, E, nu, u_el, isPlaneStress, xi_points, eta_points, weights, H_new, exx_i, eyy_i, exy_i, tr_e_i);
            
            // Update global H with irreversibility (max)
            std::vector<double>& H_old = H[i];
            for(int k = 0; k < num_points; ++k) {
                H_old[k] = std::max(H_old[k], H_new[k]);
            }

            // Store strains
            exx[i] = exx_i;
            eyy[i] = eyy_i;
            exy[i] = exy_i;
            tr_e[i] = tr_e_i;
        }
        
        if((i + 1) % 10000 == 0) std::cout << i + 1 << ", " << std::flush;
        else if(i == elements.size() - 1) std::cout << i + 1 << std::endl;
    }
}

// System reduction
void ApplyDirichletBC(const Eigen::SparseMatrix<double> &K, const Eigen::VectorXd &F, const std::vector<bool> &isConstrained, const std::unordered_map<int, double*> &constrainedMap, const std::vector<int> &freeMap, Eigen::SparseMatrix<double> &KR, Eigen::VectorXd &FR)
{
    FR.setZero();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(K.nonZeros());

    // assemble reduced matrix and force vector
    for(int col = 0; col < K.outerSize(); ++col)
	{
        for(Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it)
		{
            int row = it.row();
            double val = it.value();
            int r = freeMap[row];
            int c = freeMap[col];

            if(r != -1 && c != -1)
			{
				// free-free
                triplets.emplace_back(r, c, val);
            }
			else if(r != -1 && c == -1)
			{
                // free row, constrained column -> subtract K_rc*prescribed_value
				FR[r] -= val*(*constrainedMap.at(col));
            }
            // if r == -1 (row constrained) -> nothing to add for KR row
        }
    }

    KR.setFromTriplets(triplets.begin(), triplets.end());
	
    // now add free-part of F
    for(int fullRow = 0; fullRow < (int)freeMap.size(); ++fullRow)
	{
        int r = freeMap[fullRow];
        if (r != -1) FR[r] += F[fullRow];
    }
}

template<typename SparseMatrixType, typename PreconditionerType>
double estimateConditionNumber(
    const SparseMatrixType& A,
    const PreconditionerType& M,
    int maxIter)
{
    using Vec = Eigen::VectorXd;
    const int n = A.rows();

    if (A.rows() != A.cols()) {
        std::cerr << "[estimateConditionNumber] Matrix must be square!" << std::endl;
        return -1.0;
    }

    // --- Helper lambda: apply preconditioned operator ---
    auto applyPrecondA = [&](const Vec& x) -> Vec {
        // Apply M^{-1/2} * A * M^{-1/2} * x
        // Approximated by: M.solve(A * M.solve(x))
        Vec y = M.solve(x);
        Vec z = A * y;
        return M.solve(z);
    };

    // --- LLT factorization (robust SPD solver) ---
    Eigen::SimplicialLLT<SparseMatrixType> llt;
    llt.compute(A);
    if (llt.info() != Eigen::Success) {
        std::cerr << "[estimateConditionNumber] LLT factorization failed!" << std::endl;
        return -1.0;
    }

    // --- Power iteration: estimate λ_max ---
    Vec x = Vec::Random(n).normalized();
    double lambda_max = 0.0;
    for (int i = 0; i < maxIter; ++i) {
        Vec y = applyPrecondA(x);
        lambda_max = x.dot(y);
        x = y.normalized();
    }

    // --- Inverse iteration: estimate λ_min ---
    Vec v = Vec::Random(n).normalized();
    double lambda_min = 0.0;
    for (int i = 0; i < maxIter; ++i) {
        Vec y = applyPrecondA(v);
        Vec w = llt.solve(y); // stable solve
        if (llt.info() != Eigen::Success) {
            std::cerr << "[estimateConditionNumber] LLT solve failed in inverse iteration." << std::endl;
            return -1.0;
        }
        v = w.normalized();
        lambda_min = v.dot(applyPrecondA(v));
    }

    double cond_est = lambda_max / std::max(lambda_min, 1e-16);
    return cond_est;
}
