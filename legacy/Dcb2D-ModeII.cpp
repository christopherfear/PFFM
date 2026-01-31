#include "Dcb2D-ModeII.h"

#include <cholmod.h>
#include <Eigen/CholmodSupport>

int main()
{
	// dimensions
	double L = 0.1; // length of DCB (m)
	double h1 = 2e-3; // height of upper DCB arm (m)
	double h2 = 2e-3; // height of lower DCB arm (m)
	double hi = 1e-4; // height of interface (m)
	double a0 = 0.03; // initial crack length (m)
	double ls = 20e-6; // length scale
	double t = 1; // thickness
	double a_min = a0 + 3.7e-5; // minimum crack length before compliance calibration or dPi/dA fitting is started
	
	// material properties
	double E_int = 126.e9; // interface Young's modulus
	double nu_int = 0.3; // interface Poisson's ratio
	double Gc_int = 562.; // interface fracture toughness for beam theory calculation only
	double Gc_eff = 497.63; // effective interface fracture toughness for FEM only
	double E_bulk = 126.e9; // bulk Young's modulus
	double nu_bulk = 0.3; // bulk Poisson's ratio
	double Gc_bulk = 10.*Gc_int; // bulk fracture toughness
	bool isPlaneStress = true; // if false, plane strain
	
	// mesh properties
	double Delta_y = L/3000; // initial element length at crack tip
	double Delta_y_max_ahead = L/10.; // max element length ahead of crack tip
	double Delta_y_max_behind = 1e-4; // max element length behind crack tip
	double GR_y_ahead = 1.1; // growth rate of element length ahead of crack tip
	double GR_y_behind = 1.2; // growth rate of element length behind crack tip
	double Delta = hi/24; // interface element height
	double Delta_x = Delta; // initial element height in bulk
	double Delta_x_max = h1/10.; // max element height in bulk
	double yminus = 4e-3; // length of refined region ahead of crack tip
	double yplus = 1e-3; // length of refined region behind crack tip
	double xplus = hi/2; // length of refined region above interface
	double xminus = hi/2; // length of refined region below interface
	double GR_x = 1.5; // growth rate of element height in bulk
	bool isQuadratic = true; // if false, use linear elements
	int quadratureDegree = 5; // degree of polynomial that can be fully integrated
	
	// \/run parameters
	double tol = 1e-7; // Newton's method convergence tolerance
	double relax = 1.; // fraction of the increment to apply each iteration
	int save_freq = 500; // number of increments per files getting saved (1 = every increment, 2 = even increments, etc.)
	int restart_from = 0; // iteration number to restart from
	
	// output files
	std::string Hfilename = "H,inc=";
	std::string usfilename = "inc=";
	std::string history_filename = "history.csv";
	int prec = 16; // number of digits after the decimal point
	
	// loading
	double Etheory = isPlaneStress ? E_bulk : E_bulk/(1. - nu_bulk*nu_bulk);
	double h1i = h1 + hi/2; // thickness of upper arm including half interface
	double h2i = h2 + hi/2; // thickness of lower arm including half interface
	double Dinit = 8.48e-3 ; // initial displacement
	double Dend = 100e-3; // applied displacement
	double Dinc = 2.5e-9; // initial displacement increment
	
	
	/***************/
	/* CREATE MESH */
	/***************/
	// x coordinates from support to refined region
	std::vector<double> x1 = linspace2(0, L - a0 - yminus, Delta_y, GR_y_ahead, Delta_y_max_ahead);
	std::reverse(x1.begin(), x1.end());
	for(auto& elem : x1)
	{
		elem *= -1;
		elem += (L - a0 - yminus);
	}
	
	std::vector<double> x2 = linspace(L - a0 - yminus, L - a0, Delta_y); // x coordinates in refined region ahead of crack tip
	std::vector<double> x3 = linspace(L - a0, L - a0 + yplus, Delta_y); // x coordinates in refined region behind of crack tip
	std::vector<double> x4 = linspace2(L - a0 + yplus, L, Delta_y, GR_y_behind, Delta_y_max_behind); // x coordinates in cracked region
	
	std::vector<double> y6 = linspace2(h2 + hi + xplus, h2 + hi + h1, Delta_x, GR_x, Delta_x_max); // y coordinates in upper arm
	std::vector<double> y5 = linspace(h2 + hi, h2 + hi + xplus, Delta); // y coordinates in refined region above interface
	std::vector<double> y4 = linspace(h2 + hi/2, h2 + hi, Delta); // y coordinates in upper interface
	std::vector<double> y3 = linspace(h2, h2 + hi/2, Delta); // y coordinates in lower interface
	std::vector<double> y2 = linspace(h2 - xplus, h2, Delta); // y coordinates in refined region below interface
	
	// y coordinates in lower arm
	std::vector<double> y1 = linspace2(0, h2 - xplus, Delta_x, GR_x, Delta_x_max); // y coordinates in lower arm, mirrored about y = h2/2
	std::reverse(y1.begin(), y1.end());
	for(auto& elem : y1)
	{
		elem *= -1;
		elem += (h2 - xplus);
	}
	
	// merged x coordinates
	std::vector<double> x = x1;
	x.insert(x.end(), x2.begin() + 1, x2.end());
	x.insert(x.end(), x3.begin() + 1, x3.end());
	x.insert(x.end(), x4.begin() + 1, x4.end());
	
	// merged y coordinates
	std::vector<double> y = y1;
	y.insert(y.end(), y2.begin() + 1, y2.end());
	y.insert(y.end(), y3.begin() + 1, y3.end());
	y.insert(y.end(), y4.begin() + 1, y4.end());
	y.insert(y.end(), y5.begin() + 1, y5.end());
	y.insert(y.end(), y6.begin() + 1, y6.end());
	
	// insert mid-points if using quadratic elements
	if(isQuadratic)
	{
		insert_midpoints(x);
		insert_midpoints(y);
	}

	// create grid of coordinates
	std::vector<std::vector<double>> X, Y;
	meshgrid(x, y, X, Y);
	
	// create nodes matrix
	std::vector<std::vector<double>> nodes;
	create_nodes_matrix(x, y, nodes);
	
	// create elements matrix
	std::vector<std::vector<int>> elements;
	create_elements_matrix(X, Y, isQuadratic, elements);
	
	// split mesh
	double findtol = Delta/10;
	add_crack(nodes, elements, L - a0, h2 + hi/2, findtol);
	
	int numElements = elements.size();
	int numNodes = nodes.size();
	std::cout << "Number of elements: " << numElements << std::endl;
	std::cout << "Number of nodes: " << numNodes << std::endl << std::endl;
	
	
	/************************/
	/* SAVE GEOMETRY TO CSV */
	/************************/
	std::ofstream outfile("geometry.csv"); // overwrite
	outfile << "L,h1,h2,hi,a0" << std::endl;
	outfile << std::fixed << std::setprecision(prec) << L
		<< "," << h1
		<< "," << h2
		<< "," << hi
		<< "," << a0
		<< std::endl;
	outfile.close();
	
	
	/*************/
	/* NODE SETS */
	/*************/
	std::vector<int> initial_crack_nodes;
	std::vector<int> centreline_nodes; // excluding initial crack
	std::vector<int> interface_nodes; // excluding initial crack
	std::vector<int> upper_arm_end_nodes;
	std::vector<int> lower_arm_end_nodes;
	std::vector<int> support_nodes;
	for(size_t i = 0; i < nodes.size(); ++i)
	{
		double x = nodes[i][0];
		double y = nodes[i][1];
		
		// set of nodes forming the initial crack
		if((x >= (L - a0) - findtol) && (std::abs(y - (h2 + hi/2)) < findtol))
		{
			initial_crack_nodes.push_back(i);
		}
		
		// set of nodes forming the centreline of the interface
		if((x <= (L - a0) + findtol) && (std::abs(y - (h2 + hi/2)) < findtol))
		{
			centreline_nodes.push_back(i);
		}
		
		// set of all nodes forming the interface
		if((x <= (L - a0) + findtol) && (y >= h2 - findtol) && (y <= h2 + hi + findtol))
		{
			interface_nodes.push_back(i);
		}
		
		// set of nodes on end of upper arm
		if((y >= (h2 + hi) - findtol) && (std::abs(x - L) < findtol))
		{
			upper_arm_end_nodes.push_back(i);
		}
		
		// set of nodes on end of lower arm
		if((y <= h2 + findtol) && (std::abs(x - L) < findtol))
		{
			lower_arm_end_nodes.push_back(i);
		}
		
		// set of nodes on support end
		if((std::abs(x) < findtol))
		{
			support_nodes.push_back(i);
		}
	}

	/**********************/
	/* INITIAL CONDITIONS */
	/**********************/
	Eigen::VectorXd u = Eigen::VectorXd::Zero(2*numNodes);
	Eigen::VectorXd s = Eigen::VectorXd::Ones(numNodes);
	std::vector<std::vector<double>> H; // H[i][j] = H for element i at integration point j
	std::vector<std::vector<double>> exx; std::vector<std::vector<double>> eyy; std::vector<std::vector<double>> exy; std::vector<std::vector<double>> tr_e;
	std::vector<double> F1_history, D1_history, F2_history, D2_history, a_history, GcInit_CC_history;
	std::vector<double> F1_beamtheory, F2_beamtheory, G_beamtheory, a_theory;
	std::vector<double> U_history, Gc_dPidA_history;
	int incr_count = 0;
	if(restart_from > 0) // load in starting condition from saved files
	{
		incr_count = restart_from;
		
		// lambda to split a line by commas
		auto splitByComma = [](const std::string& line) -> std::vector<std::string>
		{
			std::vector<std::string> result;
			std::stringstream ss(line);
			std::string item;
			while (std::getline(ss, item, ','))
			{
				result.push_back(item);
			}
			return result;
		};

		// read u and s from usfilename
		std::string filename = usfilename;
		filename += std::to_string(incr_count);
		filename += ".csv";
		std::cout << "Reading in u, v and s from " << filename << " ... ";
		std::ifstream file(filename);
		std::string line;
		if(!file.is_open())
		{
			std::cerr << "Could not open " << filename << "!" << std::endl;
			return 1;
		}
		std::getline(file, line); // skip header row
		unsigned int node_index = 0;
		while(std::getline(file, line))
		{
			std::vector<std::string> fields = splitByComma(line);
			for(size_t i = 0; i < fields.size(); i++)
			{
				if(i == 2) u(2*node_index) = std::stod(fields[i]);
				else if(i == 3) u(2*node_index + 1) = std::stod(fields[i]);
				else if(i == 4) s(node_index) = std::stod(fields[i]);
			}
			node_index++;
		}
		file.close();
		std::cout << "Done!" << std::endl;
		
		// read H from Hfilename
		filename = Hfilename;
		filename += std::to_string(incr_count);
		filename += ".csv";
		std::cout << "Reading in H from " << filename << " ... ";
		file = std::ifstream(filename);
		if(!file.is_open())
		{
			std::cerr << "Could not open " << filename << "!" << std::endl;
			return 1;
		}
		std::getline(file, line); // skip header row
		Eigen::VectorXd xi_points, eta_points, weights; // integration points and weightings
		Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
		while(std::getline(file, line))
		{
			std::vector<double> Hi;
			for(unsigned int i = 0; i < weights.size(); ++i)
			{
				if(i > 0) std::getline(file, line);
				std::vector<std::string> fields = splitByComma(line);
				Hi.push_back(std::stod(fields[2]));
			}
			H.push_back(Hi);
		}
		file.close();
		std::cout << "Done!" << std::endl;
		
		// read from history
		std::cout << "Reading in history from " << history_filename << " ... ";
		file = std::ifstream(history_filename);
		if(!file.is_open())
		{
			std::cerr << "Could not open " << filename << "!" << std::endl;
			return 1;
		}
		std::getline(file, line); // skip header row
		while(std::getline(file, line))
		{
			std::vector<std::string> fields = splitByComma(line);
			for(size_t i = 0; i < fields.size(); i++)
			{
				if(i == 0)
				{
					int inc = std::stoi(fields[i]);
					if(inc > restart_from) break;
				}
				
				double val = std::stod(fields[i]);
				
				if(i == 1)
				{
					D1_history.push_back(val);
					Dinit = val;
				}
				else if(i == 2) F1_history.push_back(val);
				else if(i == 3) D2_history.push_back(val);
				else if(i == 4) F2_history.push_back(val);
				else if(i == 5) a_history.push_back(val);
				else if(i == 7) U_history.push_back(val);
				
				// recalculate Gc_initiation from compliance calibration based on history to this point and value of a_min
				if(i == 5)
				{
					GcInit_CC_history.push_back(GetGcInitiationFromComplianceCalibration(D1_history, F1_history, F2_history, a_history, a_min));
				}
			}
		}
		file.close();
		std::cout << "Done!" << std::endl;
		Dinit += Dinc; // start on next increment
		incr_count++; // start on next increment
	}
	else
	{
		for(auto& node : initial_crack_nodes)
		{
			s[node] = 0;
		}
		
		Eigen::VectorXd xi_points, eta_points, weights; // integration points and weightings
		Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
		for(size_t i = 0; i < numElements; ++i)
		{
			std::vector<double> Hi(weights.size(), 0.); // vector of zeros to allocate memory
			H.push_back(Hi);
		}
		
	}
	
	// Allocated memory for strains
	Eigen::VectorXd xi_points, eta_points, weights; // integration points and weightings
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	for(size_t i = 0; i < numElements; ++i)
	{
		std::vector<double> Hi(weights.size(), 0.); // vector of zeros to allocate memory
		exx.push_back(Hi);
		eyy.push_back(Hi);
		exy.push_back(Hi);
		tr_e.push_back(Hi);
	}
	

	/***********************/
	/* SOLVE INCREMENTALLY */
	/***********************/
	// persistent solver instances
	static Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solverCM;
	static Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> solverLLT;
	static Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> solverLDLT;
	static bool firstLoop = true;

	// Pre-allocation of full vectors and matrices
	Eigen::SparseMatrix<double> K(3*numNodes, 3*numNodes); // global stiffness matrix
	Eigen::VectorXd F(K.rows()); // global force vector
	Eigen::VectorXd d(K.rows()); // combined u and s vector
	Eigen::VectorXd unm1(u.size()); // solution from previous increment
	Eigen::VectorXd Dd(K.rows()); // change in solution vector
	// boundary conditions
	std::vector<int> constrained; // constrained DOF indices
	for(int node : initial_crack_nodes)
	{
		constrained.push_back(2*numNodes + node); // damage of the initial crack nodes
	}
	for(int node : support_nodes)
	{
		constrained.push_back(2*node); // horizontal displacement of the support nodes
		constrained.push_back(2*node + 1); // vertical displacement of the support nodes
	}
	for(int node : upper_arm_end_nodes)
	{
		constrained.push_back(2*node + 1); // vertical displacement of upper arm nodes
	}
	for(int node : lower_arm_end_nodes)
	{
		constrained.push_back(2*node + 1); // vertical displacement of lower arm nodes
	}
	std::vector<double> constrainedValues(constrained.size(), 0.); // all constrained DOFs have zero value (applied displacements are updated later)
	
	// build map of constrained DOF vs pointer to value in constrainedValues
	std::unordered_map<int, double*> constrainedMap;
	constrainedMap.reserve(constrained.size());
	for(size_t k = 0; k < constrained.size(); ++k)
	{
		constrainedMap[constrained[k]] = &constrainedValues[k];
	}
	
	// build lookup indicating if a DOF is constrained / unconstrained
    std::vector<bool> isConstrained(K.rows(), false);
    for(const auto &pair : constrainedMap)
	{
        isConstrained[pair.first] = true;
	}
	
    // build lookup translating DOF in full system to DOF in reduced system
	std::vector<int> freeMap(K.rows(), -1);
    int nFree = 0;
    for(int i = 0; i < K.rows(); ++i)
	{
        if(!isConstrained[i]) freeMap[i] = nFree++;
	}
	
	// Pre-allocation of reduced vectors and matrices
	Eigen::SparseMatrix<double> KR(nFree, nFree); // global stiffness matrix
	Eigen::VectorXd d_reduced(nFree); // reduced d vector
	Eigen::VectorXd FR(nFree); // reduced force vector
	Eigen::VectorXd Dd_reduced(nFree); // change in reduced solution vector
	Eigen::VectorXd residual(nFree); // residual vector in Newton's method (KR*d)

	// apply displacement incrementally
	for(double D = Dinit; D <= Dend + Dinc/10; D += Dinc)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		printf("Increment number: %d\n", incr_count);
		printf(("Applied displacement: %." + std::to_string(prec) + "f\n").c_str(), D);
		
		// update boundary conditions
		size_t startIdx = initial_crack_nodes.size() + 2*support_nodes.size();
		std::fill(constrainedValues.begin() + startIdx, constrainedValues.end(), D);
		
		// combine u and s into combined solution vector
		d.head(u.size()) = u; // assign u to first part of d
		d.tail(s.size()) = s; // assign s to second part of d
		
		unm1 = u; // update previous increment's solution u_{inc-1} = u_inc
		int iter_count = 1;
		while(true)
		{
			printf("Iteration: %d\n", iter_count);
			
			// assemble the problem
			AssembleAll(nodes, elements, u, s, H, E_int, E_bulk, nu_int, nu_bulk, Gc_eff, Gc_bulk, t, h2, hi, ls, K, F, isPlaneStress, isQuadratic, quadratureDegree);

			// apply boundary conditions
			ApplyDirichletBC(K, F, isConstrained, constrainedMap, freeMap, KR, FR);
			
			// build reduced d vector (values at free DOFs)
			for(int i = 0; i < K.rows(); ++i)
			{
				int idx = freeMap[i];
				if(idx != -1) d_reduced[idx] = d[i];
			}
			
			// calculate the next Newton's method iteration
			std::cout << "Solving... " << std::flush;
			residual = FR - KR*d_reduced;
			bool solved = false;
			
			// reuse symbolic analysis
			if(firstLoop)
			{
				solverCM.analyzePattern(KR);
				solverLLT.analyzePattern(KR);
				solverLDLT.analyzePattern(KR);
				firstLoop = false;
			}
			
			// Try CholMod first (should be fastest for SPD)
			solverCM.factorize(KR);
			if(solverCM.info() == Eigen::Success)
			{
				Dd_reduced = solverCM.solve(residual);
				if(solverCM.info() == Eigen::Success)
				{
					std::cout << "[CHOLMOD] succeeded. " << std::flush;
					solved = true;
				}
				else
				{
					std::cerr << "[CHOLMOD] solve failed. " << std::flush;
				}
			}
			else
			{
				std::cerr << "[CHOLMOD] decomposition failed. " << std::flush;
			}

			// try SimplicialLLT (very fast for SPD)
			if(!solved)
			{
				solverLLT.factorize(KR);
				if(solverLLT.info() == Eigen::Success)
				{
					Dd_reduced = solverLLT.solve(residual);
					if(solverLLT.info() == Eigen::Success)
					{
						std::cout << "[LLT] succeeded. " << std::flush;
						solved = true;
					}
					else
					{
						std::cerr << "[LLT] solve failed. " << std::flush;
					}
				}
				else
				{
					std::cerr << "[LLT] decomposition failed. " << std::flush;
				}
			}
			
			// if LLT failed, try SimplicialLDLT (robust SPD/indefinite)
			if(!solved)
			{
				solverLDLT.factorize(KR);
				if(solverLDLT.info() == Eigen::Success)
				{
					Dd_reduced = solverLDLT.solve(residual);
					if (solverLDLT.info() == Eigen::Success)
					{
						std::cout << "[LDLT] succeeded. " << std::flush;
						solved = true;
					}
					else
					{
						std::cerr << "[LDLT] solve failed. " << std::flush;
					}
				}
				else
				{
					std::cerr << "[LDLT] decomposition failed. " << std::flush;
				}
			}
			
			// if other solvers both failed, fall back to SparseLU
			if(!solved)
			{
				Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int>> solverLU;
				solverLU.compute(KR);
				if (solverLU.info() == Eigen::Success)
				{
					Dd_reduced = solverLU.solve(residual);
					if (solverLU.info() == Eigen::Success)
					{
						std::cout << "[SparseLU] succeeded. " << std::flush;
						solved = true;
					}
					else
					{
						std::cerr << "[SparseLU] solve failed!" << std::endl;
					}
				}
				else
				{
					std::cerr << "[SparseLU] decomposition failed!" << std::endl;
				}
			}
			
			// final check
			if(!solved)
			{
				std::cerr << "All solvers failed! Aborting." << std::endl;
				return -1;
			}
			
			// Reconstruct full Dd (increments) of length N
			for(int i = 0; i < K.rows(); ++i)
			{
				int idx = freeMap[i];
				if (idx != -1) Dd[i] = Dd_reduced[idx];
				else Dd[i] = 0.0; // constrained DOFs get zero increment
			}

			d += relax*Dd; // apply the displacement/damage increment
			for(const auto& pair : constrainedMap)
			{
				d[pair.first] = *(pair.second); // apply prescribed Dirichlet values
			}
			u = d.head(u.size());
			s = d.tail(s.size());
			std::cout << "Done! " << std::flush;
			
			// convergence checks
			double Rk = residual.norm();
			double R_normalised = Rk/FR.norm();
			printf("|R|/|F| = %.3e\n", R_normalised);
			
			// stop iterating if norm within tolerance
            if(R_normalised < tol)
			{
				/*******************/
				/* POST PROCESSING */
				/*******************/
				
				// compute force vector
				F = K*d;
				
				// sum vertical forces on upper arm end nodes
				double F1 = 0.0;
				for(int node : upper_arm_end_nodes) F1 += F(2*node + 1);
				F1_history.push_back(F1);
				D1_history.push_back(D);
				
				// sum vertical forces on lower arm end nodes
				double F2 = 0.0;
				for(int node : lower_arm_end_nodes) F2 += F(2*node + 1);
				F2_history.push_back(F2);
				D2_history.push_back(D);

				// crack length
				double intact_length = L;
				double s_threshold = 0.01;
				std::vector<double> processed_y; // rows already processed
				std::vector<double> processed_x; // corresponding x where calculation was done
				for(size_t i = 0; i < interface_nodes.size(); ++i)
				{
					int n2 = interface_nodes[i];
					double s2 = d(2*numNodes + n2);
					double x2 = nodes[n2][0];
					double y2 = nodes[n2][1];

					// find if this row was processed before
					int row_index = -1;
					for(size_t j = 0; j < processed_y.size(); ++j)
					{
						if(std::fabs(y2 - processed_y[j]) < Delta/10) // almost equal to within 10% of the interface element height
						{
							row_index = j;
							break;
						}
					}

					// case 1: row was already processed - check x
					if(row_index != -1)
					{
						// if current x >= stored x, skip (we already have an earlier crossing)
						if(x2 >= processed_x[row_index]) continue;
						// otherwise, allow recalculation since this is further left
					}

					// case 2: perform calculation if threshold is met
					if(s2 <= s_threshold)
					{
						int n1 = interface_nodes[i-1];
						double x1 = nodes[n1][0];
						double y1 = nodes[n1][1];
						double s1 = d(2*numNodes + n1);

						double current_intact_length = x1 + (s_threshold - s1)/(s2 - s1)*(x2 - x1);
						if(current_intact_length < intact_length) intact_length = current_intact_length;

						// if row was already seen, update stored x; else add a new entry
						if(row_index != -1) processed_x[row_index] = x2;
						else
						{
							processed_y.push_back(y2);
							processed_x.push_back(x2);
						}
					}
				}

				a_history.push_back(L - intact_length);
				// std::cout << "crack length: " << (L - intact_length) << std::endl;
			
				// calculate Gc from compliance calibration based on history to this point and value of a_min
				double GcInit_CC = GetGcInitiationFromComplianceCalibration(D1_history, F1_history, F2_history, a_history, a_min);
				GcInit_CC_history.push_back(GcInit_CC);
				std::cout << "Gc from compliance calibration = " << GcInit_CC << std::endl;
				
				// calculate strain energy
				double U = 0.5*u.dot(K.topLeftCorner(u.size(), u.size())*u);
				U_history.push_back(U);

				// save increment to CSV
				std::cout << "Saving history.csv" << std::endl;
				std::ofstream outfile(history_filename); // overwrite
				outfile << "inc,D1,F1,D2,F2,a,G_CC,U" << std::endl;
				for(size_t i = 0; i < F1_history.size(); ++i)
				{
					outfile << i
						<< "," << std::fixed << std::setprecision(prec) << D1_history[i]
						<< "," << F1_history[i]
						<< "," << D2_history[i]
						<< "," << F2_history[i]
						<< "," << a_history[i]
						<< "," << GcInit_CC_history[i]
						<< "," << U_history[i]
						<< std::endl;
				}
				outfile.close();
				
				break;
			}
			
			++iter_count;
		}
		
		// update H ready for next increment
		UpdateH(nodes, elements, u, H, E_int, E_bulk, nu_int, nu_bulk, h2, hi, ls, isPlaneStress, isQuadratic, quadratureDegree, exx, eyy, exy, tr_e);
		
		// save H output file
		if(incr_count % save_freq == 0)
		{
			std::string filename = Hfilename;
			filename += std::to_string(incr_count);
			filename += ".csv";
			std::cout << "Saving " << filename << std::endl;
			SaveHField(nodes, elements, H, isQuadratic, quadratureDegree, filename, prec, exx, eyy, exy, tr_e);
		}
		
		// save u and s output file
		if(incr_count % save_freq == 0)
		{

			std::string filename = usfilename;
			filename += std::to_string(incr_count);
			filename += ".csv";
			std::cout << "Saving " << filename << std::endl;
			std::ofstream outfile(filename); // overwrite
			outfile << "x,y,u,v,s" << std::endl;
			for(size_t i = 0; i < nodes.size(); ++i)
			{
				outfile << std::fixed << std::setprecision(prec) << nodes[i][0]
					<< "," << nodes[i][1]
					<< "," << u[2*i]
					<< "," << u[2*i + 1]
					<< "," << s[i]
					<< std::endl;
			}
			outfile.close();
		}

		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> seconds = t2 - t1;
		printf("Increment %d took %.3f seconds\n", incr_count, seconds.count());
		std::cout << std::endl;
		
		++incr_count;
	}

	return 0;
}

std::vector<double> linspace(double start, double end, double initial_size, double growth_rate, double tol)
{
	std::vector<double> values;
	
	if(start == end || initial_size <= 0 || growth_rate <= 0 || initial_size > end - start)
	{
		throw std::invalid_argument("Invalid parameters: Ensure start != end, initial_size > 0 and initial_size < end - start, and growth_rate > 0.");
	}

	values.push_back(start);
	
	// compute number of steps for growth_rate = 1
	int num_steps;
	double adjusted_growth;
	if(growth_rate != 1)
	{
		// compute number of steps required based on growth rate
		num_steps = std::round(std::log(1. - (end - start)/initial_size*(1. - growth_rate))/std::log(growth_rate));
		
		// use Newton Raphson method to adjust growth rate to ensure final value matches exactly
		adjusted_growth = growth_rate;
		for(int i = 0; i < 1e6; ++i)
		{
			double f = initial_size*(1. - std::pow(adjusted_growth, num_steps))/(1. - adjusted_growth) - (end - start);
			double f_deriv = (initial_size*(-num_steps*std::pow(adjusted_growth, num_steps - 1.)*(1. - adjusted_growth) - (1. - std::pow(adjusted_growth, num_steps))))/std::pow(1. - adjusted_growth, 2.);
			adjusted_growth -= f/f_deriv;
			if(std::fabs(f) < tol) break; // converged
		}
	}
	else
	{
		num_steps = std::round((end - start)/initial_size);
		initial_size = (end - start)/num_steps;
		adjusted_growth = 1;
	}

	// generate the values
	double current_size = initial_size;
	for(int i = 0; i < num_steps - 1; ++i)
	{
		double next_value = values.back() + current_size;
		values.push_back(next_value);
		current_size = current_size*adjusted_growth; // update for next iteration
	}
	
	values.push_back(end); // avoid cumulative floating point error by adding precise value of end

	return values;
}

std::vector<double> linspace2(double start, double end, double initial_size, double growth_rate, double max_size, double tol)
{
	// initial pass before any tuning
	double x = start;
	double current_size = initial_size;
	int count_max_size = 0;
	while(x < end)
	{
		if(current_size <= max_size)
		{
			x += current_size;
			current_size = current_size*growth_rate;
		}
		else
		{
			count_max_size++;
			x += max_size;
		}
	}
	if(count_max_size > 0) count_max_size--; // remove overshoot
	double newend = end - count_max_size*max_size; // reduce total length by number of max_size elements
	
	std::vector<double> values;
	
	if(start == end || initial_size <= 0 || growth_rate <= 0 || initial_size > end - start)
	{
		throw std::invalid_argument("Invalid parameters: Ensure start != end, initial_size > 0 and initial_size < end - start, and growth_rate > 0.");
	}

	values.push_back(start);
	
	// compute number of steps for growth_rate = 1
	int num_steps;
	double adjusted_growth;
	if(growth_rate != 1)
	{
		// compute number of steps required based on growth rate
		num_steps = std::round(std::log(1. - (newend - start)/initial_size*(1. - growth_rate))/std::log(growth_rate));
		
		// use Newton Raphson method to adjust growth rate to ensure final value matches exactly
		adjusted_growth = growth_rate;
		for(int i = 0; i < 1e6; ++i)
		{
			double f = initial_size*(1. - std::pow(adjusted_growth, num_steps))/(1. - adjusted_growth) - (newend - start);
			double f_deriv = (initial_size*(-num_steps*std::pow(adjusted_growth, num_steps - 1.)*(1. - adjusted_growth) - (1. - std::pow(adjusted_growth, num_steps))))/std::pow(1. - adjusted_growth, 2.);
			adjusted_growth -= f/f_deriv;
			if(std::fabs(f) < tol) break; // converged
		}
	}
	else
	{
		num_steps = std::round((newend - start)/initial_size);
		initial_size = (newend - start)/num_steps;
		adjusted_growth = 1;
	}

	// generate the values
	current_size = initial_size;
	for(int i = 0; i < num_steps - 1; ++i)
	{
		double next_value = values.back() + current_size;
		values.push_back(next_value);
		current_size = current_size*adjusted_growth; // update for next iteration
	}
	
	values.push_back(newend);
	
	// add on max_size elements
	for(int i = 0; i < count_max_size - 1; ++i)
	{
		double next_value = values.back() + max_size;
		values.push_back(next_value);
	}
	
	if(count_max_size > 0) values.push_back(end); // avoid cumulative floating point error by adding precise value of end

	return values;
}

void insert_midpoints(std::vector<double>& vec)
{
	int original_size = vec.size();

	// iterate over the vector and insert midpoints between consecutive points
	for(int i = original_size - 1; i > 0; --i)
	{
		double midpoint = (vec[i - 1] + vec[i])/2.0;

		// insert the midpoint into the vector
		vec.insert(vec.begin() + i, midpoint);
	}
}

void meshgrid(const std::vector<double>& x, const std::vector<double>& y, std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y)
{
	int m = x.size();
	int n = y.size();
	
	// resize the output matrices to the correct size
	X.resize(n, std::vector<double>(m));
	Y.resize(n, std::vector<double>(m));

	// create the meshgrid
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < m; ++j)
		{
			X[i][j] = x[j]; // set X to repeat x values along rows
			Y[i][j] = y[i]; // set Y to repeat y values along columns
		}
	}
}

void create_elements_matrix(const std::vector<std::vector<double>>& X, const std::vector<std::vector<double>>& Y, bool isQuadratic, std::vector<std::vector<int>>& elements)
{
	int rows = X.size();
	int cols = X[0].size();

	int step;
	if(isQuadratic) step = 2;
	else step = 1;
	

	for (int i = 0; i < rows - 1; i += step)
	{
		for (int j = 0; j < cols - 1; j += step)
		{
			if(isQuadratic)
			{
				// quadratic elements (9 nodes) with midpoints already in X, Y
				int node1 = i*cols + j;
				int mid1 = i*cols + j + 1;
				int node2 = i*cols + j + 2;
				int mid2 = (i + 1)*cols + j + 2;
				int node3 = (i + 2)*cols + j + 2;
				int mid3 = (i + 2)*cols + j + 1;
				int node4 = (i + 2)*cols + j;
				int mid4 = (i + 1)*cols + j;
				int cen = (i + 1)*cols + j + 1;

				// store the 8-node quadratic element
				elements.push_back({node1, mid1, node2, mid2, node3, mid3, node4, mid4, cen});
			}
			else // linear elements (4 nodes)
			{
				int node1 = i*cols + j;
				int node2 = node1 + 1;
				int node3 = (i + 1)*cols + j + 1;
				int node4 = node3 - 1;
				
				elements.push_back({node1, node2, node3, node4}); // store the 4-node linear element
			}
		}
	}
}

void create_nodes_matrix(const std::vector<double>& x, const std::vector<double>& y, std::vector<std::vector<double>>& nodes)
{
	int m = x.size();
	int n = y.size();
	
	// create the nodes matrix with coordinates (x, y) values
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < m; ++j)
		{
			nodes.push_back({x[j], y[i]});
		}
	}
}

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

Eigen::MatrixXd LinearSolidElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double E, double nu, double t, Eigen::VectorXd u, Eigen::VectorXd s, bool plane_stress, int quadratureDegree)
{
	// constitutive properties matrix
	Eigen::MatrixXd D(3, 3);
	Eigen::MatrixXd Dvol(3, 3);
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
	Eigen::MatrixXd Ddev = D - Dvol;
	
	// integration points and weightings
	Eigen::VectorXd xi_points;
	Eigen::VectorXd eta_points;
	Eigen::VectorXd weights;
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	
	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(8, 8);
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
		double J_det = J.determinant();
		
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
		
		Eigen::Vector3d sigma = D*B*u;
		double sigma_vol = (sigma(0) + sigma(1))/2.; // volumetric stress (positive in tension)
		bool tensile = (sigma_vol > -std::numeric_limits<double>::epsilon()); // sigma_vol >= 0 with some tolerance for floating point error
		
		// integrate damage-affected stiffness matrix with quadrature
		if(tensile)
		{
			k += t*B_T*(N.transpose()*s*s.transpose()*N)*D*B*weight*J_det;
		}
		else // compressive
		{
			k += t*B_T*((N.transpose()*s*s.transpose()*N)*Ddev + Dvol)*B*weight*J_det;
		}
	}
	return k;
}

Eigen::MatrixXd QuadraticSolidElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, double E, double nu, double t, Eigen::VectorXd u, Eigen::VectorXd s, bool plane_stress, int quadratureDegree)
{
	// constitutive properties matrix
	Eigen::MatrixXd D(3, 3);
	Eigen::MatrixXd Dvol(3, 3);
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
	Eigen::MatrixXd Ddev = D - Dvol;
	
	// integration points and weightings
	Eigen::VectorXd xi_points;
	Eigen::VectorXd eta_points;
	Eigen::VectorXd weights;
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	
	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(18, 18);
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
		double J_det = J.determinant();
		
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
		
		Eigen::Vector3d sigma = D*B*u;
		double sigma_vol = (sigma(0) + sigma(1))/2.; // volumetric stress (positive in tension)
		bool tensile = (sigma_vol > -std::numeric_limits<double>::epsilon()); // sigma_vol >= 0 with some tolerance for floating point error
		
		// integrate damage-affected stiffness matrix with quadrature
		if(tensile)
		{
			k += t*B_T*(N.transpose()*s*s.transpose()*N)*D*B*weight*J_det;
		}
		else // compressive
		{
			k += t*B_T*((N.transpose()*s*s.transpose()*N)*Ddev + Dvol)*B*weight*J_det;
		}
	}
	return k;
}

void LinearPFFMElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double Gc, double ls, Eigen::VectorXd s, std::vector<double> H, int quadratureDegree, Eigen::MatrixXd& k, Eigen::VectorXd& f)
{
	// integration points and weightings
	Eigen::VectorXd xi_points;
	Eigen::VectorXd eta_points;
	Eigen::VectorXd weights;
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	
	k = Eigen::MatrixXd::Zero(4, 4);
	f = Eigen::VectorXd::Zero(4);
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
		double J_det = J.determinant();
		
		// shape function derivatives in global coordinates
		Eigen::MatrixXd dN_local(2, 4);
		dN_local << dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, 
			dN1_deta, dN2_deta, dN3_deta, dN4_deta;
		Eigen::MatrixXd dN_global = J_inv*dN_local; // rows = dNi/dx, dNi/dy; columns = dN1/d(x,y), dN2/d(x,y), ...
		
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

void QuadraticPFFMElement(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, double Gc, double ls, Eigen::VectorXd s, std::vector<double> H, int quadratureDegree, Eigen::MatrixXd& k, Eigen::VectorXd& f)
{
	// integration points and weightings
	Eigen::VectorXd xi_points;
	Eigen::VectorXd eta_points;
	Eigen::VectorXd weights;
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	
	k = Eigen::MatrixXd::Zero(9, 9);
	f = Eigen::VectorXd::Zero(9);
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
		double J_det = J.determinant();
		
		// shape function derivatives in global coordinates
		Eigen::MatrixXd dN_local(2, 9);
		dN_local << dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi, dN5_dxi, dN6_dxi, dN7_dxi, dN8_dxi, dN9_dxi, 
			dN1_deta, dN2_deta, dN3_deta, dN4_deta, dN5_deta, dN6_deta, dN7_deta, dN8_deta, dN9_deta;
		Eigen::MatrixXd dN_global = J_inv*dN_local;
		
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
	
	std::cout << "Assembling element " << std::flush;
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
			
			double Gc = (y9 > h2 && y9 < h2 + hi) ? Gc_eff : Gc_bulk;
			double E = (y9 > h2 && y9 < h2 + hi) ? E_int : E_bulk;
			double nu = (y9 > h2 && y9 < h2 + hi) ? nu_int : nu_bulk;
			
			Eigen::VectorXd s_el(9); s_el << s(n1), s(n2), s(n3), s(n4), s(n5), s(n6), s(n7), s(n8), s(n9);
			Eigen::VectorXd u_el(18); u_el << u(2*n1), u(2*n1 + 1), u(2*n2), u(2*n2 + 1), u(2*n3), u(2*n3 + 1), u(2*n4), u(2*n4 + 1), u(2*n5), u(2*n5 + 1), u(2*n6), u(2*n6 + 1), u(2*n7), u(2*n7 + 1), u(2*n8), u(2*n8 + 1), u(2*n9), u(2*n9 + 1);
			std::vector<double> H_el = H[i];
			
			Eigen::MatrixXd k_solid = QuadraticSolidElement(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8, x9, y9, E, nu, t, u_el, s_el, isPlaneStress, quadratureDegree); // solid element stiffness matrix
			AssembleQuadraticSolidMatrix(triplets, k_solid, n1, n2, n3, n4, n5, n6, n7, n8, n9);
			
			Eigen::MatrixXd k_pffm; // PFFM element stiffness matrix
			Eigen::VectorXd f_pffm; // PFFM element force vector
			QuadraticPFFMElement(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8, x9, y9, Gc, ls, s_el, H_el, quadratureDegree, k_pffm, f_pffm); // PFFM element stiffness matrix and force vector
			int offset = 2*nodes.size(); // assemble PFFM matrix and force vector to lower-right and bottom, respectively
			AssembleQuadraticPFFMMatrix(triplets, k_pffm, offset + n1, offset + n2, offset + n3, offset + n4, offset + n5, offset + n6, offset + n7, offset + n8, offset + n9);
			AssembleQuadraticPFFMVector(F, f_pffm, offset + n1, offset + n2, offset + n3, offset + n4, offset + n5, offset + n6, offset + n7, offset + n8, offset + n9);
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
			
			double yc = (y1 + y2 + y3 + y4)/4; // y coordinate of centre
			double Gc = (yc > h2 && yc < h2 + hi) ? Gc_eff : Gc_bulk;
			double E = (yc > h2 && yc < h2 + hi) ? E_int : E_bulk;
			double nu = (yc > h2 && yc < h2 + hi) ? nu_int : nu_bulk;
			
			Eigen::VectorXd s_el(4); s_el << s(n1), s(n2), s(n3), s(n4);
			Eigen::VectorXd u_el(8); u_el << u(2*n1), u(2*n1 + 1), u(2*n2), u(2*n2 + 1), u(2*n3), u(2*n3 + 1), u(2*n4), u(2*n4 + 1);
			std::vector<double> H_el = H[i];
			
			Eigen::MatrixXd k_solid = LinearSolidElement(x1, y1, x2, y2, x3, y3, x4, y4, E, nu, t, u_el, s_el, isPlaneStress, quadratureDegree); // solid element stiffness matrix
			AssembleLinearSolidMatrix(triplets, k_solid, n1, n2, n3, n4);
			
			Eigen::MatrixXd k_pffm; // PFFM element stiffness matrix
			Eigen::VectorXd f_pffm; // PFFM element force vector
			LinearPFFMElement(x1, y1, x2, y2, x3, y3, x4, y4, Gc, ls, s_el, H_el, quadratureDegree, k_pffm, f_pffm); // PFFM element stiffness matrix and force vector
			int offset = 2*nodes.size(); // assemble PFFM matrix and force vector to lower-right and bottom, respectively
			AssembleLinearPFFMMatrix(triplets, k_pffm, offset + n1, offset + n2, offset + n3, offset + n4);
			AssembleLinearPFFMVector(F, f_pffm, offset + n1, offset + n2, offset + n3, offset + n4);
		}
		
		if((i + 1) % 10000 == 0)
		{
			std::cout << i + 1 << ", ";
		}
		else if(i == elements.size() - 1)
		{
			std::cout << i + 1 << std::endl;
		}
	}
	K.setFromTriplets(triplets.begin(), triplets.end());
}

void UpdateH(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const Eigen::VectorXd& u, std::vector<std::vector<double>>& H, double E_int, double E_bulk, double nu_int, double nu_bulk, double h2, double hi, double ls, bool isPlaneStress, bool isQuadratic, int quadratureDegree, std::vector<std::vector<double>>&exx, std::vector<std::vector<double>>&eyy, std::vector<std::vector<double>>&exy, std::vector<std::vector<double>>&tr_e)
{
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
			
			double E = (y9 > h2 && y9 < h2 + hi) ? E_int : E_bulk;
			double nu = (y9 > h2 && y9 < h2 + hi) ? nu_int : nu_bulk;
			
			Eigen::VectorXd u_el(18); u_el << u(2*n1), u(2*n1 + 1), u(2*n2), u(2*n2 + 1), u(2*n3), u(2*n3 + 1), u(2*n4), u(2*n4 + 1), u(2*n5), u(2*n5 + 1), u(2*n6), u(2*n6 + 1), u(2*n7), u(2*n7 + 1), u(2*n8), u(2*n8 + 1), u(2*n9), u(2*n9 + 1);
			
			std::vector<double> H_old = H[i];
			std::vector<double> exx_i, eyy_i, exy_i, tr_e_i;
			
			std::vector<double> U_tensile = GetQuadraticElementTensileElasticStrainEnergy(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8, x9, y9, E, nu, u_el, isPlaneStress, quadratureDegree, exx_i, eyy_i, exy_i, tr_e_i);
			
			std::vector<double> H_new(H_old.size());
			std::transform(H_old.begin(), H_old.end(), U_tensile.begin(), H_new.begin(), [](double a, double b) { return std::max(a, b); }); // compute element-wise max
			H[i] = H_new;
			
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
			double x1 = nodes[n1][0];
			double y1 = nodes[n1][1];
			double x2 = nodes[n2][0];
			double y2 = nodes[n2][1];
			double x3 = nodes[n3][0];
			double y3 = nodes[n3][1];
			double x4 = nodes[n4][0];
			double y4 = nodes[n4][1];
			
			double yc = (y1 + y2 + y3 + y4)/4; // y coordinate of centre
			double E = (yc > h2 && yc < h2 + hi) ? E_int : E_bulk;
			double nu = (yc > h2 && yc < h2 + hi) ? nu_int : nu_bulk;
			
			Eigen::VectorXd u_el(8); u_el << u(2*n1), u(2*n1 + 1), u(2*n2), u(2*n2 + 1), u(2*n3), u(2*n3 + 1), u(2*n4), u(2*n4 + 1);
			
			std::vector<double> H_old = H[i];
			std::vector<double> U_tensile = GetLinearElementTensileElasticStrainEnergy(x1, y1, x2, y2, x3, y3, x4, y4, E, nu, u_el, isPlaneStress, quadratureDegree);
			
			std::vector<double> H_new(H_old.size());
			std::transform(H_old.begin(), H_old.end(), U_tensile.begin(), H_new.begin(), [](double a, double b) { return std::max(a, b); }); // compute element-wise max
			H[i] = H_new;
		}
		
		if((i + 1) % 10000 == 0)
		{
			std::cout << i + 1 << ", " << std::flush;
		}
		else if(i == elements.size() - 1)
		{
			std::cout << i + 1 << std::endl;
		}
	}
}

std::vector<double> GetLinearElementTensileElasticStrainEnergy(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double E, double nu, Eigen::VectorXd u, bool plane_stress, int quadratureDegree)
{
	// constitutive properties matrix
	Eigen::MatrixXd D(3, 3);
	double K;
	if(plane_stress == true)
	{
		D <<
			1., nu, 0.,
			nu, 1., 0.,
			0., 0., (1. - nu)/2.;
		D *= E/(1. - nu*nu);
		K = E/2./(1. - nu);
	}
	else
	{
		D <<
			1. - nu, nu, 0.,
			nu, 1. - nu, 0.,
			0., 0., (1. - 2.*nu)/2.;
		D *= E/(1. + nu)/(1. - 2.*nu);
		K = E/2./(1. + nu)/(1. - 2.*nu);
	}
	
	// integration points and weightings
	Eigen::VectorXd xi_points;
	Eigen::VectorXd eta_points;
	Eigen::VectorXd weights;
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	
	std::vector<double> H(weights.size());
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
		double epsilon_vol = epsilon(0) + epsilon(1); // volumetric strain
		double sigma_vol = (sigma(0) + sigma(1))/2.; // volumetric stress (positive in tension)
		bool tensile = (sigma_vol > -std::numeric_limits<double>::epsilon()); // sigma_vol >= 0 with some tolerance for floating point error
		double U = 0.5*sigma.transpose()*epsilon;
		double U_vol = 0.5*K*epsilon_vol*epsilon_vol;
		double U_dev = U - U_vol;
		
		H[p] = U_dev;
	}
	return H;
}

std::vector<double> GetQuadraticElementTensileElasticStrainEnergy(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double x5, double y5, double x6, double y6, double x7, double y7, double x8, double y8, double x9, double y9, double E, double nu, Eigen::VectorXd u, bool plane_stress, int quadratureDegree, std::vector<double>& exx_i, std::vector<double>& eyy_i, std::vector<double>& exy_i, std::vector<double>& tr_e_i)
{
	// constitutive properties matrix
	Eigen::MatrixXd D(3, 3);
	double K;
	if(plane_stress == true)
	{
		D <<
			1., nu, 0.,
			nu, 1., 0.,
			0., 0., (1. - nu)/2.;
		D *= E/(1. - nu*nu);
		K = E/2./(1. - nu);
	}
	else
	{
		D <<
			1. - nu, nu, 0.,
			nu, 1. - nu, 0.,
			0., 0., (1. - 2.*nu)/2.;
		D *= E/(1. + nu)/(1. - 2.*nu);
		K = E/2./(1. + nu)/(1. - 2.*nu);
	}
	
	// integration points and weightings
	Eigen::VectorXd xi_points;
	Eigen::VectorXd eta_points;
	Eigen::VectorXd weights;
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	
	std::vector<double> H(weights.size());
	exx_i.resize(weights.size());
	eyy_i.resize(weights.size());
	exy_i.resize(weights.size());
	tr_e_i.resize(weights.size());
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
		double epsilon_vol = epsilon(0) + epsilon(1); // volumetric strain
		double sigma_vol = (sigma(0) + sigma(1))/2.; // volumetric stress (positive in tension)
		bool tensile = (sigma_vol > -std::numeric_limits<double>::epsilon()); // sigma_vol >= 0 with some tolerance for floating point error
		double U = 0.5*sigma.transpose()*epsilon;
		double U_vol = 0.5*K*epsilon_vol*epsilon_vol;
		double U_dev = U - U_vol;
		if(tensile)
		{
			H[p] = U;
		}
		else
		{
			H[p]=U_dev;
		}
		exx_i[p] = epsilon(0);
		eyy_i[p] = epsilon(1);
		exy_i[p] = epsilon(2);
		tr_e_i[p] = epsilon_vol;
	}
	return H;
}

void add_crack(std::vector<std::vector<double>>& nodes, std::vector<std::vector<int>>& elements, double xct, double yc, double tol)
{
	std::map<int, int> nodeDupMap;

	// find and duplicate nodes along crack line
	int numNodes_initial = nodes.size();
	for(int i = 0; i < numNodes_initial; ++i)
	{
		double x = nodes[i][0];
		double y = nodes[i][1];

		if(x >= xct - tol && std::abs(y - yc) < tol)
		{
			nodes.push_back(nodes[i]);
			int newIndex = nodes.size() - 1;
			nodeDupMap[i] = newIndex;
		}
	}

	// update elements on one side of the crack
	for(auto& elem : elements)
	{
		// calculate the average y-coordinate of the element
		double avg_y = 0.0;
		for(int idx : elem)
		{
			avg_y += nodes[idx][1];
		}
		avg_y /= elem.size();
		
		// if the average y-coordinate is below the crack line, no need to update the element
		if(avg_y < yc)
		{
			continue; // this element is below the crack line and remains unchanged
		}
		
		// otherwise, the element is above the crack line and needs updating
		for(int& idx : elem)
		{
			double x = nodes[idx][0];
			double y = nodes[idx][1];

			if(x >= xct - tol && std::abs(y - yc) < tol)
			{
				if(nodeDupMap.count(idx))
				{
					idx = nodeDupMap[idx];
				}
			}
		}
	}
}

double GetGcInitiationFromComplianceCalibration(const std::vector<double>& D, const std::vector<double>& F1, const std::vector<double>& F2, const std::vector<double>& a, double a_min)
{
	size_t size = 0;
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
	double Gc_init = 0.5*P_init*P_init*dC_da_init; // Gc at initiation
	return Gc_init;
}

void SaveHField(const std::vector<std::vector<double>>& nodes, const std::vector<std::vector<int>>& elements, const std::vector<std::vector<double>>& H, bool isQuadratic, int quadratureDegree, std::string filename, int prec, const std::vector<std::vector<double>>& exx, const std::vector<std::vector<double>>& eyy, const std::vector<std::vector<double>>& exy, const std::vector<std::vector<double>>& tr_e)
{
	Eigen::VectorXd xi_points;
	Eigen::VectorXd eta_points;
	Eigen::VectorXd weights;
	Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
	
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
				double Hval = H[i][p]; // H[i][j] = H for element i at integration point j
				double exxval = exx[i][p]; // H[i][j] = H for element i at integration point j
				double eyyval = eyy[i][p]; // H[i][j] = H for element i at integration point j
				double exyval = exy[i][p]; // H[i][j] = H for element i at integration point j
				double tr_eval = tr_e[i][p]; // H[i][j] = H for element i at integration point j
				
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
				double Hval = H[i][p]; // H[i][j] = H for element i at integration point j
				
				outfile << xp
					<< "," << yp
					<< "," << Hval
					<< std::endl;
			}
		}
	}
	outfile.close();
}

// remove constrained rows/cols from K and update F for nonzero prescribed displacements
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
