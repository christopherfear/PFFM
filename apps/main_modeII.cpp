// =========================================================================
// MODE II DRIVER (end-loaded split)
// =========================================================================

// built in modules
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <unordered_map>

// custom models
#include "mesh.h"
#include "elements.h"
#include "assembly.h"
#include "postprocess.h"
#include "linearsolver.h"

// external libraries
#include <Eigen/Core>
#include <Eigen/Sparse>

int main()
{
	// ---------------------------------------------------
	// 1. PARAMETERS (Specific to Mode II ELS)
	// ---------------------------------------------------

	// dimensions
	DCBMeshParameters mp; // Container for specimen geometry and mesh refinement settings

	mp.L = 0.1; // length of ELS (m)
	mp.h1 = 2e-3; // height of upper ELS arm (m)
	mp.h2 = 2e-3; // height of lower ELS arm (m)
	mp.hi = 1e-4; // height of interface (m)
	mp.a0 = 0.03; // initial crack length (m)
	double ls = 20e-6; // length scale
	double t = 1; // thickness
	double a_min = mp.a0 + 3.7e-5; // minimum crack length before compliance calibration fitting is started
	
	// material properties
	double E_int = 126.e9; // interface Young's modulus
	double nu_int = 0.3; // interface Poisson's ratio
	double Gc_int = 562.; // interface fracture toughness
	double Gc_eff = 497.63; // effective interface fracture toughness for FEM only
	double E_bulk = 126.e9; // bulk Young's modulus
	double nu_bulk = 0.3; // bulk Poisson's ratio
	double Gc_bulk = 10.*Gc_int; // bulk fracture toughness
	bool isPlaneStress = true; // if false, plane strain
	
	// mesh properties
	mp.Delta_y = mp.L/3000; // initial element length at crack tip
	mp.Delta_y_max_ahead = mp.L/10.; // max element length ahead of crack tip
	mp.Delta_y_max_behind = 1e-4; // max element length behind crack tip
	mp.GR_y_ahead = 1.1; // growth rate of element length ahead of crack tip
	mp.GR_y_behind = 1.2; // growth rate of element length behind crack tip
	mp.Delta = mp.hi/24; // interface element height
	mp.Delta_x = mp.Delta; // initial element height in bulk
	mp.Delta_x_max = mp.h1/10.; // max element height in bulk
	mp.yminus = 4e-3; // length of refined region ahead of crack tip
	mp.yplus = 1e-3; // length of refined region behind crack tip
	mp.xplus = mp.hi/2; // length of refined region above interface
	mp.xminus = mp.hi/2; // length of refined region below interface
	mp.GR_x = 1.5; // growth rate of element height in bulk
	mp.isQuadratic = true; // if false, use linear elements
	int quadratureDegree = 5; // degree of polynomial that can be fully integrated
	
	// Run parameters
	double tol = 1e-7; // Newton's method convergence tolerance
	double relax = 1.; // fraction of the increment to apply each iteration
	double s_threshold = 0.01; // damage value for crack tip definition 
	int save_freq = 500; // number of increments per files getting saved (1 = every increment, 2 = even increments, etc.)
	int restart_from = 0; // iteration number to restart from
	
	
	// Output files
	std::string Hfilename = "H,inc=";
	std::string usfilename = "inc=";
	std::string history_filename = "history.csv";
	int prec = 16; 
	
	// Loading 
	double Dinit = 8.48e-3 ; // initial displacement
	double Dend = 100e-3; // applied displacement
	double Dinc = 2.5e-9; // initial displacement increment
	
	
	// ---------------------------------------------------
    // 2. MESH GENERATION
    // ---------------------------------------------------

    // Create containers and generate the mesh
    std::vector<std::vector<double>> nodes;
    std::vector<std::vector<int>> elements;
    GenerateDCBMesh(mp, nodes, elements);
	
	int numElements = elements.size();
	int numNodes = nodes.size();
	std::cout << "Number of elements: " << numElements << std::endl;
	std::cout << "Number of nodes: " << numNodes << std::endl << std::endl;
	
	// Save Geometry
	std::ofstream outfile("geometry.csv");
	outfile << "L,h1,h2,hi,a0" << std::endl;
	outfile << std::fixed << std::setprecision(prec) << mp.L
		<< "," << mp.h1
		<< "," << mp.h2
		<< "," << mp.hi
		<< "," << mp.a0
		<< std::endl;
	outfile.close();
	
	// ---------------------------------------------------
    // 3. NODE SETS
    // ---------------------------------------------------

    NodeSets sets; // struct from mesh.h
	double findtol = mp.Delta/10;
	for(size_t i = 0; i < nodes.size(); ++i)
	{
		double x = nodes[i][0];
		double y = nodes[i][1];
		
		// set of nodes forming the initial crack
		if((x >= (mp.L - mp.a0) - findtol) && (std::abs(y - (mp.h2 + mp.hi/2)) < findtol))
		{
			sets.initial_crack_nodes.push_back(i);
		}
		
		// set of nodes forming the centreline of the interface
		if((x <= (mp.L - mp.a0) + findtol) && (std::abs(y - (mp.h2 + mp.hi/2)) < findtol))
		{
			sets.centreline_nodes.push_back(i);
		}
		
		// set of all nodes forming the interface
		if((x <= (mp.L - mp.a0) + findtol) && (y >= mp.h2 - findtol) && (y <= mp.h2 + mp.hi + findtol))
		{
			sets.interface_nodes.push_back(i);
		}
		
		// set of nodes on end of upper arm
		if((y >= (mp.h2 + mp.hi) - findtol) && (std::abs(x - mp.L) < findtol))
		{
			sets.upper_arm_end_nodes.push_back(i);
		}
		
		// set of nodes on end of lower arm
		if((y <= mp.h2 + findtol) && (std::abs(x - mp.L) < findtol))
		{
			sets.lower_arm_end_nodes.push_back(i);
		}
		
		// set of nodes on support end
		if((std::abs(x) < findtol))
		{
			sets.support_nodes.push_back(i);
		}
	}


	// ---------------------------------------------------
	// 4. INITIALISATION & MEMORY ALLOCATION
	// ---------------------------------------------------

	Eigen::VectorXd u = Eigen::VectorXd::Zero(2*numNodes);
	Eigen::VectorXd s = Eigen::VectorXd::Ones(numNodes);
	std::vector<std::vector<double>> H, exx, eyy, exy, tr_e;
	std::vector<double> F1_history, D1_history, F2_history, D2_history, a_history, Gc_CC_history, U_history;
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
                // else if(i == 6) GcInit_CC_history.push_back(val); // don't store logged value - recalculate below based on current a_min
                else if(i == 7) U_history.push_back(val);
                
                // recalculate Gc_initiation from compliance calibration based on history to this point and value of a_min
                if(i == 5)
                {
                    Gc_CC_history.push_back(GetGcfromCC(D1_history, F1_history, F2_history, a_history, a_min, t));
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
        // initial crack damage
		for(auto& node : sets.initial_crack_nodes)
        {
            s[node] = 0;
        }
        
		// fresh allocation for field variables
        Eigen::VectorXd xi_points, eta_points, weights; // integration points and weightings
        Get2DQuadrature(quadratureDegree, xi_points, eta_points, weights);
        for(size_t i = 0; i < numElements; ++i)
        {
			// vector of zeros to allocate memeory
        	std::vector<double> zeros(weights.size(), 0.);
			H.push_back(zeros);
			exx.push_back(zeros); eyy.push_back(zeros); exy.push_back(zeros); tr_e.push_back(zeros);
        }
        
    }

	// ---------------------------------------------------
	// 5. SOLVER LOOP & BOUNDARY CONDITIONS
	// ---------------------------------------------------
	LinearSolver solver;
	static bool firstLoop = true;

	Eigen::SparseMatrix<double> K(3*numNodes, 3*numNodes);
	Eigen::VectorXd F(K.rows());
	Eigen::VectorXd d(K.rows());
	Eigen::VectorXd unm1(u.size());
	Eigen::VectorXd Dd(K.rows());

	// BCs Setup
	std::vector<int> constrained;
	for(int node : sets.initial_crack_nodes) constrained.push_back(2*numNodes + node);
	for(int node : sets.support_nodes) { constrained.push_back(2*node); constrained.push_back(2*node + 1); }
	for(int node : sets.lower_arm_end_nodes) constrained.push_back(2*node + 1);
	for(int node : sets.upper_arm_end_nodes) constrained.push_back(2*node + 1);
	
	std::vector<double> constrainedValues(constrained.size(), 0.);
	std::unordered_map<int, double*> constrainedMap;
	constrainedMap.reserve(constrained.size());
	for(size_t k = 0; k < constrained.size(); ++k) constrainedMap[constrained[k]] = &constrainedValues[k];
	
	std::vector<bool> isConstrained(K.rows(), false);
	for(const auto &pair : constrainedMap) isConstrained[pair.first] = true;
	
	std::vector<int> freeMap(K.rows(), -1);
	int nFree = 0;
	for(int i = 0; i < K.rows(); ++i) if(!isConstrained[i]) freeMap[i] = nFree++;
	
	Eigen::SparseMatrix<double> KR(nFree, nFree);
	Eigen::VectorXd d_reduced(nFree), FR(nFree), Dd_reduced(nFree), residual(nFree);

	// Incremental Loading
	for(double D = Dinit; D <= Dend + Dinc/10; D += Dinc)
	{
		auto t1 = std::chrono::high_resolution_clock::now();
		printf("Increment number: %d\n", incr_count);
		printf(("Applied displacement: %." + std::to_string(prec) + "f\n").c_str(), D);
		
		// Update BCs: Apply displacement D to the loading ends
        for(int node : sets.upper_arm_end_nodes) {
            *constrainedMap[2*node + 1] = D;
        }
        for(int node : sets.lower_arm_end_nodes) {
            *constrainedMap[2*node + 1] = D;
        }
		
		// combine u and s into combined solution vector
		d.head(u.size()) = u; // assign u to first part of d
		d.tail(s.size()) = s; // assign s to second part of d
		unm1 = u; 
		int iter_count = 1;

		while(true)
		{
			printf("Iteration: %d\n", iter_count);
			
			// Call Assembly
			AssembleAll(nodes, elements, u, s, H, E_int, E_bulk, nu_int, nu_bulk, Gc_eff, Gc_bulk, t, mp.h2, mp.hi, ls, K, F, isPlaneStress, mp.isQuadratic, quadratureDegree);

			// Call BCs
			ApplyDirichletBC(K, F, isConstrained, constrainedMap, freeMap, KR, FR);
			
			for(int i = 0; i < K.rows(); ++i) {
				int idx = freeMap[i];
				if(idx != -1) d_reduced[idx] = d[i];
			}
			
			std::cout << "Solving... " << std::flush;
            residual = FR - KR*d_reduced;
            
            // 1. Analyze Pattern (Only on first iteration or if mesh changes)
            if(firstLoop) {
                solver.analyzePattern(KR);
                firstLoop = false;
            }
            
            // 2. Solve System
            // The LinearSolver class handles the waterfall (Cholmod -> LLT -> LDLT -> LU)
            bool success = solver.solve(KR, residual, Dd_reduced);
            
            if(!success) {
                std::cerr << "Aborting simulation due to solver failure." << std::endl;
                return -1;
            }
			
			for(int i = 0; i < K.rows(); ++i) {
				int idx = freeMap[i];
				if (idx != -1) Dd[i] = Dd_reduced[idx];
				else Dd[i] = 0.0;
			}

			d += relax*Dd;
			for(const auto& pair : constrainedMap) d[pair.first] = *(pair.second);
			u = d.head(u.size());
			s = d.tail(s.size());
			std::cout << "Done! ";
			
			double Rk = residual.norm();
			double R_normalised = Rk/FR.norm();
			printf("|R|/|F| = %.3e\n", R_normalised);
			
			if(R_normalised < tol)
			{
				// Post Processing
				F = K*d;
				double F1 = 0.0; for(int node : sets.upper_arm_end_nodes) F1 += F(2*node + 1);
				F1_history.push_back(F1); D1_history.push_back(D);
				
				double F2 = 0.0; for(int node : sets.lower_arm_end_nodes) F2 += F(2*node + 1);
				F2_history.push_back(F2); D2_history.push_back(D);

				// Crack Length Calculation
				double crack_length = CalculateCrackLength(nodes, sets.interface_nodes, d, mp.L, s_threshold, mp.Delta);
				a_history.push_back(crack_length);

				// Strain Energy
				double U = 0.5*u.dot(K.topLeftCorner(u.size(), u.size())*u);
				U_history.push_back(U);
				
				// Gc Calculation (Mode II Specific)
				double Gc_CC = GetGcfromCC(D1_history, F1_history, F2_history, a_history, a_min, t);
				Gc_CC_history.push_back(Gc_CC);
				std::cout << "Gc from CC = " << Gc_CC << std::endl;

				// Save History
				std::cout << "Saving history.csv" << std::endl;
				std::ofstream outfile(history_filename);
				outfile << "inc,D1,F1,D2,F2,a,G_c,U" << std::endl;
				outfile << std::fixed << std::setprecision(prec);
				for(size_t i = 0; i < F1_history.size(); ++i) {
					outfile << i << "," << D1_history[i] << "," << F1_history[i] << "," << D2_history[i] << "," 
					        << F2_history[i] << "," << a_history[i] << "," << Gc_CC_history[i] << "," 
							<< U_history[i] << std::endl;
				}
				outfile.close();
				break;
			}
			++iter_count;
		}
		
		// Update H (Calls Assembly.cpp)
		UpdateH(nodes, elements, u, H, E_int, E_bulk, nu_int, nu_bulk, mp.h2, mp.hi, ls, isPlaneStress, mp.isQuadratic, quadratureDegree, exx, eyy, exy, tr_e);
		
		if(incr_count % save_freq == 0) {
			std::string filename = Hfilename + std::to_string(incr_count) + ".csv";
			SaveHField(nodes, elements, H, mp.isQuadratic, quadratureDegree, filename, prec, exx, eyy, exy, tr_e);
			
			filename = usfilename + std::to_string(incr_count) + ".csv";
			std::ofstream usfile(filename);
			usfile << "x,y,u,v,s" << std::endl;
			for(size_t i = 0; i < nodes.size(); ++i) {
				usfile << std::fixed << std::setprecision(prec) << nodes[i][0] << "," << nodes[i][1] << "," 
				       << u[2*i] << "," << u[2*i + 1] << "," << s[i] << std::endl;
			}
			usfile.close();
		}
		
		auto t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> seconds = t2 - t1;
		printf("Increment %d took %.3f seconds\n", incr_count, seconds.count());
		std::cout << std::endl;
		
		++incr_count;
	}

	return 0;
}