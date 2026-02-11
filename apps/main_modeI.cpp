#include "modeI.h"
#include "inputhandling.h"
#include <iostream>

int main(int argc, char** argv) {

	ModeIInputData data;

	data.L = 0.1; // length of DCB (m)
	data.h1 = 2e-3; // height of upper DCB arm (m)
	data.h2 = 2e-3; // height of lower DCB arm (m)
	data.hi = 1e-4; // height of interface (m)
	data.a0 = 0.03; // initial crack length (m)
	data.ls = 20e-6; // length scale
	data.t = 1; // thickness
	data.a_min = data.a0 + 1e-3; // minimum crack length before dPi/dA fitting is started
	
	// material properties
	data.E_int = 126.e9; // interface Young's modulus
	data.nu_int = 0.3;  // interface Poisson's ratio
	data.Gc_int = 281.; // interface fracture toughness
	data.Gc_eff = 252.96; // effective interface fracture toughness for FEM only
	data.E_bulk = 126.e9; // bulk Young's modulus
	data.nu_bulk = 0.3; // bulk Poisson's ratio
	data.Gc_bulk = 10.*data.Gc_int;  // bulk fracture toughness
	data.isPlaneStress = true; // if false, plane strain
		
	// Loading
	data.Dinit = 0.5875e-3 ; // initial displacement
	data.Dend = 100e-3; // applied displacement
	data.Dinc = 6.25e-9; // initial displacement increment

	// Run parameters
	data.tol = 1e-7; // Newton's method convergence tolerance
	data.relax = 1.; // fraction of the increment to apply each iteration
	data.save_freq = 500; // number of increments per files getting saved (1 = every increment, 2 = even increments, etc.)
	data.restart_from = 0; // iteration number to restart from

	// Mesh properties
	data.Delta_y = data.L/4000; // initial element length at crack tip
	data.Delta_y_max_ahead = data.L/10.; // max element length ahead of crack tip
	data.Delta_y_max_behind = 1e-4; // max element length behind crack tip
	data.GR_y_ahead = 1.1; // growth rate of element length ahead of crack tip
	data.GR_y_behind = 1.2; // growth rate of element length behind crack tip
	data.Delta = data.hi/24; // interface element height
	data.Delta_x = data.Delta; // initial element height in bulk
	data.Delta_x_max = data.h1/10.; // max element height in bulk
	data.yminus = 2e-3; // length of refined region ahead of crack tip
	data.yplus = 1e-3; // length of refined region behind crack tip
	data.xplus = data.hi/2; // length of refined region above interface
	data.xminus = data.hi/2; // length of refined region below interface
	data.GR_x = 1.5; // growth rate of element height in bulk
	data.isQuadratic = true; // if false, use linear elements
	data.quadratureDegree = 5; // degree of polynomial that can be fully integrated


	if (argc >= 3) {
        // RESTART: Load from file
        std::string restart_dir = argv[1];
        int restart_inc = std::stoi(argv[2]);
        
        load_parameters(restart_dir + "/parameters.txt", data);
        data.restart_from = restart_inc;

    } else {
        // FRESH: Save to file
        if(argc >= 2) {
            std::string output_dir = argv[1];
            save_parameters(output_dir + "/parameters.txt", data);
        }
    }

    ModeISimulation sim(data);

    try{
        sim.Run();
    }
    catch (const std::exception& e) {
        std::cerr << "Simulation failed:" << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}