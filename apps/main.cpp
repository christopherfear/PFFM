#include <iostream>
#include <memory>
#include <string>

#include "inputhandling.h" 
#include "modeI.h"
#include "modeII.h"

int main(int argc, char** argv) {

    InputData data;

    // Placeholder data aligning with mode-I simulation
    // Used for fresh runs, unless overridden by a parameters results file
    
    data.L = 0.1; // Length of beams
    data.h1 = 2e-3; // Height of upper arm
    data.h2 = 2e-3; // Height of lower arm
    data.hi = 1e-4; // Height of interface
    data.a0 = 0.03; // Initial crack length

    data.ls = 20e-6; // length scale
    data.t = 1.0; // thickness
    data.a_min = data.a0 + 1e-3; // minimum crack length for fracture toughness calculation

    data.E_int = 126.e9; // Young's modulus (interface)
    data.nu_int = 0.3; // Poisson's ratio (interface)
    data.Gc_int = 281.0; // Fracture toughness (interface)
    data.Gc_eff = 252.96; // Effective fracture toughness
    data.E_bulk = 126.e9; // Young's modulus (bulk)
    data.nu_bulk = 0.3; // Poisson's ratio (bulk)
    data.Gc_bulk = 10. * data.Gc_int; // Fracture toughness (bulk)
    data.isPlaneStress = true; // if false, plane strain

    data.Dinit = 0.5875e-3; // Initial applied displacement
    data.Dend = 100e-3; // Final applied displacement
    data.Dinc = 6.25e-9; // Displacement increment

    data.tol = 1e-7; // tolerance
    data.relax = 1.0; // fraction of the incremement applied at each iteration
    data.save_freq = 2; // number of increments field data is saved 
    data.restart_from = 0; // iteration to restart from

    data.Delta_y = data.L / 4000.0; // initial element length at crack tip
    data.Delta_y_max_ahead = data.L / 10.0; // max element length ahead of crack tip
    data.Delta_y_max_behind = 1e-4; // max element length behind crack tip
    data.GR_y_ahead = 1.1; // growth rate of element length ahead of crack tip
    data.GR_y_behind = 1.2; // growth rate of element length behind crack tip
    data.Delta = data.hi / 24.0; // element height in interface
    data.Delta_x = data.Delta; // initial element height in bulk
    data.Delta_x_max = data.h1 / 10.0; // max element height in bulk
    data.yminus = 2e-3; // length of refined region ahead of crack tip
    data.yplus = 1e-3; // length of refined region behind crack tip
    data.xplus = data.hi / 2.0; // length of refined region above interface
    data.xminus = data.hi / 2.0; // length of refined region below interface
    data.GR_x = 1.5; // growth rate of element height in the bulk
    data.isQuadratic = true; // if false, use linear elements
    data.quadratureDegree = 5; // degree of polynomial that can be fully integrated

    std::string mode = "I";
    std::string restart_dir = "";
    
    if (argc >= 2) mode = argv[1];
    if (argc >= 4) { 
        restart_dir = argv[2];
        load_parameters(restart_dir + "/parameters.txt", data);
        data.restart_from = std::stoi(argv[3]);
    } else {
        if(argc >= 3) {
            std::string output_dir = argv[2];
            save_parameters(output_dir + "/parameters.txt", data);
        }
    }

    std::unique_ptr<simulation> sim = nullptr;

    if (mode == "I") {
        std::cout << ">> Selecting mode-I (DCB) solver..." << std::endl;
        sim = std::make_unique<ModeISimulation>(data, "I");
    } 
    else if (mode == "II") {
        std::cout << ">> Selecting mode-II (ELS) solver..." << std::endl;
        sim = std::make_unique<ModeIISimulation>(data, "II");
    } 
    else {
        std::cerr << "Error: Unknown mode '" << mode << "'" << std::endl;
        return 1;
    }

    try {
        if (sim) sim->Run(mode);
    } 
    catch (const std::exception& e) {
        std::cerr << "Simulation failed: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}