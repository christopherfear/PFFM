#ifndef MODEI_H
#define MODEI_H

#include "simulation.h"

class ModeISimulation : public simulation { // inherit everything from simulation class

public:
    ModeISimulation(){
        std::cout << "Initialising mode-I DCB simulation" << std::endl;
    }

protected:

    void Setup() override {

        // dimensions
        mp.L = 0.1; // length of DCB (m)
        mp.h1 = 2e-3; // height of upper DCB arm (m)
        mp.h2 = 2e-3; // height of lower DCB arm (m)
        mp.hi = 1e-4; // height of interface (m)
        mp.a0 = 0.03; // initial crack length (m)

        ls = 20e-6; // length scale (m)
        t = 1.; // thickness
        a_min = mp.a0 + 1e-3; // minimum crack length for Gc calculation

        // mesh parameters
        mp.Delta_y = mp.L/4000; // initial element length at crack tip
        mp.Delta_y_max_ahead = mp.L/10.; // max element length ahead of crack tip
        mp.Delta_y_max_behind = 1e-4; // max element length behind crack tip
        mp.GR_y_ahead = 1.1; // growth rate of element length ahead of crack tip
        mp.GR_y_behind = 1.2; // growth rate of element length behind crack tip
        mp.Delta = mp.hi/24; // interface element height
        mp.Delta_x = mp.Delta; // initial element height in bulk
        mp.Delta_x_max = mp.h1/10.; // max element height in bulk
        mp.yminus = 2e-3; // length of refined region ahead of crack tip
        mp.yplus = 1e-3; // length of refined region behind crack tip
        mp.xplus = mp.hi/2; // length of refined region above interface
        mp.xminus = mp.hi/2; // length of refined region below interface
        mp.GR_x = 1.5; // growth rate of element height in bulk

        // material properties
        E_int = 126e9; // Young's modulus (interface)
        E_bulk = E_int; // Young's modulus (bulk)
        nu_int = 0.3; // Poisson's ratio (interface)
        nu_bulk = nu_int; // Poisson's ratio (bulk)

        Gc_int = 281.; // Fracture toughness (interface)
        Gc_bulk = 10.*Gc_int; // Fracture toughness (bulk)
        Gc_eff = 252.96; // Effective fracture toughness

        // run controls
        D = 0.6875e-3;        // Dinit (Starting displacement)
        Dinc = 6.25e-9;       // Step size
        Dend = 100e-3;        // End displacement
        tol = 1e-7;           // Tolerance
        save_freq = 500;      // Save full results every 500 steps
        restart_from = 0;     // Set to > 0 to resume a crash
        relax = 1.0;          // Relaxation
        s_threshold = 0.01;   // Crack tip definition
        prec = 16;            // CSV Decimal precision

        isPlaneStress = true;
        isQuadratic = true;
        quadratureDegree = 5;

        GenerateDCBMesh(mp, nodes, elements);

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
    }

    void ApplyBoundaryConditions() override {
        
        for(int node : sets.upper_arm_end_nodes){*constrainedMap[2*node+1] = D;}
        for(int node : sets.lower_arm_end_nodes){*constrainedMap[2*node+1] = -D;}

        
    }

    double CalculateGcForRestart() override {

        return GetGcFromdPidA(D1_history, F1_history, a_history, U_history, a_min);
    }

    void SaveResults() override {
        // Post processomg
        F = K*d;
        double F1 = 0.0; for(int node : sets.upper_arm_end_nodes) F1 += F(2*node + 1);
        F1_history.push_back(F1); D1_history.push_back(D);
				
        double F2 = 0.0; for(int node : sets.lower_arm_end_nodes) F2 += F(2*node + 1);
        F2_history.push_back(F2); D2_history.push_back(D);

        // Crack length calculation
        double crack_length = CalculateCrackLength(nodes, sets.interface_nodes, d, mp.L, s_threshold, mp.Delta);
        a_history.push_back(crack_length);

        // Strain energy
        double U = 0.5*u.dot(K.topLeftCorner(u.size(), u.size())*u);
        U_history.push_back(U);

        // Gc calculation specific to mode-I
        double current_Gc = GetGcFromdPidA(D1_history, F1_history, a_history, U_history, a_min);
        Gc_history.push_back(current_Gc);
        std::cout << "Gc from -dPi/dA = " << current_Gc << std::endl;

        // history.csv save
        std::cout << "Saving history.csv" << std::endl;
        history_file << incr_count << "," << D << "," << F1 << "," << D << "," << F2 << "," 
                     << crack_length << "," << current_Gc << "," << U << std::endl;

        // periodic full-field save
        if (incr_count % save_freq == 0) {
            std::string filename = usfilename + std::to_string(incr_count) + ".csv";
            std::ofstream usfile(filename);
            usfile << "x,y,u,v,s\n";
            for(size_t i = 0; i < nodes.size(); ++i) {
                usfile << std::fixed << std::setprecision(prec) << nodes[i][0] << "," << nodes[i][1] << "," 
                       << u[2*i] << "," << u[2*i+1] << "," << s[i] << "\n";
            }
            usfile.close();

            std::string h_out = Hfilename + std::to_string(incr_count) + ".csv";
            SaveHField(nodes, elements, H, isQuadratic, quadratureDegree, h_out, prec, exx, eyy, exy, tr_e);
        }
         
    }



};

#endif