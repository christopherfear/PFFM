#ifndef MODEI_H
#define MODEI_H

#include "simulation.h"

struct ModeIInputData {
    // Geometry
    double L, h1, h2, hi, a0;
    
    // Physics & Materials
    double ls, t, a_min;
    double E_int, nu_int, Gc_int, Gc_eff;
    double E_bulk, nu_bulk, Gc_bulk;
    bool isPlaneStress;

    // Mesh Refinement
    double Delta_y, Delta_y_max_ahead, Delta_y_max_behind;
    double GR_y_ahead, GR_y_behind;
    double Delta, Delta_x, Delta_x_max;
    double yminus, yplus, xplus, xminus;
    double GR_x;
    bool isQuadratic;
    int quadratureDegree;

    // Solver Controls
    double tol, relax;
    int save_freq, restart_from;
    
    // Loading
    double Dinit, Dend, Dinc;
};


class ModeISimulation : public simulation { // inherit everything from simulation class

    ModeIInputData input; //store the inputs

public:
    ModeISimulation(const ModeIInputData data_in) : input(data_in) {
        std::cout << "Initialising mode-I DCB simulation with custom inputs." << std::endl;
    }

protected:

    void Setup() override {

        // Geometry
        mp.L = input.L; mp.h1 = input.h1; mp.h2 = input.h2; 
        mp.hi = input.hi; mp.a0 = input.a0;

        // Mesh Refinement
        mp.Delta_y = input.Delta_y;
        mp.Delta_y_max_ahead = input.Delta_y_max_ahead;
        mp.Delta_y_max_behind = input.Delta_y_max_behind;
        mp.GR_y_ahead = input.GR_y_ahead;
        mp.GR_y_behind = input.GR_y_behind;
        mp.Delta = input.Delta;
        mp.Delta_x = input.Delta_x;
        mp.Delta_x_max = input.Delta_x_max;
        mp.yminus = input.yminus;
        mp.yplus = input.yplus;
        mp.xplus = input.xplus;
        mp.xminus = input.xminus;
        mp.GR_x = input.GR_x;
        mp.isQuadratic = input.isQuadratic;
        
        // Materials & Physics
        ls = input.ls;
        t = input.t;
        a_min = input.a_min;
        E_int = input.E_int; nu_int = input.nu_int;
        E_bulk = input.E_bulk; nu_bulk = input.nu_bulk;
        Gc_int = input.Gc_int; Gc_bulk = input.Gc_bulk; Gc_eff = input.Gc_eff;
        isPlaneStress = input.isPlaneStress;

        // Solver
        tol = input.tol;
        relax = input.relax;
        save_freq = input.save_freq;
        restart_from = input.restart_from;
        quadratureDegree = input.quadratureDegree;
        
        // Loading
        D = input.Dinit;
        Dend = input.Dend;
        Dinc = input.Dinc;

        s_threshold = 0.01;
        prec = 16;

        GenerateDCBMesh(mp, nodes, elements);
        IdentifyNodeSets();
        SaveGeometryData();
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
            SaveHField(nodes, elements, H, mp.isQuadratic, quadratureDegree, h_out, prec, exx, eyy, exy, tr_e);
        }
         
    }



};

#endif