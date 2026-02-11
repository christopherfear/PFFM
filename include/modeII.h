#ifndef MODEII_H
#define MODEII_H

#include "simulation.h"
#include "inputhandling.h"

class ModeIISimulation : public simulation { // inherit everything from simulation class

    InputData input; //store the inputs

public:
    ModeIISimulation(const InputData& data_in, std::string mode) : input(data_in) {
    this->mode_name = mode;
}

protected:

    void Setup() override {

        LoadParameters(input);
        GenerateDCBMesh(mp, nodes, elements);
        IdentifyNodeSets();
        SaveGeometryData();
    
    }

    void ApplyBoundaryConditions() override {
        
        for(int node : sets.upper_arm_end_nodes){*constrainedMap[2*node+1] = D;}
        for(int node : sets.lower_arm_end_nodes){*constrainedMap[2*node+1] = D;}

    }

    double CalculateGcForRestart() override {

        return GetGcfromCC(D1_history, F1_history, F2_history, a_history, a_min, t);
    }

    void ExtractResults() override {
        // Post processing
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

        // Gc calculation specific to mode-II
        double current_Gc = GetGcfromCC(D1_history, F1_history, F2_history, a_history, a_min, t);
        Gc_history.push_back(current_Gc);
        std::cout << "Gc from CC = " << current_Gc << std::endl;

        // history.csv save
        std::cout << "Saving history.csv" << std::endl;
        history_file << incr_count << "," << D << "," << F1 << "," << D << "," << F2 << "," 
                     << crack_length << "," << current_Gc << "," << U << std::endl;
         
    }

};

#endif