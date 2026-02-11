#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include <chrono>
#include <iomanip> 
#include <unordered_map>
#include <sstream>
#include <Eigen/Core>
#include <Eigen/Sparse>


#include "mesh.h"
#include "elements.h"
#include "assembly.h"
#include "postprocess.h"
#include "linearsolver.h"

class simulation {

    protected:
        DCBMeshParameters mp;

        std::vector<std::vector<double>> nodes;
        std::vector<std::vector<int>> elements;
        NodeSets sets;
        
        Eigen::VectorXd u,s;
        std::vector<double> H, exx, eyy, exy, tr_e;
        std::vector<Eigen::Triplet<double>> triplets;

        Eigen::SparseMatrix<double> K;
        Eigen::VectorXd F, d, Dd, residual;

        std::vector<double> F1_history, F2_history, D1_history, D2_history, a_history, U_history;
        std::vector<double> Gc_history;

        // solver tools
        LinearSolver solver;
        bool firstLoop = true;

        // geometry and material parameter
        double ls, t, a_min, E_int, E_bulk, nu_int, nu_bulk, Gc_int, Gc_bulk, Gc_eff;

        // run parameters
        double tol, relax, s_threshold, D, Dinc, Dend;
        int quadratureDegree, incr_count, save_freq, restart_from, prec;
        bool isPlaneStress;

        // filenames
        std::string Hfilename = "H,inc=", usfilename = "inc=", history_filename = "history.csv";
        std::ofstream history_file;

        // boundary condition mapping
        std::vector<int> constrained;
        std::vector<double> constrainedValues;
        std::vector<bool> isConstrained;
        std::vector<int> freeMap;
        std::unordered_map<int, double*> constrainedMap;

    public:
        virtual ~simulation() = default;

        void Run() {
            Setup();
            Initialise(); 

            if (restart_from > 0) {
                HandleRestart(); 
            } else {
                
                for(int node : sets.initial_crack_nodes) s[node] = 0.0;
            }

            while (D <= Dend + Dinc/10.0) {
                auto t1 = std::chrono::high_resolution_clock::now();
                printf("Increment number: %d, Applied displacement: %.*f\n", incr_count, prec, D);

                ApplyBoundaryConditions();

                if (!SolveStep()) { // solve and check if not solved
                    std::cerr << "Newton-Raphson failed to converge!" << std::endl;
                    break;
                }

                SaveResults();

                auto t2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> seconds = t2 - t1;
                printf("Increment %d took %.3f seconds\n", incr_count, seconds.count());
                std::cout << std::endl;

                D += Dinc;
                incr_count++;

            }

        }
    
    protected:
        virtual void Setup() = 0;
        virtual void ApplyBoundaryConditions() = 0;
        virtual void SaveResults() = 0;
        virtual double CalculateGcForRestart() = 0;
        // = 0 leaves open to intepretation depending on run parameters 

        std::vector<std::string> split(const std::string& line, char delimiter) {
            std::vector<std::string> tokens;
            std::string token;
            std::istringstream tokenStream(line);
            while (std::getline(tokenStream, token, delimiter)) tokens.push_back(token);
            return tokens;
        }

        void IdentifyNodeSets() {
            sets = NodeSets();

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

        void SaveGeometryData() {
            std::ofstream outfile("geometry.csv");
            if(outfile.is_open()) {
                outfile << "L,h1,h2,hi,a0\n";
                outfile << std::fixed << std::setprecision(prec) 
                        << mp.L << "," << mp.h1 << "," 
                        << mp.h2 << "," << mp.hi << "," << mp.a0 << std::endl;
                outfile.close();
            }
        }

        void Initialise() {
            int numNodes = nodes.size();
            int totalDOFs = 3 * numNodes;
            

            u = Eigen::VectorXd::Zero(2*numNodes);
            s = Eigen::VectorXd::Ones(numNodes);

            K.resize(totalDOFs, totalDOFs);
            F.resize(totalDOFs);
            d.resize(totalDOFs);
            Dd.resize(totalDOFs);
            
            Eigen::VectorXd xi, eta, w;
            Get2DQuadrature(quadratureDegree, xi, eta, w);
            size_t totalIPs = elements.size() * w.size();

            H.assign(totalIPs, 0.0);
            exx.assign(totalIPs, 0.0);
            eyy.assign(totalIPs, 0.0);
            exy.assign(totalIPs, 0.0);
            tr_e.assign(totalIPs, 0.0);

            constrained.clear();
            for(int node : sets.initial_crack_nodes) constrained.push_back(2*numNodes + node);
            for(int node : sets.support_nodes) { constrained.push_back(2*node); constrained.push_back(2*node + 1); };
            for(int node : sets.lower_arm_end_nodes) constrained.push_back(2*node + 1);
            for(int node : sets.upper_arm_end_nodes) constrained.push_back(2*node + 1);
    
            constrainedValues.assign(constrained.size(), 0.0);
            isConstrained.assign(3*numNodes, false);
            freeMap.assign(3*numNodes, -1);

            for(size_t k = 0; k < constrained.size(); ++k) {
                constrainedMap[constrained[k]] = &constrainedValues[k];
                isConstrained[constrained[k]] = true;
            }

            int nFree = 0;
            for(int i = 0; i < 3*numNodes; ++i) if(!isConstrained[i]) freeMap[i] = nFree++;

            if (restart_from == 0) {
                incr_count = 0;
                history_file.open(history_filename);
                history_file << "inc,D1,F1,D2,F2,a,G_c,U" << std::endl;
                history_file << std::fixed << std::setprecision(prec);  
            }   
        }
        void HandleRestart() {
            std::cout << "Handling restart, loading past data" << std::endl;

            // read history.csv
            std::ifstream file_in(history_filename);
            std::vector<std::string> valid_lines;
            std::string line;
            if (file_in.is_open()) {
                if (std::getline(file_in, line)) valid_lines.push_back(line);
                    while (std::getline(file_in, line)) {
                    auto fields = split(line, ',');
                    if (fields.size() < 8) continue;

                    int inc = std::stoi(fields[0]);
                    if (inc > restart_from) break; // ignore 'the future'

                    valid_lines.push_back(line);

                    // refill the history vectors one line at a time
                    D1_history.push_back(std::stod(fields[1]));
                    F1_history.push_back(std::stod(fields[2]));
                    D2_history.push_back(std::stod(fields[3]));
                    F2_history.push_back(std::stod(fields[4]));
                    a_history.push_back(std::stod(fields[5]));
                    U_history.push_back(std::stod(fields[7]));

                    Gc_history.push_back(CalculateGcForRestart()); // recalculate Gc depending on a_min
                }
                file_in.close();
            }

            // truncate saved file (wipe and write to ensure data safety)
            history_file.open(history_filename, std::ios::out);
            for (const auto& l : valid_lines) history_file << l << "\n";
            history_file.close();

            // reopen for current run
            history_file.open(history_filename, std::ios::app);
            history_file << std::fixed << std::setprecision(prec);

            // LOAD NODAL DATA (u, v, s) from "inc=X.csv"
            // Structure: x[0], y[1], u[2], v[3], s[4]
            std::string node_file = usfilename + std::to_string(restart_from) + ".csv";
            std::ifstream f_node(node_file);
            if (f_node.is_open()) {
                std::getline(f_node, line); // Header
                int i = 0;
                while (std::getline(f_node, line)) {
                    auto fields = split(line, ',');
                    u(2*i)     = std::stod(fields[2]); // u (x-disp)
                    u(2*i + 1) = std::stod(fields[3]); // v (y-disp)
                    s(i)       = std::stod(fields[4]); // damage
                    i++;
                }
                f_node.close();
            }

            // LOAD INT-POINT DATA (H) from "H,inc=X.csv"
            // Structure: x[0], y[1], H[2], exx[3], eyy[4], exy[5], tre[6]
            std::string h_file = Hfilename + std::to_string(restart_from) + ".csv";
            std::ifstream f_h(h_file);
            if (f_h.is_open()) {
                std::getline(f_h, line); // Header
                int idx = 0;
                while (std::getline(f_h, line) && idx < H.size()) {
                    auto fields = split(line, ',');
                    H[idx] = std::stod(fields[2]); // Grab H, ignore strains
                    idx++;
                }
                f_h.close();
            }

            D = D1_history.back() + Dinc;
            incr_count = restart_from + 1;

        }

        bool SolveStep(){
                int numNodes = nodes.size();
                int totalDOFs = 3*numNodes;

                int nFree = 0;
                for(int m : freeMap) if(m != -1) nFree = std::max(nFree, m + 1);
                Eigen::SparseMatrix<double> KR(nFree, nFree);
                Eigen::VectorXd d_reduced(nFree), FR(nFree), Dd_reduced(nFree);

                residual.resize(nFree);

                d.head(u.size()) = u; d.tail(s.size()) = s;

                int iter_count = 1;
                while(true) {
                    AssembleAll(nodes, elements, u, s, H, E_int, E_bulk, nu_int, nu_bulk, 
                                Gc_eff, Gc_bulk, t, mp.h2, mp.hi, ls, K, F, triplets, 
                                isPlaneStress, mp.isQuadratic, quadratureDegree);

                    ApplyDirichletBC(K, F, isConstrained, constrainedMap, freeMap, KR, FR);

                    for(int i = 0; i < totalDOFs; ++i) if(freeMap[i] != -1) d_reduced[freeMap[i]] = d[i];
                    residual = FR - KR * d_reduced;
                    double R_norm = residual.norm() / FR.norm();
                    printf("Iteration: %d, |R|/|F| = %.3e\n", iter_count, R_norm);

                    if(R_norm < tol) break;
                    if(iter_count > 50) return false;

                    if(firstLoop) { solver.analyzePattern(KR); firstLoop = false; }
                    if(!solver.solve(KR, residual, Dd_reduced)) return false;

                    for(int i = 0; i < totalDOFs; ++i) Dd[i] = (freeMap[i] != -1) ? Dd_reduced[freeMap[i]] : 0.0;
                    d += relax * Dd;
                    for(auto const& pair : constrainedMap) d[pair.first] = *pair.second;

                    u = d.head(u.size()); s = d.tail(s.size());
                    iter_count++;
                }
                UpdateH(nodes, elements, u, H, E_int, E_bulk, nu_int, nu_bulk, mp.h2, mp.hi, 
                    ls, isPlaneStress, mp.isQuadratic, quadratureDegree, exx, eyy, exy, tr_e);
            return true;
        }
        

};

#endif