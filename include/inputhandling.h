#ifndef INPUTHANDLING_H
#define INPUTHANDLING_H

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>


struct InputData {
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

inline double find_val(const std::string& filename, const std::string& key) {
    std::ifstream file(filename);
    std::string line;
    std::string search_key = key + "=";
    
    if (file.is_open()) {
        while (std::getline(file, line)) {
            // Remove whitespace/newlines for safety
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            
            // Check if line starts with "L=" (for example)
            if (line.find(search_key) == 0) { 
                size_t start = line.find('=') + 1;
                size_t end = line.find_first_of(";#"); // Stop at comments
                if (end == std::string::npos) end = line.length();
                
                try {
                    return std::stod(line.substr(start, end - start));
                } catch (...) {
                    return 0.0;
                }
            }
        }
    }
    std::cerr << "Warning: Could not find parameter '" << key << "' in snapshot." << std::endl;
    return 0.0; 
}


#define SAVE_PARAM(STREAM, DATA, VAR) STREAM << #VAR << "=" << DATA.VAR << "\n"
#define LOAD_PARAM(DATA, VAR, PATH) DATA.VAR = (decltype(DATA.VAR))find_val(PATH, #VAR)

// Function to SAVE parameters
inline void save_parameters(const std::string& filename, const InputData& data) {
    std::ofstream out(filename);
    out.precision(16);

    SAVE_PARAM(out, data, L);       SAVE_PARAM(out, data, h1);      SAVE_PARAM(out, data, h2);
    SAVE_PARAM(out, data, hi);      SAVE_PARAM(out, data, a0);      SAVE_PARAM(out, data, ls);
    SAVE_PARAM(out, data, t);       SAVE_PARAM(out, data, a_min);
    
    SAVE_PARAM(out, data, E_int);   SAVE_PARAM(out, data, nu_int);  SAVE_PARAM(out, data, Gc_int);
    SAVE_PARAM(out, data, Gc_eff);  SAVE_PARAM(out, data, E_bulk);  SAVE_PARAM(out, data, nu_bulk);
    SAVE_PARAM(out, data, Gc_bulk); SAVE_PARAM(out, data, isPlaneStress);
    
    SAVE_PARAM(out, data, Dinit);   SAVE_PARAM(out, data, Dend);    SAVE_PARAM(out, data, Dinc);
    SAVE_PARAM(out, data, tol);     SAVE_PARAM(out, data, relax);   SAVE_PARAM(out, data, save_freq);
    
    SAVE_PARAM(out, data, Delta_y); SAVE_PARAM(out, data, Delta_y_max_ahead); 
    SAVE_PARAM(out, data, Delta_y_max_behind); SAVE_PARAM(out, data, GR_y_ahead);
    SAVE_PARAM(out, data, GR_y_behind); SAVE_PARAM(out, data, Delta);
    SAVE_PARAM(out, data, Delta_x); SAVE_PARAM(out, data, Delta_x_max);
    SAVE_PARAM(out, data, yminus);  SAVE_PARAM(out, data, yplus);
    SAVE_PARAM(out, data, xplus);   SAVE_PARAM(out, data, xminus);
    SAVE_PARAM(out, data, GR_x);    SAVE_PARAM(out, data, isQuadratic); 
    SAVE_PARAM(out, data, quadratureDegree);

    out.close();
    std::cout << ">> Parameters saved to " << filename << std::endl;
}

// Function to LOAD parameters
inline void load_parameters(const std::string& filename, InputData& data) {
    std::cout << ">> Loading snapshot from: " << filename << std::endl;

    LOAD_PARAM(data, L, filename);       LOAD_PARAM(data, h1, filename);      LOAD_PARAM(data, h2, filename);
    LOAD_PARAM(data, hi, filename);      LOAD_PARAM(data, a0, filename);      LOAD_PARAM(data, ls, filename);
    LOAD_PARAM(data, t, filename);       LOAD_PARAM(data, a_min, filename);
    
    LOAD_PARAM(data, E_int, filename);   LOAD_PARAM(data, nu_int, filename);  LOAD_PARAM(data, Gc_int, filename);
    LOAD_PARAM(data, Gc_eff, filename);  LOAD_PARAM(data, E_bulk, filename);  LOAD_PARAM(data, nu_bulk, filename);
    LOAD_PARAM(data, Gc_bulk, filename); LOAD_PARAM(data, isPlaneStress, filename);
    
    LOAD_PARAM(data, Dinit, filename);   LOAD_PARAM(data, Dend, filename);    LOAD_PARAM(data, Dinc, filename);
    LOAD_PARAM(data, tol, filename);     LOAD_PARAM(data, relax, filename);   LOAD_PARAM(data, save_freq, filename);
    
    LOAD_PARAM(data, Delta_y, filename); LOAD_PARAM(data, Delta_y_max_ahead, filename); 
    LOAD_PARAM(data, Delta_y_max_behind, filename); LOAD_PARAM(data, GR_y_ahead, filename);
    LOAD_PARAM(data, GR_y_behind, filename); LOAD_PARAM(data, Delta, filename);
    LOAD_PARAM(data, Delta_x, filename); LOAD_PARAM(data, Delta_x_max, filename);
    LOAD_PARAM(data, yminus, filename);  LOAD_PARAM(data, yplus, filename);
    LOAD_PARAM(data, xplus, filename);   LOAD_PARAM(data, xminus, filename);
    LOAD_PARAM(data, GR_x, filename);    LOAD_PARAM(data, isQuadratic, filename); 
    LOAD_PARAM(data, quadratureDegree, filename);
}

#endif