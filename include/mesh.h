#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <map>

/**
 * Comprehensive mesh generation and management for DCB and ELS * 
*/

// mesh operations
void insert_midpoints(std::vector<double>& vec);
void meshgrid(const std::vector<double>& x, const std::vector<double>& y, std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y);
void create_nodes_matrix(const std::vector<double>& x, const std::vector<double>& y, std::vector<std::vector<double>>& nodes);
void create_elements_matrix(const std::vector<std::vector<double>>& X, const std::vector<std::vector<double>>& Y, bool isQuadratic, std::vector<std::vector<int>>& elements);
void add_crack(std::vector<std::vector<double>>& nodes, std::vector<std::vector<int>>& elements, double xct, double yc, double tol);

// container for geometric and mesh settings
struct DCBMeshParameters {
    // Geometry dimensions
    double L, h1, h2, hi, a0;
    
    // Mesh density settings
    double Delta_y, Delta_y_max_ahead, Delta_y_max_behind;
    double GR_y_ahead, GR_y_behind;
    double Delta, Delta_x, Delta_x_max;
    double yminus, yplus, xplus, xminus;
    double GR_x;
    
    // Element type
    bool isQuadratic;
};

// containers for node sets that will have boundary conditions applied
struct NodeSets {
	std::vector<int> initial_crack_nodes;
	std::vector<int> centreline_nodes; 
	std::vector<int> interface_nodes; 
	std::vector<int> upper_arm_end_nodes;
	std::vector<int> lower_arm_end_nodes;
	std::vector<int> support_nodes;
};

// builds the specific DCB/ELS geometry
void GenerateDCBMesh
(
    const DCBMeshParameters& p, std::vector<std::vector<double>>& nodes, std::vector<std::vector<int>>& elements
);

NodeSets IdentifyNodeSets
(
    const std::vector<std::vector<double>>& nodes, const DCBMeshParameters& p
);

#endif