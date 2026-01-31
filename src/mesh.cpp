#include "mesh.h"


// =================
// Helper functions
// =================

std::vector<double> linspace(double start, double end, double initial_size, double growth_rate = 1.0, double tol = 1e-9)
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

std::vector<double> linspace2(double start, double end, double initial_size, double growth_rate, double max_size, double tol = 1e-9)
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

// Crack insertion logic
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

// ========================
// Main DCB mesh generator
// ========================
// Specimen layout:
// Y-Direction: y1(Bulk) -> y2(Refined) -> y3/y4(Interface) -> y5(Refined) -> y6(Bulk)
// X-Direction: x1(Bulk) -> x2/x3(Refinement around Tip) -> x4(Bulk/Cracked)
void GenerateDCBMesh(const DCBMeshParameters& p, std::vector<std::vector<double>>& nodes, std::vector<std::vector<int>>& elements)
{
    // Clear outputs just in case
    nodes.clear();
    elements.clear();

    // 1. Generate X Coordinates

    // x1: From support (0) to start of refined region
    std::vector<double> x1 = linspace2(0, p.L - p.a0 - p.yminus, p.Delta_y, p.GR_y_ahead, p.Delta_y_max_ahead);
    std::reverse(x1.begin(), x1.end());
    for(auto& elem : x1) { elem *= -1; elem += (p.L - p.a0 - p.yminus); }
    
    // x2: Refined region ahead of crack tip
    std::vector<double> x2 = linspace(p.L - p.a0 - p.yminus, p.L - p.a0, p.Delta_y); 
    
    // x3: Refined region behind crack tip
    std::vector<double> x3 = linspace(p.L - p.a0, p.L - p.a0 + p.yplus, p.Delta_y); 
    
    // x4: Cracked region (rest of the beam)
    std::vector<double> x4 = linspace2(p.L - p.a0 + p.yplus, p.L, p.Delta_y, p.GR_y_behind, p.Delta_y_max_behind); 

    // Merge X
    std::vector<double> x = x1;
    x.insert(x.end(), x2.begin() + 1, x2.end());
    x.insert(x.end(), x3.begin() + 1, x3.end());
    x.insert(x.end(), x4.begin() + 1, x4.end());


    // 2. Generate Y Coordinates

    // y6: Upper arm
    std::vector<double> y6 = linspace2(p.h2 + p.hi + p.xplus, p.h2 + p.hi + p.h1, p.Delta_x, p.GR_x, p.Delta_x_max); 
    
    // y5: Refined region above interface
    std::vector<double> y5 = linspace(p.h2 + p.hi, p.h2 + p.hi + p.xplus, p.Delta); 
    
    // y4: Upper interface
    std::vector<double> y4 = linspace(p.h2 + p.hi/2, p.h2 + p.hi, p.Delta); 
    
    // y3: Lower interface
    std::vector<double> y3 = linspace(p.h2, p.h2 + p.hi/2, p.Delta); 
    
    // y2: Refined region below interface
    std::vector<double> y2 = linspace(p.h2 - p.xplus, p.h2, p.Delta); 
    
    // y1: Lower arm (mirrored about y = h2/2 roughly logic)
    std::vector<double> y1 = linspace2(0, p.h2 - p.xplus, p.Delta_x, p.GR_x, p.Delta_x_max); 
    std::reverse(y1.begin(), y1.end());
    for(auto& elem : y1) { elem *= -1; elem += (p.h2 - p.xplus); }
    
    // Merge Y
    std::vector<double> y = y1;
    y.insert(y.end(), y2.begin() + 1, y2.end());
    y.insert(y.end(), y3.begin() + 1, y3.end());
    y.insert(y.end(), y4.begin() + 1, y4.end());
    y.insert(y.end(), y5.begin() + 1, y5.end());
    y.insert(y.end(), y6.begin() + 1, y6.end());


    // 3. Handle Quadratic Midpoints
    if(p.isQuadratic) {
        insert_midpoints(x);
        insert_midpoints(y);
    }

    // 4. Create Grid (Meshgrid)
    std::vector<std::vector<double>> X, Y;
    meshgrid(x, y, X, Y);

    // 5. Create Nodes and Elements Matrices
    create_nodes_matrix(x, y, nodes);
    create_elements_matrix(X, Y, p.isQuadratic, elements);

    // 6. Split nodes for the initial crack
    double findtol = p.Delta/10;
    add_crack(nodes, elements, p.L - p.a0, p.h2 + p.hi/2, findtol);
}

// Find nodes to group into prexisting structure for boudary conditions
NodeSets IdentifyNodeSets(const std::vector<std::vector<double>>& nodes, const DCBMeshParameters& p) {
    NodeSets sets;
    double findtol = p.Delta / 10.0;
	for(size_t i = 0; i < nodes.size(); ++i)
	{
		double x = nodes[i][0];
		double y = nodes[i][1];
		
		// set of nodes forming the initial crack
		if((x >= (p.L - p.a0) - findtol) && (std::abs(y - (p.h2 + p.hi/2)) < findtol))
		{
			sets.initial_crack_nodes.push_back(i);
		}
		
		// set of nodes forming the centreline of the interface
		if((x <= (p.L - p.a0) + findtol) && (std::abs(y - (p.h2 + p.hi/2)) < findtol))
		{
			sets.centreline_nodes.push_back(i);
		}
		
		// set of all nodes forming the interface
		if((x <= (p.L - p.a0) + findtol) && (y >= p.h2 - findtol) && (y <= p.h2 + p.hi + findtol))
		{
			sets.interface_nodes.push_back(i);
		}
		
		// set of nodes on end of upper arm
		if((y >= (p.h2 + p.hi) - findtol) && (std::abs(x - p.L) < findtol))
		{
			sets.upper_arm_end_nodes.push_back(i);
		}
		
		// set of nodes on end of lower arm
		if((y <= p.h2 + findtol) && (std::abs(x - p.L) < findtol))
		{
			sets.lower_arm_end_nodes.push_back(i);
		}
		
		// set of nodes on support end
		if((std::abs(x) < findtol))
		{
			sets.support_nodes.push_back(i);
		}
	}

    return sets;
}


