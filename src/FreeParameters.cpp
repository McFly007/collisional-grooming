#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Grid3D.h"
#include "Outputs.h"
#include "Settings.h"
#include "Records.h"
#include "FreeParameters.h"



std::vector<float> unique_diameters;	 // Contains unique diameters from all datafiles - this is done in the init() function
std::vector<double> SFD_profile;		 // Contains unique diameters from all datafiles - this is done in the init() function
std::vector<double> SFD_profile_Kessler; // Contains unique diameters from all datafiles - this is done in the init() function
std::vector<vector<double>> SFD_array;
std::vector<double> temp_diameters; // Mid-values of SFD diameters 
std::vector<vector<double>> Integration_Setup_array;


// Default values - these are changed using the input file
double alpha = 4.1;
double beta = 3.7;
double Dmid = 60.0;
double As = 1e6 * 1e-7 * 1e3; // J kg-1 - Using Krivov et al. (2006)
double Bs = -0.24;
double pericenter_index = -1.3;
double Dmax = 3000.0;
double N0 = 1.0;
double particle_density = 2000; // kg m^-3
double As_modifier = 1;

// Global variables
int total_number_of_particles_integration;
double initial_weight = 1.0;
double intended_cloud_area = 2e17; // This is the total cloud cross-section we want to achieve