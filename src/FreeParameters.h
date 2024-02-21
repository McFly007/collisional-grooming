#ifndef FreeParameters_h_
#define FreeParameters_h_

#include "Grid3D.h"
#include "Settings.h"

extern double As; 
extern double Bs;
extern double alpha; // Size-frequency differential index alpha
extern double beta; // Size-frequency differential index beta
extern double Dmid; // In micrometers
extern double N0; // Normalization for the size-frequency distribution
extern double Dmax; // Maximum diameter for the size-frequency distribution
extern double particle_density; // In kg m^-3
extern double pericenter_index;
extern double As_modifier;

extern int total_number_of_particles_integration;
extern double initial_weight;
extern double intended_cloud_area; // This is the total cloud cross-section we want to achieve


extern std::vector<float> unique_diameters;	 // Contains unique diameters from all datafiles - this is done in the init() function
extern std::vector<double> SFD_profile;		 // Contains unique diameters from all datafiles - this is done in the init() function
extern std::vector<double> SFD_profile_Kessler; // Contains unique diameters from all datafiles - this is done in the init() function
extern std::vector<vector<double>> SFD_array;
extern std::vector<double> temp_diameters; // Mid-values of SFD diameters 
extern std::vector<vector<double>> Integration_Setup_array;

#endif