
#ifndef Outputs_h_
#define Outputs_h_

#include "Grid3D.h"
#include "Settings.h"

extern	std::vector<double> r_prof_init;		 // Vector that keeps the radial profile for the 1st iteration
extern	std::vector<double> r_prof_iter;		 // Vector that keeps the radial profile for the 2nd iteration
extern	std::vector<double> r_prof_brightness; // Vector that keeps the radial profile for the 2nd iteration
extern	std::vector<double> ecl_profile;		 // Vector that keeps the ecliptic latitude profile
extern	std::vector<double> ecl_profile_mass;	 // Vector that keeps the ecliptic latitude profile

struct Outputs
{
	double mass_accreted_at_Earth;
	std::vector<vector<double>> disk_profile_face_on_number; // x,y, flux
	std::vector<vector<double>> disk_profile_face_on_brightness;
	std::vector<vector<double>> disk_profile_face_on_impact_velocity;

	std::vector<vector<double>> disk_profile_edge_on_number;
	std::vector<vector<double>> disk_profile_edge_on_brightness;
	std::vector<vector<double>> disk_profile_edge_on_impact_velocity;

	std::vector<vector<double>> disk_profile_face_on_massloss;
	std::vector<vector<double>> disk_profile_edge_on_massloss;
	std::vector<vector<double>> disk_profile_radial;

	std::vector<double> disk_profile_radial_massloss;
	std::vector<vector<double>> semimajor_axis_profile_number_SFD;			// dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> semimajor_axis_profile_brightness_SFD;		// dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> semimajor_axis_profile_cross_section_SFD;	// dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> semimajor_axis_profile_mass_SFD;			// dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> collisional_velocity_radial_profile_SFD;	// dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> collisional_lifetime_radial_profile_SFD;	// dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> collisional_velocity_semimajor_profile_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> collisional_lifetime_semimajor_profile_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
};

extern Outputs My_Outputs;

void Initiate_Output_Files();
void Initialize_1D_Array(std::vector<double> &init_array, int xdim);
void Initialize_2D_Array(std::vector<std::vector<double>> &init_array, int xdim, int ydim);
void Initialize_Outputs(petrpokorny::Grid3D *grid3Ddata);
void Initiate_Output_Files();
void Print_Outputs_Big(petrpokorny::Grid3D *grid3Ddata, int file_counter, double refactor);
void Print_Outputs_Small(petrpokorny::Grid3D *grid3Ddata, int file_counter, double refactor);
void Zero_Outputs(petrpokorny::Grid3D *grid3Ddata);


#endif
