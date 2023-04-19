
#ifndef Outputs_h_
#define Outputs_h_

struct Outputs {

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
	std::vector<vector<double>> semimajor_axis_profile_number_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> semimajor_axis_profile_brightness_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> semimajor_axis_profile_cross_section_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> semimajor_axis_profile_mass_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> collisional_velocity_radial_profile_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> collisional_lifetime_radial_profile_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> collisional_velocity_semimajor_profile_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index
	std::vector<vector<double>> collisional_lifetime_semimajor_profile_SFD; // dimension 1 is the semimajor axis, dimension 2 is the particle diameter index

};

#endif


