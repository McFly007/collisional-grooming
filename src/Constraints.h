
#ifndef Constraints_h_
#define Constraints_h_

struct Constraints {
	std::vector<float> CMOR_AH_a_x;
	std::vector<float> CMOR_AH_a_y;
	std::vector<float> CMOR_AH_vg_x;
	std::vector<float> CMOR_AH_vg_y;
	std::vector<float> CMOR_AH_ecc_x;
	std::vector<float> CMOR_AH_ecc_y;
	std::vector<float> CMOR_AH_inc_x;
	std::vector<float> CMOR_AH_inc_y;

	std::vector<float> AMOR_AH_a_x;
	std::vector<float> AMOR_AH_a_y;
	std::vector<float> AMOR_AH_vg_x;
	std::vector<float> AMOR_AH_vg_y;
	std::vector<float> AMOR_AH_ecc_x;
	std::vector<float> AMOR_AH_ecc_y;
	std::vector<float> AMOR_AH_inc_x;
	std::vector<float> AMOR_AH_inc_y;
};

struct Constraints_Double {
	std::vector<double> CMOR_AH_a_x;
	std::vector<double> CMOR_AH_a_y;
	std::vector<double> CMOR_AH_vg_x;
	std::vector<double> CMOR_AH_vg_y;
	std::vector<double> CMOR_AH_ecc_x;
	std::vector<double> CMOR_AH_ecc_y;
	std::vector<double> CMOR_AH_inc_x;
	std::vector<double> CMOR_AH_inc_y;

	std::vector<double> AMOR_AH_a_x;
	std::vector<double> AMOR_AH_a_y;
	std::vector<double> AMOR_AH_vg_x;
	std::vector<double> AMOR_AH_vg_y;
	std::vector<double> AMOR_AH_ecc_x;
	std::vector<double> AMOR_AH_ecc_y;
	std::vector<double> AMOR_AH_inc_x;
	std::vector<double> AMOR_AH_inc_y;

	std::vector<double> SFD_profile;
	double mass_accreted_at_Earth_Kessler;
	double zodiacal_cloud_total_cross_section;
	double zodiacal_cloud_total_cross_section_beyond_1au;
};

#endif


