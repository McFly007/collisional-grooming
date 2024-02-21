#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Grid3D.h"
#include "Outputs.h"
#include "Settings.h"
#include "Records.h"
#include "Elementary_Functions.h"
#include "FreeParameters.h"


// Define namespaces
using namespace petrpokorny;
using namespace std;

void record_disk_profiles(Grid3D *grid3Ddata, Vec3 point, float weight, float diameter, float collision_velocity)
{
	std::vector<int> grid_index = grid3Ddata->getIndexes(point); // x,y, and z

	int x = grid_index[0];
	int y = grid_index[1];
	int z = grid_index[2];

	float heliocentric_distance = sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
	float flux_at_distance = 1.0 / (heliocentric_distance * heliocentric_distance);
	float cross_section = pow(diameter, 2.0) * 0.25 * M_PI * 1e-12; // diameter is in micrometers

	My_Outputs.disk_profile_face_on_number[x][y] += weight;
	My_Outputs.disk_profile_edge_on_number[x][z] += weight;

	My_Outputs.disk_profile_face_on_brightness[x][y] += weight * flux_at_distance * cross_section;
	My_Outputs.disk_profile_edge_on_brightness[x][z] += weight * flux_at_distance * cross_section;

	My_Outputs.disk_profile_face_on_impact_velocity[x][y] += weight * collision_velocity;
	My_Outputs.disk_profile_edge_on_impact_velocity[x][z] += weight * collision_velocity;
}

void record_disk_mass_loss(Grid3D *grid3Ddata, Vec3 point, float mass_loss)
{
	std::vector<int> grid_index = grid3Ddata->getIndexes(point); // x,y, and z

	int x = grid_index[0];
	int y = grid_index[1];
	int z = grid_index[2];

	float ecliptic_distance = sqrt(point.x * point.x + point.y * point.y);

	My_Outputs.disk_profile_face_on_massloss[x][y] += mass_loss;
	My_Outputs.disk_profile_edge_on_massloss[x][z] += mass_loss;

	// For central bins only (near the ecliptic)
	////////////
	// int central_index = floor(octree_limit/dr.z * 0.5); // getting the ecliptic index
	// if( z == central_index) {My_Outputs.disk_profile_face_on_massloss[x][y] += mass_loss;} // in the ecliptic
	// if( y == central_index) {My_Outputs.disk_profile_edge_on_massloss[x][z] += mass_loss;} // in the y-plane		


	int gridPoint_int = floor(ecliptic_distance / (settings_grid_half_size) * (int)My_Outputs.disk_profile_radial_massloss.size());
	if (gridPoint_int < (int)My_Outputs.disk_profile_radial_massloss.size())
	{
		My_Outputs.disk_profile_radial_massloss[gridPoint_int] += mass_loss;
	}
	// if(gridPoint_int < (int)My_Outputs.disk_profile_radial_massloss.size() && z == central_index ) {My_Outputs.disk_profile_radial_massloss[gridPoint_int] += mass_loss;}
	// We count only the points at the ecliptic; i.e. the central index
}

void record_radialProfile(std::vector<double> &r_prof, Vec3 point, float weight)
{
	//Records the ecliptic distance profile of particle weights
	float gridPoint_rad;
	int gridPoint_int;
	gridPoint_rad = sqrt(point.x * point.x + point.y * point.y);
	gridPoint_int = floor(gridPoint_rad / (settings_dr.x));
	if (gridPoint_int < (int)r_prof.size())
	{
		r_prof[gridPoint_int] += weight;
	}
}

void record_radialProfileCross(std::vector<double> &r_prof, Vec3 point, float weight, float cross_section)
{
	//Records the ecliptic distance profile of particle weights * cross_section
	float gridPoint_rad;
	int gridPoint_int;
	gridPoint_rad = sqrt(point.x * point.x + point.y * point.y);
	gridPoint_int = floor(gridPoint_rad / (settings_dr.x));
	if (gridPoint_int < (int)r_prof.size())
	{
		r_prof[gridPoint_int] += weight * cross_section;
	}
}

void record_eclipticLatitudeProfile(std::vector<double> &ecl_prof, std::vector<double> &ecl_prof_mass, Vec3 point, float weight, float diameter, float particle_density)
{
	// Ecliptic latitude profile beyond 1 au - we calculate the cross-section per degree of the ecliptic latitude + We also calculate the mass for the same profile
	// These two quantities are easy to compare to Nesvorny et al. (2011) etc.
	float ecliptic_latitude;
	int gridPoint_int;
	double r_eclipt = sqrt(point.x * point.x + point.y * point.y);
	double d_phi = 180.0 / float(ecl_prof.size());

	if (r_eclipt > 1.0)
	{
		ecliptic_latitude = atan(point.z / (r_eclipt - 1.0)) * 180.0 / M_PI; // We are assuming that the particles in the cloud are not in 1:1 resonance with Earth and thus can be randomly oriented with respect to Earth
		gridPoint_int = floor((ecliptic_latitude + 90.0) * d_phi);
		ecl_prof[gridPoint_int] += weight * pow(diameter, 2.0) * 0.25 * M_PI * 1e-12;						 // Weight * Particle Cross-section - diameters are in micrometers so we have *1e-12
		ecl_prof_mass[gridPoint_int] += weight * pow(diameter, 3.0) / 6.0 * M_PI * 1e-18 * particle_density; // Weight * Particle Cross-section - diameters are in micrometers so we have *1e-12
	}
}


void record_SFDProfile(std::vector<double> &rec_profile, float weight, float diameter)
{
	int dia = find_diameter_index(diameter);
	rec_profile[dia] += weight;
}