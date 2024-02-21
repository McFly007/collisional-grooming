
#ifndef Records_h_
#define Records_h_

#include "Grid3D.h"
#include "Settings.h"
#include "FreeParameters.h"

void record_disk_profiles(petrpokorny::Grid3D *grid3Ddata, Vec3 point, float weight, float diameter, float collision_velocity); 
void record_disk_mass_loss(petrpokorny::Grid3D *grid3Ddata, Vec3 point, float mass_loss);
void record_radialProfile(std::vector<double> &r_prof, Vec3 point, float weight);
void record_radialProfileCross(std::vector<double> &r_prof, Vec3 point, float weight, float cross_section);
void record_eclipticLatitudeProfile(std::vector<double> &ecl_prof, std::vector<double> &ecl_prof_mass, Vec3 point, float weight, float diameter, float particle_density);
void record_SFDProfile(std::vector<double> &rec_profile, float weight, float diameter);
#endif
