

#ifndef Settings_h_
#define Settings_h_

#include "math.h"
#include "Constants.h"

const double angle_diff = cos(10.0 / 180.0 * M_PI); // 10 degrees - the angular difference between two particle velocity vectors - for merging
const double ratio_diff = 0.10;                     // 10% - the velocity magnitude difference between two vectors - for merging
const Vec3 settings_dr = Vec3(.05, .05, .05);             // Size of the 3D bin
const double volume_unit = settings_dr.x * settings_dr.y * settings_dr.z * pow(au, 3.0); // Volume of one binning box in m3
const double settings_grid_half_size = 10.0; // maximum distance from the grid center

#endif