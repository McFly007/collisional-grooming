

#ifndef Constants_h_
#define Constants_h_

// Definitions of constants
const double auday2ms = double(149597900000e0 / 86400.0); // au per day to meters per second conversion
const double au = double(149597900000e0);				  // au in meters
const double GM_sun = double(1.32712440018e20);			  // Graviational parameter of the Sun
const double GM_earth = double(3.986004418e14);			  // Graviational parameter of the Earth
const double deg2rad = double(0.0174532925199e0);		  // Conversion from degrees to radians
const double rad2deg = double(1.0 / deg2rad);

#endif