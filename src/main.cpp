#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cfenv>
#define _USE_MATH_DEFINES
#include "math.h"
#include "Octree.h"
#include "Grid3D.h"
#include "Stopwatch.h"
#include "Constraints.h"
#include "Outputs.h"
#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split

// Definitions of constants
const double auday2ms = double(149597900000e0/86400.0); // au per day to meters per second conversion
const double au = double(149597900000e0); // au in meters
const double GM_sun = double(1.32712440018e20); // Graviational parameter of the Sun
const double GM_earth = double(3.986004418e14); // Graviational parameter of the Earth
const double deg2rad = double(0.0174532925199e0); // Conversion from degrees to radians
const double rad2deg = double(1.0/deg2rad);
using namespace brandonpelfrey;
using namespace petrpokorny;
using namespace std;

// Global variables
double angle_diff = cos(10.0/180.0*M_PI); // 10 degrees - the angular difference between two particle velocity vectors - for merging
double ratio_diff = 0.10; // 10% - the velocity magnitude difference between two vectors - for merging
Vec3 dr = Vec3(.025,.025,.025);
double volume_unit = dr.x*dr.y*dr.z*8.0*pow(au,3.0); // Volume of one binning box in m3
double octree_limit = 10.0;
double initial_weight = 1.0;
double intended_cloud_area = 2e17; // This is the total cloud cross-section we want to achieve
double mass_accreted_at_Earth =0.0;
int total_number_of_particles_integration;
std::vector<double> TEST_probabilities;

std::vector<float> unique_diameters; // Contains unique diameters from all datafiles - this is done in the init() function
std::vector<double> SFD_profile; // Contains unique diameters from all datafiles - this is done in the init() function
std::vector<double> SFD_profile_Kessler; // Contains unique diameters from all datafiles - this is done in the init() function
std::vector<vector<double>> SFD_array;
std::vector<double> temp_diameters; // Mid-values of SFD diameters - FINISH THIS
std::vector<vector<double>> Integration_Setup_array;

std::string work_directory;

Constraints_Double My_Constraints;
Constraints_Double My_Measurements;

Outputs My_Outputs;

int output_switch = 0;

// Free parameters:
double As;
double Bs;
double alpha;
double beta;
double Dmid; // In micrometers
double N0;
double Dmax;
double comet_pericenter_power;
double particle_density; // In kg m^-3
double pericenter_index;

string radial_file = "Outputs_radial_output";
string ecl_file = "Outputs_ecliptic_profile_output";
const char * cradial_file = radial_file.c_str();
const char * cecl_file = ecl_file.c_str();
FILE * pFile;
FILE * pFile1;
FILE * pFile2;
FILE * pFile3;
FILE * pFile4;

FILE * pFile666;

// Used for testing
std::vector<Vec3> points;
Grid3D *peterito;
Octree *octree;
Octree *octree_iteration;
OctreePoint *octreePoints;
Vec3 qmin, qmax, r_scan;
Vec3 point;

vector<double> r_prof_init(10000); // Vector that keeps the radial profile for the 1st iteration
vector<double> r_prof_iter(10000); // Vector that keeps the radial profile for the 2nd iteration
vector<double> r_prof_brightness(10000); // Vector that keeps the radial profile for the 2nd iteration
vector<double> ecl_profile(180); // Vector that keeps the ecliptic latitude profile
vector<double> ecl_profile_mass(180); // Vector that keeps the ecliptic latitude profile




std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

void Initialize_1D_Array(std::vector<double> &init_array, int xdim){
	init_array.resize(xdim);
}

void Initialize_2D_Array(std::vector<vector<double>> &init_array, int xdim, int ydim){
	init_array.resize(xdim);
	for(int i=0; i<xdim; i++){
		init_array[i].resize(ydim);
	}
}

void Initialize_Outputs(){
	// CHANGE - unfinished and dirty - redo it
	int xdim = peterito->data.size();
	int ydim = peterito->data[0].size();

	Initialize_2D_Array(My_Outputs.disk_profile_face_on_number,xdim,ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_face_on_brightness,xdim,ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_face_on_impact_velocity,xdim,ydim);

	Initialize_2D_Array(My_Outputs.disk_profile_edge_on_number,xdim,ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_edge_on_brightness,xdim,ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_edge_on_impact_velocity,xdim,ydim);

	Initialize_2D_Array(My_Outputs.disk_profile_face_on_massloss,xdim,ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_edge_on_massloss,xdim,ydim);
	Initialize_1D_Array(My_Outputs.disk_profile_radial_massloss,10000);

	cout << "Check Output Size:" << std::endl;
 	cout << "Profile 1 size: " << My_Outputs.disk_profile_face_on_number.size() 			<< "\t" << My_Outputs.disk_profile_face_on_number[0].size()			<< std::endl;
	cout << "Profile 2 size: " << My_Outputs.disk_profile_face_on_brightness.size()		<< "\t" << My_Outputs.disk_profile_face_on_brightness[0].size()		<< std::endl;
	cout << "Profile 3 size: " << My_Outputs.disk_profile_face_on_impact_velocity.size()	<< "\t" << My_Outputs.disk_profile_face_on_impact_velocity[0].size()	<< std::endl;
	cout << "Profile 4 size: " << My_Outputs.disk_profile_edge_on_number.size()			<< "\t" << My_Outputs.disk_profile_edge_on_number[0].size()			<< std::endl;
	cout << "Profile 5 size: " << My_Outputs.disk_profile_edge_on_brightness.size()		<< "\t" << My_Outputs.disk_profile_edge_on_brightness[0].size()		<< std::endl;
	cout << "Profile 6 size: " << My_Outputs.disk_profile_edge_on_impact_velocity.size()	<< "\t" << My_Outputs.disk_profile_edge_on_impact_velocity[0].size()	<< std::endl;

}

void Initiate_Output_Files(){
	FILE * pFile_temp;
	pFile_temp = fopen("Outputs_Disk_Profile_Face_On_Number","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Disk_Profile_Face_On_Brightness","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Disk_Profile_Face_On_Impact_Velocity","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Disk_Profile_Edge_On_Number","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Disk_Profile_Edge_On_Brightness","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Disk_Profile_Edge_On_Impact_Velocity","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Mass_Accreted_Kessler","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Disk_Profile_Face_On_MassLoss","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Disk_Profile_Edge_On_MassLoss","w+"); fclose(pFile_temp);
	pFile_temp = fopen("Outputs_Disk_Profile_Radial_MassLoss","w+"); fclose(pFile_temp);
}

void Print_Outputs_Big(int file_counter, double refactor = 1.0){

std::vector<int> grid_index = peterito->getIndexes(point); //x,y, and z
int x = peterito->data.size();
int y = peterito->data[0].size();
int z = peterito->data[0][0].size();

FILE * pFile_temp;

pFile_temp = fopen("Outputs_Disk_Profile_Face_On_Number","a+");
for (int i=0;i<x; ++i){	for (int j=0;j<y; ++j){	if(My_Outputs.disk_profile_face_on_number[i][j] >0.0) {fprintf(pFile_temp,"%.5f %.5f %13.5e %d \n", dr[0]+i*2.0*dr[0]-octree_limit,dr[1]+j*2.0*dr[1]-octree_limit, My_Outputs.disk_profile_face_on_number[i][j]*refactor, file_counter );}	} };  fprintf(pFile_temp,"\n");	fclose(pFile_temp);

pFile_temp = fopen("Outputs_Disk_Profile_Face_On_Brightness","a+");
for (int i=0;i<x; ++i){	for (int j=0;j<y; ++j){	if(My_Outputs.disk_profile_face_on_brightness[i][j] >0.0) {fprintf(pFile_temp,"%.5f %.5f %13.5e %d \n", dr[0]+i*2.0*dr[0]-octree_limit,dr[1]+j*2.0*dr[1]-octree_limit, My_Outputs.disk_profile_face_on_brightness[i][j]*refactor, file_counter );}	} };  fprintf(pFile_temp,"\n");	fclose(pFile_temp);

pFile_temp = fopen("Outputs_Disk_Profile_Face_On_Impact_Velocity","a+");
for (int i=0;i<x; ++i){	for (int j=0;j<y; ++j){	if(My_Outputs.disk_profile_face_on_number[i][j] >0.0) {fprintf(pFile_temp,"%.5f %.5f %13.5e %d \n", dr[0]+i*2.0*dr[0]-octree_limit,dr[1]+j*2.0*dr[1]-octree_limit, My_Outputs.disk_profile_face_on_impact_velocity[i][j]*refactor, file_counter );}	} };  fprintf(pFile_temp,"\n");	fclose(pFile_temp);	

pFile_temp = fopen("Outputs_Disk_Profile_Edge_On_Number","a+");
for (int i=0;i<x; ++i){	for (int j=0;j<z; ++j){	if(My_Outputs.disk_profile_edge_on_number[i][j] >0.0) {fprintf(pFile_temp,"%.5f %.5f %13.5e %d \n", dr[0]+i*2.0*dr[0]-octree_limit,dr[1]+j*2.0*dr[1]-octree_limit, My_Outputs.disk_profile_edge_on_number[i][j]*refactor, file_counter );	} }; } fprintf(pFile_temp,"\n");	fclose(pFile_temp);

pFile_temp = fopen("Outputs_Disk_Profile_Edge_On_Brightness","a+");
for (int i=0;i<x; ++i){	for (int j=0;j<z; ++j){	if(My_Outputs.disk_profile_edge_on_brightness[i][j] >0.0) {fprintf(pFile_temp,"%.5f %.5f %13.5e %d \n", dr[0]+i*2.0*dr[0]-octree_limit,dr[1]+j*2.0*dr[1]-octree_limit, My_Outputs.disk_profile_edge_on_brightness[i][j]*refactor, file_counter );}	} };  fprintf(pFile_temp,"\n");	fclose(pFile_temp);

pFile_temp = fopen("Outputs_Disk_Profile_Edge_On_Impact_Velocity","a+");
for (int i=0;i<x; ++i){	for (int j=0;j<z; ++j){	if(My_Outputs.disk_profile_edge_on_number[i][j] >0.0) {fprintf(pFile_temp,"%.5f %.5f %13.5e %d \n", dr[0]+i*2.0*dr[0]-octree_limit,dr[1]+j*2.0*dr[1]-octree_limit, My_Outputs.disk_profile_edge_on_impact_velocity[i][j]*refactor, file_counter );}	} };  fprintf(pFile_temp,"\n");	fclose(pFile_temp);

/// 
pFile_temp = fopen("Outputs_Disk_Profile_Face_On_MassLoss","a+");
for (int i=0;i<x; ++i){	for (int j=0;j<y; ++j){	if(My_Outputs.disk_profile_face_on_massloss[i][j] >0.0) {fprintf(pFile_temp,"%.5f %.5f %13.5e %d \n", dr[0]+i*2.0*dr[0]-octree_limit,dr[1]+j*2.0*dr[1]-octree_limit, My_Outputs.disk_profile_face_on_massloss[i][j]*refactor, file_counter );}	} };  fprintf(pFile_temp,"\n");	fclose(pFile_temp);

pFile_temp = fopen("Outputs_Disk_Profile_Edge_On_MassLoss","a+");
for (int i=0;i<x; ++i){	for (int j=0;j<z; ++j){	if(My_Outputs.disk_profile_edge_on_massloss[i][j] >0.0) {fprintf(pFile_temp,"%.5f %.5f %13.5e %d \n", dr[0]+i*2.0*dr[0]-octree_limit,dr[1]+j*2.0*dr[1]-octree_limit, My_Outputs.disk_profile_edge_on_massloss[i][j]*refactor, file_counter );}	} };  fprintf(pFile_temp,"\n");	fclose(pFile_temp);

pFile_temp = fopen("Outputs_Disk_Profile_Radial_MassLoss","a+");
// for (int i=0;i<(int)My_Outputs.disk_profile_radial_massloss.size(); ++i){if(My_Outputs.disk_profile_radial_massloss[i] >0.0) {fprintf(pFile_temp,"%.5f %13.5e %d\n", i*dr.x*2.0+dr.x,My_Outputs.disk_profile_radial_massloss[i]*refactor, file_counter );}	} ;  fprintf(pFile_temp,"\n");	fclose(pFile_temp);
for (int i=0;i<(int)My_Outputs.disk_profile_radial_massloss.size(); ++i){if(My_Outputs.disk_profile_radial_massloss[i] >0.0) {fprintf(pFile_temp,"%.5f %13.5e %d\n", (i+0.5) * ( octree_limit) / (float)My_Outputs.disk_profile_radial_massloss.size() ,My_Outputs.disk_profile_radial_massloss[i]*refactor, file_counter );}	} ;  fprintf(pFile_temp,"\n");	fclose(pFile_temp);
}



void Print_Outputs_Small(int file_counter, double refactor = 1.0){

std::vector<int> grid_index = peterito->getIndexes(point); //x,y, and z
int x = peterito->data.size();
int y = peterito->data[0].size();
int z = peterito->data[0][0].size();

FILE * pFile_temp;

pFile_temp = fopen("Outputs_Mass_Accreted_Kessler","a+");
fprintf(pFile_temp,"%13.5e %d \n", mass_accreted_at_Earth*refactor, file_counter);fclose(pFile_temp);
}

void Zero_Outputs(){

int x = peterito->data.size();
int y = peterito->data[0].size();
int z = peterito->data[0][0].size();

for (int i=0;i<x; ++i){	for (int j=0;j<y; ++j){ My_Outputs.disk_profile_face_on_number[i][j] = 0.0;	}}
for (int i=0;i<x; ++i){	for (int j=0;j<y; ++j){ My_Outputs.disk_profile_face_on_brightness[i][j] = 0.0;	}}
for (int i=0;i<x; ++i){	for (int j=0;j<y; ++j){ My_Outputs.disk_profile_face_on_impact_velocity[i][j] = 0.0;}}
for (int i=0;i<x; ++i){	for (int j=0;j<z; ++j){ My_Outputs.disk_profile_edge_on_number[i][j] = 0.0;}}
for (int i=0;i<x; ++i){	for (int j=0;j<z; ++j){ My_Outputs.disk_profile_edge_on_brightness[i][j] = 0.0;}}
for (int i=0;i<x; ++i){	for (int j=0;j<z; ++j){ My_Outputs.disk_profile_edge_on_impact_velocity[i][j] = 0.0;}}

for (int i=0;i<x; ++i){	for (int j=0;j<y; ++j){ My_Outputs.disk_profile_face_on_massloss[i][j] = 0.0;}}
for (int i=0;i<x; ++i){	for (int j=0;j<z; ++j){ My_Outputs.disk_profile_edge_on_massloss[i][j] = 0.0;}}	
for (int i=0;i<(int)My_Outputs.disk_profile_radial_massloss.size(); ++i) {My_Outputs.disk_profile_radial_massloss[i] = 0.0;}

}

void record_disk_profiles(Vec3 point, float weight, float diameter, float collision_velocity){
std::vector<int> grid_index = peterito->getIndexes(point); //x,y, and z

int x = grid_index[0];
int y = grid_index[1];
int z = grid_index[2];

float heliocentric_distance = sqrt( point.x * point.x + point.y * point.y  + point.z * point.z);
float flux_at_distance = 1.0/(heliocentric_distance * heliocentric_distance);
float cross_section = pow(diameter,2.0)*0.25*M_PI*1e-12; // diameter is in micrometers

My_Outputs.disk_profile_face_on_number[x][y] += weight;
My_Outputs.disk_profile_edge_on_number[x][z] += weight;

My_Outputs.disk_profile_face_on_brightness[x][y] += weight * flux_at_distance * cross_section;
My_Outputs.disk_profile_edge_on_brightness[x][z] += weight * flux_at_distance * cross_section;

My_Outputs.disk_profile_face_on_impact_velocity[x][y] += weight * collision_velocity;
My_Outputs.disk_profile_edge_on_impact_velocity[x][z] += weight * collision_velocity;

}


void record_disk_mass_loss(Vec3 point, float mass_loss){
std::vector<int> grid_index = peterito->getIndexes(point); //x,y, and z

int x = grid_index[0];
int y = grid_index[1];
int z = grid_index[2];


float ecliptic_distance = sqrt( point.x * point.x + point.y * point.y);

My_Outputs.disk_profile_face_on_massloss[x][y] += mass_loss;
My_Outputs.disk_profile_edge_on_massloss[x][z] += mass_loss;

// int central_index = floor(octree_limit/dr.x * 0.5); // getting the ecliptic index

	int gridPoint_int = floor( ecliptic_distance / ( octree_limit) * (int)My_Outputs.disk_profile_radial_massloss.size() );
	if(gridPoint_int < (int)My_Outputs.disk_profile_radial_massloss.size() ) {My_Outputs.disk_profile_radial_massloss[gridPoint_int] += mass_loss;}
	// if(gridPoint_int < (int)My_Outputs.disk_profile_radial_massloss.size() && z == central_index ) {My_Outputs.disk_profile_radial_massloss[gridPoint_int] += mass_loss;}
	// We count only the points at the ecliptic; i.e. the central index
}

void radialProfile(std::vector<double> &r_prof, Vec3 point, float weight) {
	// Adds particle weight to the radial profile vector 
	float gridPoint_rad;
	int gridPoint_int;
	gridPoint_rad = sqrt( point.x * point.x + point.y * point.y );
	gridPoint_int = floor( gridPoint_rad / ( dr.x * 2.0 ) );
	if(gridPoint_int < (int)r_prof.size()) {r_prof[gridPoint_int] += weight;}
}

void radialProfileCross(std::vector<double> &r_prof, Vec3 point, float weight, float cross_section) {
	// Adds particle weight to the radial profile vector 
	float gridPoint_rad;
	int gridPoint_int;
	gridPoint_rad = sqrt( point.x * point.x + point.y * point.y );
	gridPoint_int = floor( gridPoint_rad / ( dr.x * 2.0 ) );
	if(gridPoint_int < (int)r_prof.size()) {r_prof[gridPoint_int] += weight*cross_section;}
}

void eclipticLatitudeProfile(std::vector<double> &ecl_prof, std::vector<double> &ecl_prof_mass, Vec3 point, float weight, float diameter) {
	// Ecliptic latitude profile beyond 1 au - we calculate the cross-section per degree of the ecliptic latitude + We also calculate the mass for the same profile
	// These two quantities are easy to compare to Nesvorny et al. (2011) etc.
	float ecliptic_latitude;
	int gridPoint_int;
	double r_eclipt = sqrt(point.x*point.x+point.y*point.y);
	double d_phi = 180.0/float(ecl_prof.size());

	if(r_eclipt > 1.0) {
		ecliptic_latitude = atan(point.z/(r_eclipt-1.0))*180.0/M_PI; // We are assuming that the particles in the cloud are not in 1:1 resonance with Earth and thus can be randomly oriented with respect to Earth
		gridPoint_int = floor( ( ecliptic_latitude+90.0 ) *d_phi );
		ecl_prof[gridPoint_int]      += weight*pow(diameter,2.0)*0.25*M_PI*1e-12; //Weight * Particle Cross-section - diameters are in micrometers so we have *1e-12
		ecl_prof_mass[gridPoint_int] += weight*pow(diameter,3.0)/6.0*M_PI*1e-18*particle_density; //Weight * Particle Cross-section - diameters are in micrometers so we have *1e-12
	}

}

int find_diameter_index(float diameter)	{	
	// Finds the diameter index in the SFD_array
	for(int i=0;i<(int)SFD_array.size();i++){
		if(float(SFD_array[i][0]) ==  diameter){
			return i;
		} 
	}
	cout << "Error finding dimameters, something is terribly wrong, this should not happen at all! --- function int find_diameter_index(float diameter)\n";
	exit(EXIT_FAILURE);
	return 999;
}

void SFD_1AU_Profile(std::vector<double> &SFD_profile, Vec3 point, float weight, float diameter) {
	// Ecliptic latitude profile beyond 1 au - we calculate the cross-section per degree of the ecliptic latitude + We also calculate the mass for the same profile
	// These two quantities are easy to compare to Nesvorny et al. (2011) etc.

		int dia = find_diameter_index(diameter);
		// if(point.norm() > 1.0-0.025 && point.norm() < 1.0+0.025) SFD_profile[dia] += weight;
		SFD_profile[dia] += weight;
}

void SFD_1AU_Profile_Kessler(std::vector<double> &SFD_profile, Vec3 point, float weight, float diameter) {
	// Ecliptic latitude profile beyond 1 au - we calculate the cross-section per degree of the ecliptic latitude + We also calculate the mass for the same profile
	// These two quantities are easy to compare to Nesvorny et al. (2011) etc.

		int dia = find_diameter_index(diameter);
		SFD_profile_Kessler[dia] += weight;
}

float dotProduct(Vec3 v1, Vec3 v2) { return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z; }

Vec3 crossProduct(Vec3 vector_a, Vec3 vector_b) { 
	Vec3 temp;
    temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
    temp[1] = vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0];
    temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
   return temp; }

void cosAngleVec(Vec3 v1, Vec3 v2, float &cos_angle, float &ratio){ 
	float v1_norm = v1.norm();
	float v2_norm = v2.norm();
	cos_angle = dotProduct(v1,v2)/v1_norm/v2_norm;
	ratio = v2_norm/v1_norm;
} 

void getDeltaV(Vec3 v1, Vec3 v2, float &deltaV) { Vec3 dVec = Vec3(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z); 	deltaV = dVec.norm(); } 

void readMyLine(string line, Vec3 &point, Vec3 &velvec1, int &parID, float &time_int){
	std::stringstream s(line);
	string getme;
	float x;
	float y;
	float z;
	float vx,vy,vz;

	while (getme.length() == 0) {getline(s, getme, ' ');}parID = stoi(getme)-1; getme = "";
	while (getme.length() == 0) {getline(s, getme, ' ');}time_int = stof(getme); getme = "";
	while (getme.length() == 0) {getline(s, getme, ' ');}x = stof(getme); getme = "";
	while (getme.length() == 0) {getline(s, getme, ' ');}y = stof(getme); getme = "";
	while (getme.length() == 0) {getline(s, getme, ' ');}z = stof(getme); getme = "";
	while (getme.length() == 0) {getline(s, getme, ' ');}vx = stof(getme); getme = "";
	while (getme.length() == 0) {getline(s, getme, ' ');}vy = stof(getme); getme = "";
	while (getme.length() == 0) {getline(s, getme, ' ');}vz = stof(getme); getme = "";

	point=Vec3(x,y,z);
	velvec1=Vec3(vx,vy,vz);
}

void readDatafile(string line, std::__cxx11::string &file_to_load, float &diameter, float &pericenter){
	std::stringstream s(line);
	string getme;
	pericenter = 999.999;
	
	getline(s, getme, '|');file_to_load = getme; getme = "";
	getline(s, getme, '|');diameter = stof(getme); getme = "";
	getline(s, getme, '|');if(getme.length()>0) {pericenter = stof(getme); getme = "";} // this is optional so we are outputting 999.999 instead
}

void readMyLineBinary(ifstream &file, Vec3 &point, Vec3 &velvec1, int &parID, float &time_int){
	char * memblock_id = new char [4];
	char * memblock_elms = new char [28];
	float* read_values; 

	file.read (memblock_id, 4);
	int* int_val = (int*)memblock_id;
	parID = int_val[0]-1;

	file.read (memblock_elms, 28);
	read_values = (float*)memblock_elms;//reinterpret as doubles

	time_int = read_values[0];
	point=Vec3(read_values[1],read_values[2],read_values[3]);
	velvec1=Vec3(read_values[4],read_values[5],read_values[6]);
	delete memblock_id;
	delete memblock_elms;
}

Vec3 gridVec3(Vec3 v, Vec3 dr) // Attaches the point to the 3d grid defined by dr
{ return Vec3(round(v.x*0.5/dr.x)*dr.x*2.0, round(v.y*0.5/dr.y)*dr.y*2.0, round(v.z*0.5/dr.z)*dr.z*2.0); }

void addOctreePoint(Octree *octree_insert, Vec3 point,float weight,int index,Vec3 velvec1){
		octreePoints = new OctreePoint;
		point = gridVec3(point,dr);
	    octreePoints->setPosition(point);
    	octreePoints->insertWeight(weight);
    	octreePoints->insertWeight_Iter(weight);
    	octreePoints->insertIndex(index);
    	octreePoints->insertVelocity(velvec1);
    	octreePoints->addRecord();
    	octree_insert->insert(octreePoints);
    	octreePoints = NULL;
}

void checkRadialProfile(FILE * pFile_w, std::vector<double> &r_prof,int counter, float reweighting_factor = 1.0){
	for (int i=0;i<(int)r_prof.size(); ++i){
		if(r_prof[i]>0) {fprintf(pFile_w,"%.5f %13.5e %d %13.5e\n", i*dr.x*2.0+dr.x, r_prof[i]*reweighting_factor, counter, r_prof_brightness[i]*reweighting_factor );}
	}
}

void checkEclipticProfile(FILE * pFile_w, std::vector<double> &ecl_prof, std::vector<double> &ecl_prof_mass, int counter, float reweighting_factor = 1.0){
	for (int i=0;i<(int)ecl_prof.size(); ++i){
		fprintf(pFile_w,"%.5f %13.5e %13.5e %d \n", float(i-89.5), ecl_prof[i]*reweighting_factor ,ecl_prof_mass[i]*reweighting_factor, counter );
	}
}

void checkSFDProfile(FILE * pFile_w, std::vector<double> &SFD_profile, int counter, float reweighting_factor = 1.0){
	for (int i=0;i<(int)SFD_profile.size(); ++i){
		fprintf(pFile_w,"%.5f %13.5e  %d  %13.5e\n", SFD_array[i][0], SFD_profile[i]*reweighting_factor, counter, SFD_profile_Kessler[i]*reweighting_factor );
	}
}

float getRadialDifference(std::vector<double> &r_prof1,std::vector<double> &r_prof2){
	double radial_difference = 0.0;
	double max_diff = 0.0;
	// for (int k=0;k<r_prof1.size() ;++k){ // The original version looks throughout the entire vector of radial distances, but the convergence in the bins at the origin takes long
	for (int k=3;k<(int)r_prof1.size() ;++k){
	 // radial_difference += pow(r_prof1[k]-r_prof2[k],2);
		if(r_prof1[k]>0.0 && r_prof2[k]>0.0){
			// printf("DEBUG INSIDE: iter: %d, radial_difference %.5e, temp_calc: %.5e, cntr: %.5e, array 1: %.5e, array 2: %.5e\n", k, radial_difference, r_prof2[k]/r_prof1[k]-1.0,cntr,r_prof1[k],r_prof2[k]);
			radial_difference = std::max(r_prof2[k]/r_prof1[k],r_prof1[k]/r_prof2[k]);
			if(max_diff<radial_difference) { max_diff=radial_difference;}
		}
	}
		// printf("DEBUG: radial_difference %.5e, cntr: %.5e\n", radial_difference,cntr);
		return max_diff-1.0;
}


void getCollProbWithEarth(Vec3 r_vec, Vec3 v_vec, Vec3 r_pl_vec, Vec3 v_pl_vec, float v_esc, float pl_cross, float timestep, float &a, float &e, float &inc, float &v_rel, float &prob, float &prob_unfoc) {
	// Using Kessler (1981) paper
	Vec3 H_vec;
	float H2;
	float HZ2;
	float cos2i;
	float sin2i;
	float energy;
	float GM = 1.32712440018e20;
	// float a;
	// float e;
	float sin2b;
	float peri;
	float apo;
	float cosgamma;
	float r_pl,v_pl,r,v;
	float v_par_target;
	float cosphi;
	// float v_rel;
	// float prob;
	// float prob_unfoc;
	float grav_foc;
	// float inc;

	r_pl = r_pl_vec.norm();
	v_pl = v_pl_vec.norm();
	r = r_vec.norm();
	v = v_vec.norm();
	

	//We assume that the planet is orbiting in the plane of reference, i.e. it has zero inclination that simplifies a lot of things
	sin2b = 0.0;

	//We assume that the planet is orbiting on a circular orbit -> cosgamma = 1

	H_vec = crossProduct(r_vec,v_vec);
	H2 = dotProduct(H_vec, H_vec);
	HZ2 = H_vec[2]*H_vec[2];

	cos2i = HZ2/H2;
	float cos_crit = 0.999996953828895;
	      // cos_crit = 0.996953828895;
		  // cos_crit = 0.9999999996953828895;
	if(cos2i > cos_crit) {cos2i = cos_crit;} // at least 0.1 degree
	sin2i = 1.0 - cos2i;


	energy =v*v*0.5-GM/r;
	
	if(energy<0){
		a = -0.5*GM/energy;
		e = sqrt(1.0 - H2/(GM*a));
		if(!isnormal(e)) {e = 0.0;}
	}

	peri = a*(1.0 - e);
	apo = a*(1.0 + e);

	v_par_target = sqrt(2.0*energy + 2.0*GM/(r_pl));
	cosgamma = sqrt(apo*peri/(r_pl* (2.0 * a - r_pl) ));
	cosphi = cosgamma * sqrt(cos2i);
	v_rel = sqrt(v_pl * v_pl + v_par_target * v_par_target - 2.0 * v_par_target * v_pl * cosphi);
	grav_foc = 1.0 + pow(v_esc/v_rel, 2.0);

	inc = acos(sqrt(cos2i))*180.0/M_PI;

	// printf("KESSLER TEST: v_rel = %13.5e, a = %13.5e, e = %13.5e, inc = %13.5e, v_pl = %13.5e, v_par = %13.5e, cosgamma = %13.5e, cosphi = %13.5e\n", v_rel, a/au, e, inc, v_pl, v_par_target, cosgamma, cosphi);
	// cin >> pinkevic;

// // TESTING
// 	prob = v_rel * pl_cross * grav_foc * timestep / (2.0 * pow(M_PI,3.0) * r_pl * a * sqrt( (sin2i - sin2b) *(r_pl - peri) * (apo - r_pl) ) );
// 	if(r_pl <= peri) {prob = 0.0;}
// 	if(r_pl >= apo) {prob = 0.0;}


// ORIGINAL
	// if(r_pl <= peri) {return 0.0;}
	// if(r_pl >= apo) {return 0.0;}
	prob = 0.0;
	prob_unfoc = 0.0;
	if(r_pl <= peri) {return;}
	if(r_pl >= apo) {return;}
	prob_unfoc = v_rel * pl_cross * timestep / (2.0 * pow(M_PI,3.0) * r_pl * a * sqrt( (sin2i - sin2b) *(r_pl - peri) * (apo - r_pl) ) );
	prob = prob_unfoc * grav_foc;

	// Error check
	if(!isnormal(prob)){ 
		prob = 0.0;
		prob_unfoc = 0.0;
		printf("Problem with the prob:\n"); 
		printf("Method check - a: %13.5e, e: %13.5e, inc: %13.5e, v_rel: %13.5e, grav_foc: %13.5e, energy: %13.5e, v_pl: %13.5e, v: %13.5e\n", a/1.5e11,e, acos(sqrt(cos2i))*180.0/M_PI,v_rel, grav_foc , energy, v_pl, v);
		printf("Method check - prob: %13.5e, cosphi: %13.5e, cosgamma: %13.5e, peri: %13.5e, rpl: %13.5e, apo: %13.5e, H2: %13.5e, HZ2: %13.5e,\n", prob, cosphi,cosgamma,peri,r_pl,apo, H2, HZ2);
		printf("Method check - x: %13.5e, y: %13.5e, z: %13.5e, vx: %13.5e, vy: %13.5e, vz: %13.5e, H2/GMa: %13.5e, \n", r_vec[0],r_vec[1],r_vec[2],v_vec[0],v_vec[1],v_vec[2],H2/(GM*a));
	return;
	}



	// prob can be really small in kms units, but it shouldn't go below 10e-38 in the solar system. Normally this value is above 10e-30 so we have 8 orders of magnitude cushion

	// printf("Method check - a: %13.5e, e: %13.5e, inc: %13.5e, v_rel: %13.5e, grav_foc: %13.5e, energy: %13.5e, v_pl: %13.5e, v: %13.5e\n", a/1.5e11,e, acos(sqrt(cos2i))*180.0/M_PI,v_rel, grav_foc , energy, v_pl, v);

	// TEST_probabilities.push_back(prob);

	return; 
}

std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}


   float compute_average(std::vector<double> &vi) {
     double sum = 0;
     for (double p:vi){sum = sum + p;}
     return (sum/double(vi.size()));
  }

   double compute_sum(std::vector<double> &vi) {
     double sum = 0;
     for (double p:vi){sum = sum + p;}
     return (sum);
  }

double broken_power_law(double alpha,double beta,double D,double Dmid,double Dmax,double N0){
	double a = alpha;
	double b = beta;
	if(a == 1.0) a += a+1e-6;
	if(b == 1.0) b += b+1e-6;
	double Na = N0*(a-1.0);
	double Nb = Na*pow(Dmax/Dmid,a-1.0);

	double Nout;
	if(D < Dmid){
		Nout = Nb/(b-1.0)*pow(Dmid/D,b-1.0)+Nb/(b-1.0)*(b-a)/(a-1.0);
		} 
	else {
		Nout = Nb/(a-1.0)*pow(Dmid/D,a-1.0);
	}
	return Nout; // CUMULATIVE NUMBER AT D

}


void calculate_SFD(double alpha,double beta,double Dmid,double Dmax,double N0){
	double average_log_spacing = 0.0;
	std::vector<float> differential_sizes;
	double tmp_number;
	double tmp_mass;

	if( (int)unique_diameters.size() == 1 ) 
	{
		SFD_array.push_back({unique_diameters[0],1.0,1.0});
		SFD_profile.resize(unique_diameters.size());
		SFD_profile_Kessler.resize(unique_diameters.size());
		return;
	}

	// GET THE AVERAGE SPACING BETWEEN DIAMETERS TO GET THE OFFSET FOR THE SMALLEST AND LARGEST DIAMETER - WE NEED THIS BECAUSE OUR DIAMETER IS A MID VALUE
	for(int i=1;i<(int)unique_diameters.size();i++){
		average_log_spacing += log10(unique_diameters[i])-log10(unique_diameters[i-1]);
	}
	average_log_spacing = average_log_spacing/float(unique_diameters.size()-1.0);


	// cout << "AVERAGE LOG SPACING \t" << average_log_spacing << "\n";

	 // temporary vector for the mid values
	temp_diameters.push_back(pow(10.0,log10(unique_diameters[0])-average_log_spacing*0.5)); // ADDING THE FIRST VALUE WHICH IS D-Dstep*0.5 in logscale

	// this feeds the mid-values vector
	for(int i=1;i<(int)unique_diameters.size();i++){
		tmp_number = log10(unique_diameters[i])*0.5+log10(unique_diameters[i-1])*0.5;
		temp_diameters.push_back(pow(10.0,tmp_number));
	}
	temp_diameters.push_back(pow(10.0,log10(unique_diameters[unique_diameters.size()-1])+average_log_spacing*0.5));


	// we have a cumulative distribution from the broken_power_law function so we calculate the differential between two cumulative values of N
	for(int i=0;i<(int)(temp_diameters.size()-1);i++){
		tmp_number = broken_power_law(alpha,beta,temp_diameters[i],Dmid,Dmax,N0) - broken_power_law(alpha,beta,temp_diameters[i+1],Dmid,Dmax,N0);
		tmp_mass =  pow(unique_diameters[i],3.0)/6.0*M_PI*particle_density*1e-18;
		differential_sizes.push_back(tmp_number);
		cout << i << "\t" << unique_diameters[i] << "\t" <<  tmp_number << "\t" << Integration_Setup_array[i][1] << "\n";
		SFD_array.push_back({unique_diameters[i],tmp_number/Integration_Setup_array[i][1],tmp_mass});
	   // cout << i << "\t" << unique_diameters[i] << "\t" <<  tmp_number << "\t" << temp_diameters[i] << "\t" << temp_diameters[i+1] << "\n";
	}

	SFD_profile.resize(unique_diameters.size());
	SFD_profile_Kessler.resize(unique_diameters.size());

}


void calculate_Earth_Mass_Flux(std::vector<double> &mass_vector,int iter_count){
	// CHECK FLUX AT EARTH
	Vec3 Earth_r;
	Vec3 Earth_v;
	std::vector<int> res_index;
	OctreePoint *results;

	float true_anom;
	float deltaV;
	int rec_iter;
	Vec3 velvec;



	float V0,mu,h;
	float a = 1.000e0*149597900e3; // AU in meters
	float e = 0.0e0;

	float r_pl, mu_pl, V_esc;
	float time_elapsed = 86400.0; // seconds per day
	float particle_number;
	float particle_mass;
	// float particle_density = 2000.0;
	double mass_accreted;
	float cross_section;
	

	int dia_index;


	r_pl = 6371e3;
	mu_pl = 3.986004418e14;
	V_esc = sqrt(2.0*mu_pl/r_pl);
	mu = 1.32712440018e20; // gravitational parameter of the Sun in SI units
	h = sqrt(mu*a*(1.0-e*e));
	V0 = mu/h; // Circular velocity in meters per second
	// cout << V0 << "\n";

	for(int i=0;i<1000;i++){ // steps in the true anomaly
		mass_accreted=0.0;
		true_anom=M_PI*i/500.0+1e-3; // we add a small delta so we do not start at the grid
		Earth_r.x=1.0*cos(true_anom);
		Earth_r.y=1.0*sin(true_anom);
		Earth_r.z=0.0;

		
		Earth_v.x=-V0*sin(true_anom);
		Earth_v.y=V0*cos(true_anom);
		Earth_v.z=0.0;
   		res_index = peterito->getIndexes(Earth_r); //
      results = &peterito->data[res_index[0]][res_index[1]][res_index[2]];

		if(results->getRecNum()>0){
    	rec_iter = results->getRecNum();
    	for (int k=0;k<rec_iter ;++k){

			 // dia_index = find_diameter_index(results->getDiameter(k));
    		dia_index = results->getIndex(k);

    		velvec = results->getVelocity(k) * auday2ms;;
    		particle_number = results->getWeight(k)*SFD_array[dia_index][1] / volume_unit;
    		// particle_mass = pow(results->getDiameter(k),3.0)/6.0*M_PI*rho_par*1e-18 ;
    		particle_mass = SFD_array[dia_index][2];
    		getDeltaV(velvec,Earth_v, deltaV);
    		cross_section = M_PI*pow(r_pl,2)*(1.0 + pow(V_esc/deltaV,2));
    		mass_accreted += particle_number * particle_mass * cross_section * time_elapsed * deltaV ;



    		// impactor_weight = results->getWeight(k) * SFD_array[dia_index][1] / volume_unit;
    		// // impactor_weight = results[0]->getWeight(k);
    		// tau_coll = impactor_weight * particle_cross_section * deltaV * t_record;

    		// mass_accreted += particle_number;
    		// printf("particle_number: %13.5e, diameter: %.5f, particle_mass %13.5e, cross_section %13.5e time_elapsed %13.5e, delta v %13.5e, volume_unit: %13.5e, mass accreted: %13.5e\n", particle_number, SFD_array[dia_index][0], particle_mass ,cross_section , time_elapsed, deltaV, volume_unit, particle_number * particle_mass * cross_section * time_elapsed * deltaV);

			}



    	}
    	// results.clear();
    	fprintf(pFile4,"%d %13.5e %d\n",i,mass_accreted,iter_count);
		mass_vector.push_back(mass_accreted);
	}

}


void load_2D_helper(std::string filename, std::vector<double> &x_vec, std::vector<double> &y_vec){
	std::ifstream file(filename);
	std::string str; 
    double x_axis;
    double y_axis;

    if(file.is_open()) {
    
    while (std::getline(file, str))
    {
    	std::stringstream s(str);
    	s >> x_axis >> y_axis; // Stringstream to float and float
		x_vec.push_back(x_axis);
		y_vec.push_back(y_axis);
    }

    } else {cout << "File not found! \t " << filename << "\n";}

}


void load_radar_constraints(){

	load_2D_helper(work_directory+"/./Data_Constraints/CMOR_AH_semi.txt",   My_Constraints.CMOR_AH_a_x,My_Constraints.CMOR_AH_a_y);
	load_2D_helper(work_directory+"/./Data_Constraints/CMOR_AH_inc.txt", My_Constraints.CMOR_AH_inc_x,My_Constraints.CMOR_AH_inc_y);
	load_2D_helper(work_directory+"/./Data_Constraints/CMOR_AH_ecc.txt", My_Constraints.CMOR_AH_ecc_x,My_Constraints.CMOR_AH_ecc_y);
	load_2D_helper(work_directory+"/./Data_Constraints/CMOR_AH_velocity.txt",  My_Constraints.CMOR_AH_vg_x,My_Constraints.CMOR_AH_vg_y);

	load_2D_helper(work_directory+"/./Data_Constraints/AMOR_HE_semi.txt",      My_Constraints.AMOR_AH_a_x  ,My_Constraints.AMOR_AH_a_y);
	load_2D_helper(work_directory+"/./Data_Constraints/AMOR_HE_inc.txt",       My_Constraints.AMOR_AH_inc_x,My_Constraints.AMOR_AH_inc_y);
	load_2D_helper(work_directory+"/./Data_Constraints/AMOR_HE_ecc.txt",       My_Constraints.AMOR_AH_ecc_x,My_Constraints.AMOR_AH_ecc_y);
	load_2D_helper(work_directory+"/./Data_Constraints/AMOR_HE_velocity.txt",  My_Constraints.AMOR_AH_vg_x ,My_Constraints.AMOR_AH_vg_y);

}

void copy_and_clear_vector(std::vector<double> &original_vector, std::vector<double> &zeroed_vector){
	double sum_of_elems = 0.0;


	if(zeroed_vector.size() < original_vector.size()){
		zeroed_vector.clear();
		for(int i=0; i<(int)original_vector.size();i++) {
			zeroed_vector.push_back(0.0);
		}
	}
	else {
		// printf("Putting Zeros in the vector - before: ");
			for(int j=0; j<(int)zeroed_vector.size(); j++) {
				sum_of_elems+=zeroed_vector[j];
			}
		// printf("%13.5e, after: ",sum_of_elems); sum_of_elems = 0.0;

		std::fill(zeroed_vector.begin(), zeroed_vector.end(), 0.0);

			for(int j=0; j<(int)zeroed_vector.size(); j++) {
				sum_of_elems+=zeroed_vector[j];
			}
		// printf("%13.5e \n",sum_of_elems); sum_of_elems = 0.0;
	}

	// std::copy(original_vector.begin(), original_vector.end(), zeroed_vector.begin()); std::fill(zeroed_vector.begin(), zeroed_vector.end(), 0.0);
}


void initialize_radar_measurements(){
	copy_and_clear_vector(My_Constraints.CMOR_AH_a_y,   My_Measurements.CMOR_AH_a_y);
	copy_and_clear_vector(My_Constraints.CMOR_AH_ecc_y, My_Measurements.CMOR_AH_ecc_y);
	copy_and_clear_vector(My_Constraints.CMOR_AH_inc_y, My_Measurements.CMOR_AH_inc_y);
	copy_and_clear_vector(My_Constraints.CMOR_AH_vg_y,  My_Measurements.CMOR_AH_vg_y);

	copy_and_clear_vector(My_Constraints.AMOR_AH_a_y,   My_Measurements.AMOR_AH_a_y);
	copy_and_clear_vector(My_Constraints.AMOR_AH_ecc_y, My_Measurements.AMOR_AH_ecc_y);
	copy_and_clear_vector(My_Constraints.AMOR_AH_inc_y, My_Measurements.AMOR_AH_inc_y);
	copy_and_clear_vector(My_Constraints.AMOR_AH_vg_y,  My_Measurements.AMOR_AH_vg_y);
}

// void record_radar_profile(float a, float e, float inc, float vg, float prob, float m_met, std::string radar_name){
void record_radar_profile(float a, float e, float inc, float vg, float prob, int diameter_index, std::string radar_name){

a = a/au;

float ion_velocity_comp;
float da = 0.06, de = 0.01, dinc = 1.8, dvg = 800.0;
float radar_ion_limit = 9e9;
float v_radar;
float v_esc = 11200.0;
int index_a, index_ecc, index_inc, index_vg;
// char tester;

float meteors_m_limit;
float meteors_d_limit;
float meteors_d_max;
float meteors_N_max;
float meteors_N_mid;
float meteors_N_record;
float meteors_recorded;

if(prob <1e-30) {return;}

// Radiant Position calculation assuming circular Earth
// See for example Pokorny & Vokrouhlicky 2013 or Vokrouhlicky, Pokorny, Nesvorny (2012)
float apl = 1.0;
float P_val = a*(1.0-e*e)/(apl);
float P_val_sqrt = sqrt (P_val);
// float sin_omega = sqrt(1.0 - pow((P_val-1.0) / e ,2));
// float v_apex = e*sin_omega / P_val_sqrt;
float v_apex = sqrt ( (e*e + 2.0*P_val - P_val*P_val - 1.0) / P_val );
float v_sun = 1.0 - P_val_sqrt*cos(inc*deg2rad);
float v_z = P_val_sqrt*sin(inc*deg2rad);
float v_rad = sqrt (v_apex*v_apex + v_sun*v_sun + v_z*v_z);


// There are 4 different configurations but we take only the positive ones, because they are symmetric about the origin
float ecc_lon = abs(atan2(v_apex,v_sun))*rad2deg;
float ecc_lat = asin(v_z/v_rad)*rad2deg;
// ecc_lat = atan(v_z/(sqrt(v_apex*v_apex + v_sun*v_sun)));
//
// printf("Ecc Lon: %10.3f, Ecc Lat: %10.3f, v_rad: %10.3f, v1: %10.3f, v2: %10.3f, v3: %10.3f  \n", ecc_lon, ecc_lat, v_rad*29.78, v_apex, v_sun, v_z  );
// printf("a: %10.3f,e: %10.3f,inc: %10.3f,vg: %10.3f  \n", a,e,inc,vg );


float HE_lon_center = 70.0;
float HE_lon_width = 15.0;
float HE_lat_center = 0.0;
float HE_lat_width = 10.0;

float HE_CMOR_lon_center = 67.0;
float HE_CMOR_lon_width  = 12.0;
float HE_CMOR_lat_center = 2.4;
float HE_CMOR_lat_width  = 10.0;

float HE_AMOR_lon_center = 65.0;
float HE_AMOR_lon_width  = 17.0;
float HE_AMOR_lat_center = -15.0;
float HE_AMOR_lat_width  = 23.0;

float NA_lon_center = 0.0;
float NA_lon_width = 20.0;
float NA_lat_center = 18.0;
float NA_lat_width = 18.0;

float NT_lon_center = 0.0;
float NT_lon_width = 15.0;
float NT_lat_center = 55.0;
float NT_lat_width = 10.0;

std::string radiant_source_name = "OTHER";
if( abs(ecc_lon - HE_lon_center) < HE_lon_width && abs(ecc_lat - HE_lat_center) < HE_lat_width) { radiant_source_name = "HELION";}
if( abs(ecc_lon - NA_lon_center) < NA_lon_width && abs(ecc_lat - NA_lat_center) < NA_lat_width) { radiant_source_name = "NORTH APEX";}
if( abs(ecc_lon - NT_lon_center) < NT_lon_width && abs(ecc_lat - NT_lat_center) < NT_lat_width) { radiant_source_name = "NORTH TOROIDAL";}

// printf("Radar test: ecc_lon: %13.5f,ecc_lat %13.5f, source: ", ecc_lon, ecc_lat); cout << radiant_source_name << std::endl;
v_radar = sqrt(v_esc*v_esc + vg*vg);

ion_velocity_comp = 1e7*pow((v_radar/30000.0),3.5);
// ion = m_met*1e7*pow((v_radar/30000.0),3.5);

index_a = floor(a/da);
index_ecc = floor(e/de);
index_inc = floor(inc/dinc);
index_vg = floor(vg/dvg);

// printf("Radar profile index:  %d %d %d %d %13.5e\n", index_a, index_ecc, index_inc,index_vg, ion);
// printf("Radar profile orbs:  %13.5e %13.5e %13.5e %13.5e %13.5e\n", a/au, e, inc, vg, v_radar);
// cin >> tester;

// CHANGE - QUICK HACK
radiant_source_name = "OTHER";
if( abs(ecc_lon - HE_CMOR_lon_center) < HE_CMOR_lon_width && abs(ecc_lat - HE_CMOR_lat_center) < HE_CMOR_lat_width) { radiant_source_name = "HELION";}

if(radar_name == "CMOR") {
	radar_ion_limit = 1.0;
	meteors_m_limit = radar_ion_limit/ion_velocity_comp;
	meteors_d_limit = pow(meteors_m_limit*6.0/(M_PI*particle_density), 1.0/3.0)*1e6;
	meteors_d_max = temp_diameters[diameter_index+1];
	meteors_N_max = broken_power_law(alpha,beta,meteors_d_max,Dmid,Dmax,N0);
	meteors_N_mid = broken_power_law(alpha,beta,meteors_d_limit,Dmid,Dmax,N0);
	meteors_N_record = ( meteors_N_mid - meteors_N_max ) / Integration_Setup_array[diameter_index][1];
	if(meteors_N_record > SFD_array[diameter_index][1]) {meteors_N_record = SFD_array[diameter_index][1];} // We don't want more meteors recorded than the maximum in the bin
	if(meteors_N_record < 0.0) {meteors_N_record = 0.0;} // We don't want to record negative numbers

	meteors_recorded = meteors_N_record*prob;
	// meteors_N_mid = broken_power_law(alpha,beta,temp_diameters[diameter_index],Dmid,Dmax,N0)/Integration_Setup_array[diameter_index][1];  

	// printf("Power law test: alpha: %.f, beta: %.f, Dmid: %.f, Dmax: %.f, N0: %.f\n", alpha,beta,Dmid,Dmax,N0);
	
	// printf("Radar Test: N_bin: %13.5e, d_max: %13.5e, d_mid: %13.5e, v_radar: %13.5e, N_max: %13.5e, N_mid: %13.5e, N_record: %13.5e\n", SFD_array[diameter_index][1],meteors_d_max, meteors_d_limit, v_radar, meteors_N_max, meteors_N_mid, meteors_N_record);



	// printf("Radar Test 2: N_bin: %13.5e, d_max: %13.5e, d_mid: %13.5e, v_radar: %13.5e, N_max: %13.5e, N_mid: %13.5e, N_record: %13.5e\n", SFD_array[diameter_index][1],meteors_d_max, meteors_d_limit, v_radar, meteors_N_max, meteors_N_mid, meteors_N_record);

	// cin >> tester;


	// ion = ion_velocity_comp*SFD_array[diameter_index][2];

	if(radiant_source_name == "HELION") {
		// printf("Radar test CMOR: ecc_lon: %13.5f,ecc_lat %13.5f, source: ", ecc_lon, ecc_lat); cout << radiant_source_name << std::endl;
		if(index_a < (int)My_Measurements.CMOR_AH_a_y.size()) My_Measurements.CMOR_AH_a_y[index_a] += meteors_recorded;
		if(index_ecc < (int)My_Measurements.CMOR_AH_ecc_y.size()) My_Measurements.CMOR_AH_ecc_y[index_ecc] += meteors_recorded;
		if(index_inc < (int)My_Measurements.CMOR_AH_inc_y.size()) My_Measurements.CMOR_AH_inc_y[index_inc] += meteors_recorded;
		if(index_vg < (int)My_Measurements.CMOR_AH_vg_y.size()) My_Measurements.CMOR_AH_vg_y[index_vg] += meteors_recorded;
		}
	}


// CHANGE - QUICK HACK
radiant_source_name = "OTHER";
if( abs(ecc_lon - HE_AMOR_lon_center) < HE_AMOR_lon_width && abs(ecc_lat - HE_AMOR_lat_center) < HE_AMOR_lat_width) { radiant_source_name = "HELION";}

if(radar_name == "AMOR") {
	radar_ion_limit = 0.003;
	meteors_m_limit = radar_ion_limit/ion_velocity_comp;
	meteors_d_limit = pow(meteors_m_limit*6.0/(M_PI*particle_density), 1.0/3.0)*1e6;
	meteors_d_max = temp_diameters[diameter_index+1];
	meteors_N_max = broken_power_law(alpha,beta,meteors_d_max,Dmid,Dmax,N0);
	meteors_N_mid = broken_power_law(alpha,beta,meteors_d_limit,Dmid,Dmax,N0);
	meteors_N_record = ( meteors_N_mid - meteors_N_max ) / Integration_Setup_array[diameter_index][1];
	if(meteors_N_record > SFD_array[diameter_index][1]) {meteors_N_record = SFD_array[diameter_index][1];} // We don't want more meteors recorded than the maximum in the bin
	if(meteors_N_record < 0.0) {meteors_N_record = 0.0;} // We don't want to record negative numbers

	meteors_recorded = meteors_N_record*prob;

	if(radiant_source_name == "HELION") {
		// printf("Radar test AMOR: ecc_lon: %13.5f,ecc_lat %13.5f, source: ", ecc_lon, ecc_lat); cout << radiant_source_name << std::endl;
		if(index_a < (int)My_Measurements.AMOR_AH_a_y.size()) My_Measurements.AMOR_AH_a_y[index_a] += meteors_recorded;
		if(index_ecc < (int)My_Measurements.AMOR_AH_ecc_y.size()) My_Measurements.AMOR_AH_ecc_y[index_ecc] += meteors_recorded;
		if(index_inc < (int)My_Measurements.AMOR_AH_inc_y.size()) My_Measurements.AMOR_AH_inc_y[index_inc] += meteors_recorded;
		if(index_vg < (int)My_Measurements.AMOR_AH_vg_y.size()) My_Measurements.AMOR_AH_vg_y[index_vg] += meteors_recorded;
		}
	}

if(radar_ion_limit > 1e9) {cout << "Error: Incorrect Radar Selected!" << std::endl;}
}

void normalize_vector(std::vector<double> &vector){

	double sum_of_elems = 0.0;

	for(int i=0; i<(int)vector.size(); i++) {
		sum_of_elems+=vector[i];
	}
	
	for(int i=0; i<(int)vector.size(); i++) {
		vector[i] = vector[i]/sum_of_elems;
	}
}

void normalize_radar_measurements(){
	normalize_vector(My_Measurements.CMOR_AH_a_y);
	normalize_vector(My_Measurements.CMOR_AH_ecc_y);
	normalize_vector(My_Measurements.CMOR_AH_inc_y);
	normalize_vector(My_Measurements.CMOR_AH_vg_y);

	normalize_vector(My_Measurements.AMOR_AH_a_y);
	normalize_vector(My_Measurements.AMOR_AH_ecc_y);
	normalize_vector(My_Measurements.AMOR_AH_inc_y);
	normalize_vector(My_Measurements.AMOR_AH_vg_y);
}

void write_radar_profile(FILE * pFile_w, int counter, std::string radar_name){
	if(radar_name == "CMOR") {
	for (int i=0;i<(int)My_Measurements.CMOR_AH_a_y.size(); ++i){
		fprintf(pFile_w,"%.5f %13.5e %.5f %13.5e %.5f %13.5e %.5f %13.5e %d \n", 
			My_Constraints.CMOR_AH_a_x[i], My_Measurements.CMOR_AH_a_y[i],
			My_Constraints.CMOR_AH_ecc_x[i], My_Measurements.CMOR_AH_ecc_y[i],
			My_Constraints.CMOR_AH_inc_x[i], My_Measurements.CMOR_AH_inc_y[i],
			My_Constraints.CMOR_AH_vg_x[i], My_Measurements.CMOR_AH_vg_y[i],
			counter );
	}
	}
	if(radar_name == "AMOR") {
	for (int i=0;i<(int)My_Measurements.AMOR_AH_a_y.size(); ++i){
		fprintf(pFile_w,"%.5f %13.5e %.5f %13.5e %.5f %13.5e %.5f %13.5e %d \n", 
			My_Constraints.AMOR_AH_a_x[i],   My_Measurements.AMOR_AH_a_y[i],
			My_Constraints.AMOR_AH_ecc_x[i], My_Measurements.AMOR_AH_ecc_y[i],
			My_Constraints.AMOR_AH_inc_x[i], My_Measurements.AMOR_AH_inc_y[i],
			My_Constraints.AMOR_AH_vg_x[i],  My_Measurements.AMOR_AH_vg_y[i],
			counter );
	}
	}
}

double radar_profile_RMS(){

	double RMS_a_CMOR = 0.0;
	double RMS_e_CMOR = 0.0;
	double RMS_i_CMOR = 0.0;
	double RMS_v_CMOR = 0.0;

	double RMS_a_AMOR = 0.0;
	double RMS_e_AMOR = 0.0;
	double RMS_i_AMOR = 0.0;
	double RMS_v_AMOR = 0.0;

	double RMS_total = 0.0;

	for (int i=0;i<(int)My_Measurements.CMOR_AH_a_y.size(); ++i){
			RMS_a_CMOR += pow (My_Constraints.CMOR_AH_a_y[i] -   My_Measurements.CMOR_AH_a_y[i], 2.0);
			RMS_e_CMOR += pow (My_Constraints.CMOR_AH_ecc_y[i] - My_Measurements.CMOR_AH_ecc_y[i], 2.0);
			RMS_i_CMOR += pow (My_Constraints.CMOR_AH_inc_y[i] - My_Measurements.CMOR_AH_inc_y[i], 2.0);
			RMS_v_CMOR += pow (My_Constraints.CMOR_AH_vg_y[i] -  My_Measurements.CMOR_AH_vg_y[i], 2.0);
	}

	for (int i=0;i<(int)My_Measurements.AMOR_AH_a_y.size(); ++i){
			RMS_a_AMOR += pow (My_Constraints.AMOR_AH_a_y[i] -   My_Measurements.AMOR_AH_a_y[i], 2.0);
			RMS_e_AMOR += pow (My_Constraints.AMOR_AH_ecc_y[i] - My_Measurements.AMOR_AH_ecc_y[i], 2.0);
			RMS_i_AMOR += pow (My_Constraints.AMOR_AH_inc_y[i] - My_Measurements.AMOR_AH_inc_y[i], 2.0);
			RMS_v_AMOR += pow (My_Constraints.CMOR_AH_vg_y[i] -  My_Measurements.AMOR_AH_vg_y[i], 2.0);
	}

	RMS_total += RMS_a_CMOR;
	RMS_total += RMS_e_CMOR;
	RMS_total += RMS_i_CMOR;
	RMS_total += RMS_v_CMOR;

	RMS_total += RMS_a_AMOR;
	RMS_total += RMS_e_AMOR;
	RMS_total += RMS_i_AMOR;
	RMS_total += RMS_v_AMOR;


 return RMS_total;
}

double calculate_cloud_mass(){
	double cloud_total_mass = 0.0;
	double record_total_mass;
	int record_index;
	int rec_iter;
	OctreePoint *results;


	// cout << "Test Cloud" << peterito->data.size() << "\t" << peterito->data[0].size() << "\t" << peterito->data[0][0].size() << "\n";

		for(int i=0; i<(int)peterito->data.size(); i++){
					for(int j=0; j<(int)peterito->data[0].size(); j++){
						for(int k=0; k<(int)peterito->data[0][0].size(); k++){

			results = &peterito->data[i][j][k];

		    if(results->getRecNum()>0){

    			rec_iter = results->getRecNum();
    			for (int l=0;l<rec_iter ;++l){
    		
    				record_index = results->getIndex(l);
    				record_total_mass = results->getWeight_Iter(l) * SFD_array[record_index][1] * SFD_array[record_index][2];
    				cloud_total_mass += record_total_mass;

		}
	}
}
}
}
return double(cloud_total_mass);
	
}

double calculate_cloud_area(){
	double cloud_total_area = 0.0;
	double record_total_area;
	int record_index;
	int rec_iter;
	float diameter;
	OctreePoint *results;


	// cout << "Test Cloud" << peterito->data.size() << "\t" << peterito->data[0].size() << "\t" << peterito->data[0][0].size() << "\n";

		for(int i=0; i<(int)peterito->data.size(); i++){
					for(int j=0; j<(int)peterito->data[0].size(); j++){
						for(int k=0; k<(int)peterito->data[0][0].size(); k++){

			results = &peterito->data[i][j][k];

		    if(results->getRecNum()>0){

    			rec_iter = results->getRecNum();
    			for (int l=0;l<rec_iter ;++l){
    		
    				record_index = results->getIndex(l);
    				diameter = SFD_array[record_index][0];
    				record_total_area = results->getWeight_Iter(l) * SFD_array[record_index][1] *pow(diameter,2.0)*0.25*M_PI*1e-12;

    				cloud_total_area += record_total_area;

		}
	}
}
}
}
return double(cloud_total_area);
	
}


void pre_init(){

	float diameter;

	int parID;
	float time_int;
	float pericenter;
	Vec3 velvec1,point;

	Vec3 vec1,vec2;

	const char * cver_datafile;

  string line;
  string file_to_load;


  
  ifstream datafile (work_directory+"/./COLL_datafile");
  if (datafile.is_open())
  {
    while ( getline(datafile,line,'\n') )
    {
  	readDatafile(line,file_to_load,diameter,pericenter);
  
    unique_diameters.push_back(diameter);
    Integration_Setup_array.push_back({diameter,0.0});

      file_to_load=ReplaceAll(file_to_load,std::string(" "), std::string(""));
      cver_datafile = file_to_load.c_str();

  ifstream myfile(cver_datafile, ios::in|ios::binary| ios::ate);

  
  if (myfile.is_open())
  {
    myfile.seekg (0, ios::beg);
		while (!myfile.eof())
		{
    
    readMyLineBinary(myfile, point, velvec1,parID, time_int ); // reads one line from the input file
    if(myfile.eof()) break;

    if(time_int < 0.1) {Integration_Setup_array[Integration_Setup_array.size()-1][1]+=1.0;}
    if(time_int > 0.1) {break;}

   



	
	}

    myfile.close();
  }
  else cout << "Unable to open binary file\n"; 

  // octree->initialize_Second_weights();

      }
    datafile.close();

  }
  else cout << "Unable to open file - COLL_datafile\n" ;



	// Create a vector of unique diameters for the SFD calculation
  std::sort (unique_diameters.begin(), unique_diameters.end());  // First we sort it because we do not require sorted inputs
  std::vector<float>::iterator it; // Variable containting the number of unique values
  it = std::unique (unique_diameters.begin(), unique_diameters.end());   // This finds unique values in the vector and sets "it" to their count
  unique_diameters.resize( std::distance(unique_diameters.begin(),it) ); // Then resize the original vector so we carve out empty space - We are done here


  std::vector<float> tmp_array;
  tmp_array.resize(unique_diameters.size(),0.0);

  for(int k2=0;k2<(int)Integration_Setup_array.size();k2++){
  	for(int k3=0;k3<(int)unique_diameters.size();k3++){
  		if(float(unique_diameters[k3]) == float(Integration_Setup_array[k2][0])) { tmp_array[k3] += Integration_Setup_array[k2][1];};
  	}
  }
  
  Integration_Setup_array.clear();
  for(int k3=0;k3<(int)unique_diameters.size();k3++){
  Integration_Setup_array.push_back({unique_diameters[k3],tmp_array[k3]});
  cout << "Diameter check: " << unique_diameters[k3] << "\t" <<tmp_array[k3] << "\n" ;
}



  return;



}



void init() {

	std::vector<double> particle_weight;
	float current_integration_time = 9e99;// Make this number high so it always higher than the first line 
	// Meteoroid orbital elements
	float a_met, e_met, inc_met, vg_met;
	double test_zody_mass = 0.0, test_zody_area = 0.0;
	float deltaV;

	int rec_iter,breakable;
	
	float weight;
	float cosphi;
	float ratio;
	float diameter;

	int parID;
	int dia_index;
	float time_int;
	float pericenter;
	Vec3 velvec1,point;

	std::vector<int> res_index;
	OctreePoint *results;

	total_number_of_particles_integration = 0;
	int file_cntr=0;	

	Vec3 vec1,vec2;

	const char * cver_datafile;



  string line;
  string file_to_load;
 
 // int line_width = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
std::string spaces{80 - 1, ' '};

printf("The parameters are: angle difference = %.5f (%.5f deg), velocity ratio difference = %.5f, initial weights = %.5e\n", angle_diff, acos(angle_diff)*180.0/M_PI,ratio_diff,initial_weight); fflush(stdout);

  ifstream datafile (work_directory+"/./COLL_datafile");
  if (datafile.is_open())
  {
    while ( getline(datafile,line,'\n') )
    {
    readDatafile(line,file_to_load,diameter,pericenter);
    cout << "Processing - " << file_to_load << ", Diameter: " << diameter << ", Pericenter: " << pericenter <<  ", ";
    // unique_diameters.push_back(diameter);
    // Integration_Setup_array.push_back({diameter,0.0});
    dia_index = find_diameter_index(diameter);
    weight = initial_weight * pow(pericenter, pericenter_index);

	


      file_to_load=ReplaceAll(file_to_load,std::string(" "), std::string(""));
      cver_datafile = file_to_load.c_str();
      // printf("|%s|\n|%s|",cdata_file,cver_datafile);fflush(stdout);
      // file_to_load.replace( file_to_load.begin(), file_to_load.end(), '\\', '\\\\'); //replace all occurances of 'x' with 'y'
// file_to_load=ReplaceAll(file_to_load,std::string("/"), std::string("\\\\"));
// cout << "Processing - " << file_to_load << ", Diameter: " << diameter << "\n";
  ifstream myfile(cver_datafile, ios::in|ios::binary| ios::ate);

  // ifstream myfile ("output.small");


  // diameter = 50.0;

  
  if (myfile.is_open())
  {
    myfile.seekg (0, ios::beg);
		while (!myfile.eof())
		{
    
    readMyLineBinary(myfile, point, velvec1,parID, time_int ); // reads one line from the input file
    if(myfile.eof()) break;


    // CHANGE Mar 23
    if(myfile.eof()) break;
    if(current_integration_time > time_int) { particle_weight.clear(); current_integration_time = time_int;}

    if(current_integration_time == time_int) // Here we start reading a new file so we start a new vector etc.
    {
    	particle_weight.push_back(weight);
    	// printf("Pushed weight - ID: %d, weight: %13.5e, init weight: %13.5e, vector lenght: %lu \n",parID,particle_weight[parID],initial_weight,particle_weight.size()); fflush(stdout); 
    	while((int)particle_weight.size()<parID)    	{printf("Particle %lu mysteriously missing at the very beginning of the datafile => next particle is %d \n", particle_weight.size(), parID);  particle_weight.push_back(weight);}
    }
    // ----

    // if(time_int < 0.1) {Integration_Setup_array[Integration_Setup_array.size()-1][1]+=1.0;}

    // radialProfile(r_prof_init, point,1.0); // adds the read line to the radial profile "r_prof"

    // printf("TEST: %.5f %.5f %.5f %.5f %.5f %.5f\n", point.x,point.y,point.z,velvec1.x,velvec1.y,velvec1.z); fflush(stdout);
    // cin >> pink;



    if(point.maxComponent()>=octree_limit || point.minComponent()<=-octree_limit) { continue;} // removes the points beyond the octree dimension 
	// std::cout << point.x << "\t" << point.y << "\t" << point.z << "\n";
    res_index = peterito->getIndexes(point); //
    for (int jj=0; jj<3;jj++){
    if(res_index[jj]>=400 || res_index[jj]<0) cout << "index:" << jj<< "\t is "<<res_index[jj]<< "\t" << point.x << "\t" << point.y << "\t" << point.z << "\n"; 
    }
    
    
    // octree->getPointsInsideBox({point.x-dr.x,point.y-dr.y,point.z-dr.z}, {point.x+dr.x,point.y+dr.y,point.z+dr.z}, results); //
    // cout << "First Test" << peterito->data[res_index[0]][res_index[1]][res_index[2]].getRecNum() << "\n";
    results = &peterito->data[res_index[0]][res_index[1]][res_index[2]];

    if(results->getRecNum()>0){
    	// cout << results->getRecNum() << '\n'; fflush(stdout);
    	rec_iter = results->getRecNum();
    	breakable = 1;

    	
    	
    	
    	

    	for (int k=0;k<rec_iter ;++k){
    		vec1 = velvec1;
    		vec2 = results->getVelocity(k);
    		getDeltaV(vec1,vec2, deltaV);
    		cosAngleVec(velvec1,vec2, cosphi, ratio);

    		if(cosphi > angle_diff && abs(ratio-1)< ratio_diff && (results->getIndex(k) == dia_index))
    		{
    			results->setWeight(results->getWeight(k)+weight,k);
    			results->setWeight_Iter(results->getWeight(k)+weight,k);
    			breakable = 0;
    			break;
    		}

    	}
    		if(breakable == 1)
    		{
    		results->insertWeight(weight);
    		results->insertWeight_Iter(weight);
    		results->insertIndex(dia_index);
    		results->insertVelocity(velvec1);
    		results->addRecord();
    		}
 

    } else{

    	// addOctreePoint(octree, point,weight,dia_index,velvec1);
    	// cout << "Inserted Data Before" << results->getRecNum() << "\n";
    		results->insertWeight(weight);
    		results->insertWeight_Iter(weight);
    		results->insertIndex(dia_index);
    		results->insertVelocity(velvec1);
    		results->addRecord();
		// cout << "Inserted Data After" << results->getRecNum() << "\n";


    }

    // ORIGINAL
    // radialProfile(r_prof_init, point, weight*SFD_array[dia_index][1]); // adds the read line to the radial profile "r_prof"
    // eclipticLatitudeProfile(ecl_profile, ecl_profile_mass, point, weight*SFD_array[dia_index][1], diameter);
    // SFD_1AU_Profile(SFD_profile,point,weight*SFD_array[dia_index][1], diameter);
    // ---------
    float impact_prob_with_earth;
    float impact_prob_with_probe;
	getCollProbWithEarth(point*au, velvec1*auday2ms, Vec3(1.5e11,0,0), Vec3(0,-29800.0,0), 11200e0, M_PI*6371e3*6371e3, 86400.0, a_met, e_met, inc_met, vg_met,impact_prob_with_earth,impact_prob_with_probe); // Probability per day
    mass_accreted_at_Earth += impact_prob_with_earth * particle_weight[parID] * SFD_array[dia_index][1] * SFD_array[dia_index][2];
    record_radar_profile(a_met,e_met,inc_met,vg_met, impact_prob_with_earth * particle_weight[parID], dia_index, "CMOR");
    record_radar_profile(a_met,e_met,inc_met,vg_met, impact_prob_with_earth * particle_weight[parID], dia_index, "AMOR");
	if(output_switch == 1) {record_disk_profiles(point, particle_weight[parID]*SFD_array[dia_index][1], diameter, deltaV*auday2ms);}

    radialProfile(r_prof_iter, point, particle_weight[parID]*SFD_array[dia_index][1] ); // adds the read line to the radial profile "r_prof"
    radialProfileCross(r_prof_brightness, point, particle_weight[parID]*SFD_array[dia_index][1],pow(diameter,2.0)*0.25*M_PI*1e-12 ); // adds the read line to the radial profile "r_prof"
    eclipticLatitudeProfile(ecl_profile, ecl_profile_mass, point, particle_weight[parID]*SFD_array[dia_index][1], diameter);
    SFD_1AU_Profile(SFD_profile,point,impact_prob_with_probe*particle_weight[parID]*SFD_array[dia_index][1], diameter);
    SFD_1AU_Profile_Kessler(SFD_profile_Kessler,point,impact_prob_with_earth*particle_weight[parID]*SFD_array[dia_index][1], diameter);

    test_zody_area += particle_weight[parID]*SFD_array[dia_index][1]*pow(diameter,2.0)*0.25*M_PI*1e-12;
    test_zody_mass += particle_weight[parID]*SFD_array[dia_index][1]*pow(diameter,3.0)/6.0*M_PI*1e-18*particle_density; 
    total_number_of_particles_integration+=1;
    // results.clear();

	
	}

    myfile.close();
    printf("Inserted %d points to octree\n",total_number_of_particles_integration); fflush(stdout);
    file_cntr++;
  }
  else cout << "Unable to open binary file\n"; 

  // octree->initialize_Second_weights();

      }
    datafile.close();

  }
  else cout << "Unable to open file\n" ;



  printf("Inserted Total of %d points to octree\n",total_number_of_particles_integration); fflush(stdout);


	double initial_cloud_area = calculate_cloud_area();
	double init_reweight = intended_cloud_area/initial_cloud_area;
	cout << "Total Cloud Area: " << initial_cloud_area << std::endl;
	cout << "Reweighting the area of the cloud by a factor " << init_reweight << "New weight:" <<  initial_weight * init_reweight << "\n";
	initial_weight = initial_weight * init_reweight;




  checkRadialProfile(pFile, r_prof_iter, 0, init_reweight);fprintf(pFile,"\n");
  checkEclipticProfile(pFile1, ecl_profile, ecl_profile_mass, 0,init_reweight);fprintf(pFile1,"\n");
  checkSFDProfile(pFile2, SFD_profile, 0,init_reweight);fprintf(pFile2,"\n");


  	write_radar_profile(pFile3, 0,"CMOR");fprintf(pFile3,"\n");
	write_radar_profile(pFile4, 0,"AMOR");fprintf(pFile4,"\n");
	std::copy(r_prof_iter.begin(), r_prof_iter.end(), r_prof_init.begin()); std::fill(r_prof_iter.begin(), r_prof_iter.end(), 0); // Copy the r_prof_iter to r_prof_init and put zeros to r_prof_iters
	if(output_switch == 1) {Print_Outputs_Big(0,init_reweight);}
	Print_Outputs_Small(0,init_reweight);
	Zero_Outputs();

  return;
}


/// Iteration --- CONTINUE HERE - make a vector for particle IDS and their weights and keep information about the time - if there is a discontinuity then renew the ID vector
// void iteration(Octree *octree_prev,Octree *octree_next) {
void iteration(const double radial_difference) {

	std::vector<double> particle_weight;
	
	float setting_weight;
	int rec_iter,breakable;
	
	
	// Meteoroid orbital elements
	float a_met, e_met, inc_met, vg_met;

	// float weight;
	float cosphi;
	float ratio;
	float diameter;
	float old_diameter;
	float weight;
	int parID;
	float time_int;
	float current_integration_time = 9e99;// Make this number high so it always higher than the first line 
	float deltaV;
	float impactor_weight;
	float particle_cross_section;
	float particle_mass;
	float t_record;
	float tau_coll;

	double test_zody_mass = 0.0, test_zody_area = 0.0;



	

	Vec3 velvec1,point;

  char debug_char;
  
  string file_to_load;

  const char * cver_datafile;
	std::vector<int> res_index;
	OctreePoint *results;
	// std::vector<float> results_weights;

	cout << "Checking All Weights" << std::endl;
	peterito->check_all_Weights();

	// octree_next->get_all_weights(results_weights);
	// james = 0;

	// if(results_weights.size()>0){
	// 	// printf("START DEBUG %lu\n",results_weights.size());
	// 	// cin >> pink;
	// printf("Before set\n"); fflush(stdout);
	printf("Before Setting iteration weights: Radial Difference: %14.5e\n", radial_difference);
	// setting_weight = 0.3 - tanh(log(radial_difference)*0.1)*0.3; if(setting_weight > 0.3) setting_weight = 0.3;
	// setting_weight = 0.2 - tanh(log(radial_difference)*0.1)*0.2; if(setting_weight > 0.2) setting_weight = 0.2;

	setting_weight = 0.01;

	if(isnormal(radial_difference))
	{
		setting_weight = 0.3 - tanh(log(radial_difference)*0.1)*0.3; if(setting_weight > 0.3) setting_weight = 0.3;
	}
	// setting_weight = 0.1;
	printf("Setting iteration weights: %.5f, %.5f, Radial Difference: %14.5e\n", setting_weight, 1.0-setting_weight, radial_difference);
	peterito->set_all_Weights(setting_weight); // Sets weights to 0.5*weight1 + 0.5*weight2
	peterito->zero_Weights(); //Sets Iter weights to 0.0
	// printf("Before set\n"); fflush(stdout);

	// }
	// for(int i=0;i<results_weights.size();++i){printf("weight %d: %.5e\n",i,results_weights[i]);}  fflush(stdout); cin >> pink;
	int cntr=0;
	int imp_index;
	float pericenter;
	float Qbind;
	float impactor_energy;

	float factor1;

	Vec3 vec1,vec2;


  	string line;
  	string getme;

  // ifstream datafile("COLL_datafile");
  // if (datafile.is_open())
  // {
  //   while ( getline(datafile,line,'\n') )
  //   {
  //     std::stringstream s(line);
  //     getline(s, getme, ' ');
  //     puts(getme);
  //   }
  //   datafile.close();
  // }
  // else cout << "Unable to open file"; 
  // return;




   old_diameter=666.666;

  t_record= 100.0*86400.0*365.25; //100 years in seconds
	int dia_index;
  // vector< pair <float,float> > SFD_vector;  // first are diameters, second is the SFD number for it



  ifstream datafile ("COLL_datafile");
  if (datafile.is_open())
  {
    while ( getline(datafile,line,'\n') )
    {
    readDatafile(line,file_to_load,diameter,pericenter);
    weight = initial_weight * pow(pericenter, pericenter_index);


	dia_index = find_diameter_index(diameter);

	// cout << "Sanity check - Diameter: "<<diameter <<", Found diameter: "<<SFD_array[dia_index][0] <<", Found Number of Particles: "<<SFD_array[dia_index][1] << "\n";

      if(old_diameter != diameter) cout << "Processing - " << file_to_load << ", Diameter: " << diameter << "\n";
      old_diameter=diameter;
      file_to_load=ReplaceAll(file_to_load,std::string(" "), std::string(""));
      cver_datafile = file_to_load.c_str();
      // printf("|%s|\n|%s|",cdata_file,cver_datafile);fflush(stdout);
      // file_to_load.replace( file_to_load.begin(), file_to_load.end(), '\\', '\\\\'); //replace all occurances of 'x' with 'y'
// file_to_load=ReplaceAll(file_to_load,std::string("/"), std::string("\\\\"));
// cout << "Processing - " << file_to_load << ", Diameter: " << diameter << "\n";
  ifstream myfile(cver_datafile, ios::in|ios::binary| ios::ate);
  current_integration_time = 9e99;

  // particle_cross_section = pow(diameter,2.0) * M_PI * 0.25 *1e-12; // Particle cross-section in meters
  particle_cross_section = M_PI * 0.25 *1e-12; // Particle cross-section in meters
  factor1 = particle_cross_section * auday2ms * t_record/volume_unit;
  particle_mass = pow(diameter,3.0)/6.0*M_PI*2000.0*1e-18 ; // Particle mass in kg
  Qbind = As * pow((diameter*0.5*1e-6),Bs) * particle_mass;
// printf("The parameters are: angle difference = %.5f (%.5f deg), velocity ratio difference = %.5f, initial weights = %.5f\n", angle_diff, acos(angle_diff)*180.0/PI,ratio_diff,initial_weight); fflush(stdout);
  if (myfile.is_open())
  {
    myfile.seekg (0, ios::beg);
		// while (!myfile.eof())
    while (true)
		{
    
    readMyLineBinary(myfile, point, velvec1,parID, time_int ); // reads one line from the input file
    // printf("File read weight - ID: %d, time: %13.5e \n",parID,time_int); fflush(stdout); 
    if(myfile.eof()) break;
    if(current_integration_time > time_int) { particle_weight.clear(); current_integration_time = time_int;}

    if(current_integration_time == time_int) // Here we start reading a new file so we start a new vector etc.
    {
    	particle_weight.push_back(weight);
    	// printf("Pushed weight - ID: %d, weight: %13.5e, init weight: %13.5e, vector lenght: %lu \n",parID,particle_weight[parID],initial_weight,particle_weight.size()); fflush(stdout); 
    	while((int)particle_weight.size()<parID)    	{printf("Particle %lu mysteriously missing at the very beginning of the datafile => next particle is %d \n", particle_weight.size(), parID);  particle_weight.push_back(weight);}
    }


    if(point.maxComponent()>=octree_limit || point.minComponent()<=-octree_limit) { continue;} // removes the points beyond the octree dimension 
	// std::cout << point.x << "\t" << point.y << "\t" << point.z << "\n";
    res_index = peterito->getIndexes(point); //
    for (int jj=0; jj<3;jj++){
    if(res_index[jj]>=400 || res_index[jj]<0) cout << "index:" << jj<< "\t is "<<res_index[jj]<< "\t" << point.x << "\t" << point.y << "\t" << point.z << "\n"; 
    }
    
    
    // octree->getPointsInsideBox({point.x-dr.x,point.y-dr.y,point.z-dr.z}, {point.x+dr.x,point.y+dr.y,point.z+dr.z}, results); //
    // cout << "First Test" << peterito->data[res_index[0]][res_index[1]][res_index[2]].getRecNum() << "\n";
    results = &peterito->data[res_index[0]][res_index[1]][res_index[2]];

    // octree_prev->getPointsInsideBox(point-dr, dr+point, results); //
    // octree->getPointsInsideBox({point.x-dr.x,point.y-dr.y,point.z-dr.z}, {point.x+dr.x,point.y+dr.y,point.z+dr.z}, results); //

    float collisional_lifetime = 0;
    float collisional_total_tau = 0;

    float pre_collision_weight = particle_weight[parID];



    if(results->getRecNum()>0){

    	// if(results.size()>1) {printf("DOES THIS EVER HAPPEN ACT 3 - Results number > 1 ? ParID: %d, Time: %14.5e, Position: %14.5e,%14.5e,%14.5e, Rec Iter: %d, Results Size: %lu\n",parID,time_int,point.x,point.y,point.z,rec_iter,results.size());fflush(stdout);}
    	// Comment: When the particle is exactly on the grid line then the results.size() = 2, it happens rarely and it does not cause measurable problems
    	// cout << results->getRecNum() << '\n'; fflush(stdout);
    	rec_iter = results->getRecNum();
    	breakable = 1;
    	for (int k=0;k<rec_iter ;++k){
    		// if(diameter > results->getDiameter(k)*30.0) {continue;} // Size condition for collisions - this will be the ballistic equation 

    		vec1 = velvec1;
    		vec2 = results->getVelocity(k);
    		getDeltaV(vec1,vec2, deltaV);
    		// deltaV = deltaV * auday2ms;
			deltaV = deltaV;

    		// CHANGE THIS
    		imp_index = results->getIndex(k);
    		impactor_energy = deltaV*deltaV*auday2ms*auday2ms*0.5*SFD_array[imp_index][2];
    		// if(diameter>600) {printf("Collision test: Qbind: %13.5e, Impactor Energy: %13.5e, Collision: %d,DeltaV: %13.5e, Target Diameter: %.5f, Impactor Diameter: %.5f\n", Qbind, impactor_energy, Qbind > impactor_energy, deltaV, diameter, SFD_array[imp_index][0]); fflush(stdout); cin >> pink;}
    		if(Qbind > impactor_energy) {continue;} // Binding Energy is larger than the kinetic energy of the impactor => no effect

    		///////////////

    		// impactor_weight = results->getWeight(k) * SFD_array[dia_index][1] / volume_unit;
    		// impactor_weight = results->getWeight(k);
    		// tau_coll = impactor_weight * particle_cross_section * deltaV * t_record;

    		impactor_weight = results->getWeight(k) * SFD_array[imp_index][1];
    		tau_coll = impactor_weight * deltaV * factor1 * pow(diameter+SFD_array[imp_index][0],2);

    		collisional_total_tau += tau_coll;
		    // 		printf("Check Something is wrong - ID: %d, time: %13.5e, weight: %13.5e, deltaV: %13.5e, impactor_weight: %13.5e, tau_coll: %13.5e \n",parID,time_int,particle_weight[parID],deltaV,current_integration_time,time_int); fflush(stdout); 
		    // 			    if(particle_weight[parID]>initial_weight){
		    // 	printf("BEFORE Something is wrong - ID: %d, time: %13.5e, weight: %13.5e, deltaV: %13.5e, impactor_weight: %13.5e, tau_coll: %13.5e \n",parID,time_int,particle_weight[parID],deltaV,impactor_weight,tau_coll); fflush(stdout); cin >> pink;
		    // }
    		// if(diameter>600 && SFD_array[imp_index][0]<600) {printf("Check Grinding - ID: %d, time: %13.5e, weight: %13.5e, deltaV: %13.5e, impactor dia: %.f, imp SFD: %13.5e, impactor_weight: %13.5e, tau_coll: %13.5e, decrease: %13.5e\n",parID,time_int,particle_weight[parID],deltaV,SFD_array[imp_index][0],SFD_array[imp_index][1],impactor_weight,tau_coll,exp(-tau_coll)); fflush(stdout); }
    		particle_weight[parID]=particle_weight[parID]*exp(-tau_coll);
    		 if(particle_weight[parID]>weight || tau_coll < 0 ){
    			printf("Something is wrong - ID: %d, time: %13.5e, weight: %13.5e, init_weight: %13.5e, deltaV: %13.5e, impactor_weight: %13.5e, tau_coll: %13.5e \n",parID,time_int,particle_weight[parID],weight,deltaV,impactor_weight,tau_coll);  fflush(stdout); cin >> debug_char;
   			}
    		}


    	// rec_iter = results->getRecNum(); // REMOVE IF WORKS
    	// breakable = 1;  // REMOVE IF WORKS
    	for (int k=0;k<rec_iter ;++k){
    		vec1 = velvec1;
    		vec2 = results->getVelocity(k);
    		cosAngleVec(velvec1,vec2, cosphi, ratio);


			
			
    		if(cosphi >= angle_diff && abs(ratio-1)<= ratio_diff && (results->getIndex(k) == dia_index))
    		{
    			results->setWeight_Iter(results->getWeight_Iter(k)+particle_weight[parID],k);
    			breakable = 0;
    			break;
    		}

    		}

    	if(breakable == 1) {
    			// This happens when the particle is exactly on the grid - the effect is minimal so we kill these
    		// for (int k=0;k<rec_iter ;++k){
    		// vec1 = velvec1;
    		// vec2 = results->getVelocity(k);
    		// cosAngleVec(velvec1,vec2, cosphi, ratio);
    		// int SFD_index = results->getIndex(k); 
    		// printf("DOES THIS EVER HAPPEN ACT 2.1? Vec1 %14.5e,%14.5e,%14.5e, Vec2: %14.5e,%14.5e,%14.5e, ratio: %14.5e, cosphi: %14.5e, rec: %d, grid_dia; %d, par_dia: %d\n",vec1.x,vec1.y,vec1.z,vec2.x,vec2.y,vec2.z,ratio, cosphi,rec_iter,SFD_index,dia_index);fflush(stdout);
    		// }
    		printf("Particle victim of a rounding error - ParID: %d, Time: %14.5e, Position: %14.5e,%14.5e,%14.5e, Rec Iter: %d\n",parID,time_int,point.x,point.y,point.z,rec_iter);fflush(stdout);
    		// cin >> debug_char;
    		// results->insertWeight(particle_weight[parID]);
    		// results->insertDiameter(diameter);
    		// results->insertVelocity(velvec1);
    		// results->addRecord();
    		}
    }
    else {
    	printf("DOES THIS EVER HAPPEN?\n");fflush(stdout);
    	perror ("No particles found in the grid - This should never happen, email the code author!");
    	abort();
    	// weight=particle_weight[parID];
    	// addOctreePoint(octree_prev, point,particle_weight[parID],dia_index,velvec1);
    	cntr+=1;
    }
 

 	// Collisions for one particle are done, now time to collect results

    float weight_loss = pre_collision_weight - particle_weight[parID]; // This is just particle weight difference
    weight_loss *= SFD_array[dia_index][1] * SFD_array[dia_index][2]/t_record/volume_unit; // We multiply by the absolute scaling * mass of the particle -> kg/s/m3; - I checked this and if t_record changes then also SFD_array[1] changes accordingly and conserves the mass loss regardless of t_record

    // CONTINUE HERE - just put weight_loss into an array (face on, edge and radial)
    if(output_switch == 1) {record_disk_mass_loss(point, weight_loss);}

    float impact_prob_with_earth;
    float impact_prob_with_probe;
    getCollProbWithEarth(point*au, velvec1*auday2ms, Vec3(1.5e11,0,0), Vec3(0,-29800.0,0), 11200e0, M_PI*6371e3*6371e3, 86400.0, a_met, e_met, inc_met, vg_met,impact_prob_with_earth,impact_prob_with_probe); // Probability per day
    // if(impact_prob_with_earth > 0.0) {fprintf(pFile666, "%d %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n", parID, time_int, a_met/au,e_met,inc_met,vg_met,impact_prob_with_earth);}

    collisional_lifetime = 1.0/collisional_total_tau;
    mass_accreted_at_Earth += impact_prob_with_earth * particle_weight[parID] * SFD_array[dia_index][1] * SFD_array[dia_index][2];


    // da = 0.06, de = 0.01, dinc = 1.8, dvg = 0.8
    // printf("Record particle: a: %13.5e, e: %13.5e,inc: %13.5e,vg: %13.5e, prob: %13.5e\n", a_met, e_met, inc_met, vg_met);
    // record_radar_profile(a_met,e_met,inc_met,vg_met, impact_prob_with_earth * particle_weight[parID]* SFD_array[dia_index][1], SFD_array[dia_index][2], "CMOR");
    // record_radar_profile(a_met,e_met,inc_met,vg_met, impact_prob_with_earth * particle_weight[parID]* SFD_array[dia_index][1], SFD_array[dia_index][2], "AMOR");
    record_radar_profile(a_met,e_met,inc_met,vg_met, impact_prob_with_earth * particle_weight[parID], dia_index, "CMOR");
    record_radar_profile(a_met,e_met,inc_met,vg_met, impact_prob_with_earth * particle_weight[parID], dia_index, "AMOR");
	if(output_switch == 1) {record_disk_profiles(point, particle_weight[parID]*SFD_array[dia_index][1], diameter, deltaV*auday2ms);}

    radialProfile(r_prof_iter, point, particle_weight[parID]*SFD_array[dia_index][1] ); // adds the read line to the radial profile "r_prof"
    radialProfileCross(r_prof_brightness, point, particle_weight[parID]*SFD_array[dia_index][1],pow(diameter,2.0)*0.25*M_PI*1e-12 ); // adds the read line to the radial profile "r_prof"
    eclipticLatitudeProfile(ecl_profile, ecl_profile_mass, point, particle_weight[parID]*SFD_array[dia_index][1], diameter);
    SFD_1AU_Profile(SFD_profile,point,impact_prob_with_probe*particle_weight[parID]*SFD_array[dia_index][1], diameter);
    SFD_1AU_Profile_Kessler(SFD_profile_Kessler,point,impact_prob_with_earth*particle_weight[parID]*SFD_array[dia_index][1], diameter);

    test_zody_area += particle_weight[parID]*SFD_array[dia_index][1]*pow(diameter,2.0)*0.25*M_PI*1e-12;
    test_zody_mass += particle_weight[parID]*SFD_array[dia_index][1]*pow(diameter,3.0)/6.0*M_PI*1e-18*particle_density; 
    // results.clear();

	
	}

		

    myfile.close();
  }
  else cout << "Unable to open binary file\n"; 

      }
    datafile.close();

  }
  else cout << "Unable to open file\n" ;

  // peterito->check_all_Weights();

  cout << "Testing Zone: Area" << test_zody_area << ", Mass: " << test_zody_mass << std::endl;
  // checkRadialProfile(r_prof_iter, 1);
  return;
}
/// END OF ITERATION


void calculate_PQ_vectors(float inc, float capom, float omega, Vec3 &P_vec, Vec3 &Q_vec)
{
	// inc, capom, omega are in radians!
  P_vec[0]=cos(capom)*cos(omega)-cos(inc)*sin(omega)*sin(capom);
  P_vec[1]=sin(capom)*cos(omega)+cos(inc)*sin(omega)*cos(capom);
  P_vec[2]=sin(inc)*sin(omega);

  Q_vec[0]=-cos(capom)*sin(omega)-cos(inc)*cos(omega)*sin(capom);
  Q_vec[1]=-sin(capom)*sin(omega)+cos(inc)*cos(omega)*cos(capom);
  Q_vec[2]=sin(inc)*cos(omega);
}

void calculate_radial_vector(Vec3 P_vec, Vec3 Q_vec, float a, float e, float true_anom, Vec3 &R_vec)
{
 float r_helio;
 r_helio=a*(1.0-e*e)/(1.0+e*cos(true_anom));
 R_vec = r_helio * (cos(true_anom)*P_vec+sin(true_anom)*Q_vec);
}

void calculate_velocity_vector(Vec3 P_vec, Vec3 Q_vec, float a, float e, float true_anom, Vec3 &V_vec)
{
	// to get MKS units we use au and GM_sun from the constant pameters set in the header
 float n; // mean motion
 // a=a*1.5e11
 a = a * au;
 // GM=1.32712440018d20
 n = sqrt(GM_sun*pow((a),-3.0)); // TO METERS per SECOND
 // cout << "Velocity Vector Calculation: " << a << "\t" << au << "\t" << GM_sun << std::endl;
 V_vec=n*a/sqrt(1.0-e*e)*(-sin(true_anom)*P_vec+(e+cos(true_anom))*Q_vec);
 // cout << "Velocity Vector Calculation - 2: " << V_vec.x << "\t" << V_vec.y << "\t" << V_vec.z  << std::endl;
}

void convert_orbel_to_xyz(float a, float e, float inc, float capom, float omega, float true_anom, 
						  Vec3 &R_vec, Vec3 &V_vec){
						  // float &x, float &y, float &z, float &vx, float &vy, float &vz) {


Vec3 P_vec;
Vec3 Q_vec;
// Vec3 R_vec;
// Vec3 V_vec;
calculate_PQ_vectors(inc, capom, omega, P_vec, Q_vec);;
calculate_radial_vector(P_vec, Q_vec,  a,  e,  true_anom,  R_vec);
calculate_velocity_vector(P_vec, Q_vec,  a,  e,  true_anom,  V_vec);

// x = R_vec.x;
// y = R_vec.y;
// z = R_vec.z;

// vx = V_vec.x;
// vy = V_vec.y;
// vz = V_vec.z;
}


float solve_kepler_equation_for_TA(float M, float ecc)
/// Solves Kepler Equation for given mean anomaly M, eccentricty ecc, and returns float true anomaly
{
    int i;
    int iterations_to_converge;
    int iteration_maximum = 20;
    double sgn, E, f, error, es, ec, df, ddf, dddf, d1, d2, d3;
    double max_error = 1e-10;
    double cosf;

    
    M = fmod(M, 2*M_PI);
    if (M >  M_PI) M = M - 2.0*M_PI;
    if (M < -M_PI) M = M + 2.0*M_PI;

    //Construct the initial solultion.
    sgn = 1.0;
    if (M < 0.0) sgn = -1.0;
    E = M + sgn*(0.85)*ecc;

    //Solve kepler's equation iteratively to improve the solution E.
    error = 1.0;
    for(i = 0; i < iteration_maximum; i++){
      es = ecc*sin(E);
      ec = ecc*cos(E);
      f = E - es - M;
      error = fabs(f);
      if (error < max_error) break;
      df = 1.0 - ec;
      ddf = es;
      dddf = ec;
      d1 = -f/df;
      d2 = -f/(df + d1*ddf/2.0);
      d3 = -f/(df + d2*ddf/2.0 + d2*d2*dddf/6.0);
      E = E + d3;
    }

     //Warn if solution did not converge.
     if (error > max_error)
       std::cout << "***Warning*** Orbit::keplers_eqn() failed to converge***\n";

    iterations_to_converge = i;
    cosf = (cos(E) - ecc) / (1.0 - ecc*cos(E));
    return acos(cosf);
    }

void calculate_collisional_lifetime_map(){


std::vector<double> particle_weight;
	
	float setting_weight;
	int rec_iter,breakable;
	
	
	// Meteoroid orbital elements
	float a_met, e_met, inc_met, vg_met;

	// float weight;
	float cosphi;
	float ratio;
	float diameter;
	float old_diameter;
	float weight;
	int parID;
	float time_int;
	float current_integration_time = 9e99;// Make this number high so it always higher than the first line 
	float deltaV;
	float impactor_weight;
	float particle_cross_section;
	float particle_mass;
	float t_record;
	float tau_coll;

	double test_zody_mass = 0.0, test_zody_area = 0.0;



	

	Vec3 velvec1,point;

  char debug_char;
  
  string file_to_load;

  const char * cver_datafile;
	std::vector<int> res_index;
	OctreePoint *results;
	// std::vector<float> results_weights;

	cout << "Checking All Weights" << std::endl;
	peterito->check_all_Weights();

	
	int cntr=0;
	int imp_index;
	int i1_max;
	float pericenter;
	float Qbind;
	float impactor_energy;

	float factor1;

	float mean_anomaly_met;
	float true_anomaly_met;
	Vec3 vec1,vec2;


  	string line;
  	string getme;

 

  i1_max = 30;
  // i1_max = 1;
  int omega_max = 30;
  // int omega_max = 1;

  for( int a_iter=0; a_iter<26; a_iter++) {

  diameter=pow(10.0,double( (a_iter-2)/6.0  ));
  t_record= 100.0*86400.0*365.25; //100 years in seconds

  particle_cross_section = M_PI * 0.25 *1e-12; // Particle cross-section in meters
  factor1 = particle_cross_section * auday2ms * t_record/volume_unit;
  particle_mass = pow(diameter,3.0)/6.0*M_PI*2000.0*1e-18 ; // Particle mass in kg
  Qbind = As * pow((diameter*0.5*1e-6),Bs) * particle_mass;


  FILE * pFile_temp;
  pFile_temp = fopen(  (std::string("Collisional_Map_D_met_")+std::to_string(diameter)).c_str(),"w+");
  fprintf(pFile_temp,"#Average_Collisional_Lifetime  Ecc Inc\n");

  a_met = 0.3871;
  
  int iteration_resolution = 1;
  // int iteration_resolution = 10;

  for (int ecc_int=0; ecc_int<iteration_resolution; ecc_int++) {
  	 // e_met = float(ecc_int)/float(iteration_resolution)*0.1;
  	 e_met = 0.1;
  	for (int inc_int=0; inc_int<iteration_resolution; inc_int++) {
  		// inc_met = float(inc_int)/float(iteration_resolution)*7.0*deg2rad;
  		inc_met = 7.0*deg2rad;

  double average_collisional_lifetime = 0.0;

  for (int i1=0; i1<i1_max; i1++){
  
 
  
  mean_anomaly_met = float(i1)/float(i1_max)*2.0*M_PI;

  true_anomaly_met = solve_kepler_equation_for_TA(mean_anomaly_met, e_met);


  for(int capom_int=0; capom_int<omega_max; capom_int++) {
  float capom_met = float(capom_int)/float(omega_max)*2.0*M_PI;
  for(int omega_int=0; omega_int<omega_max; omega_int++) {

  float omega_met = float(omega_int)/float(omega_max)*2.0*M_PI;
  // i1_max = 100;
  


  convert_orbel_to_xyz(a_met,e_met,inc_met,capom_met,omega_met,true_anomaly_met,point,velvec1);

    if(point.maxComponent()>=octree_limit || point.minComponent()<=-octree_limit) { continue;} // removes the points beyond the octree dimension 
	// std::cout << point.x << "\t" << point.y << "\t" << point.z << "\n";
    res_index = peterito->getIndexes(point); //
    for (int jj=0; jj<3;jj++){
    if(res_index[jj]>=400 || res_index[jj]<0) cout << "index:" << jj<< "\t is "<<res_index[jj]<< "\t" << point.x << "\t" << point.y << "\t" << point.z << "\n"; 
    }
    
    
     results = &peterito->data[res_index[0]][res_index[1]][res_index[2]];


    float collisional_lifetime = 0;
    float collisional_total_tau = 0;

    if(results->getRecNum()>0){

    	rec_iter = results->getRecNum();
    	for (int k=0;k<rec_iter ;++k){

    		vec1 = velvec1*(1.0/auday2ms);
    		vec2 = results->getVelocity(k);
    		// cout << velvec1.x << "\t" << velvec1.y << "\t" << velvec1.z  << "\t" << auday2ms << std::endl;
    		// cout << velvec1.x*(1.0/auday2ms) << "\t" << velvec1.y*(1.0/auday2ms) << "\t" << velvec1.z*(1.0/auday2ms)  << "\t - WTF ---  " << float(auday2ms)/float(auday2ms) << std::endl;
    		// cout << vec1.x << "\t" << vec1.y << "\t" << vec1.z  << std::endl;
    		// cout << vec2.x << "\t" << vec2.y << "\t" << vec2.z  << std::endl;
    		// cin >> debug_char;
    		getDeltaV(vec1,vec2, deltaV);
    		// deltaV = deltaV * auday2ms;
			deltaV = deltaV;

    		// CHANGE THIS
    		imp_index = results->getIndex(k);
    		impactor_energy = deltaV*deltaV*auday2ms*auday2ms*0.5*SFD_array[imp_index][2];
    		// if(diameter>600) {printf("Collision test: Qbind: %13.5e, Impactor Energy: %13.5e, Collision: %d,DeltaV: %13.5e, Target Diameter: %.5f, Impactor Diameter: %.5f\n", Qbind, impactor_energy, Qbind > impactor_energy, deltaV, diameter, SFD_array[imp_index][0]); fflush(stdout); cin >> pink;}
    		if(Qbind > impactor_energy) {continue;} // Binding Energy is larger than the kinetic energy of the impactor => no effect

    		impactor_weight = results->getWeight(k) * SFD_array[imp_index][1];
    		tau_coll = impactor_weight * deltaV * factor1 * pow(diameter+SFD_array[imp_index][0],2);

    		collisional_total_tau += tau_coll/t_record;
    		// printf("Collision test: Qbind: %13.5e, Impactor Energy: %13.5e, Collision: %d,DeltaV: %13.5e, Target Diameter: %.5f, Impactor Diameter: %.5f\n", Qbind, impactor_energy, Qbind > impactor_energy, deltaV, diameter, SFD_array[imp_index][0]); fflush(stdout);
    }

    
    }
    else {
    	// printf("DOES THIS EVER HAPPEN?\n");fflush(stdout);
    	// weight=particle_weight[parID];
    	// addOctreePoint(octree_prev, point,particle_weight[parID],dia_index,velvec1);
    	cntr+=1;
    }
 
    average_collisional_lifetime += collisional_total_tau*86400.0*365.25;
}	

}
}   
// cout << "Average Collisional Lifetime: " << average_collisional_lifetime/double(i1_max*i1_max*i1_max) <<" Ecc: " << e_met << " Inc: " << inc_met <<  std::endl;
cout << "Average Collisional Lifetime: " << 1.0/(average_collisional_lifetime/double(i1_max*i1_max*i1_max)) << " A_met: " << a_met << " Ecc: " << e_met << " Inc: " << inc_met << " Diameter: " << diameter << std::endl;
fprintf(pFile_temp,"%13.5e %13.5e %13.5e  %13.5e\n", 1.0/(average_collisional_lifetime/double(i1_max*i1_max*i1_max)) , e_met, inc_met, diameter);
}
}


fclose(pFile_temp);
}
  return;



}


void calculate_collisional_lifetime_plot(){


	std::vector<double> particle_weight;

	float setting_weight;
	int rec_iter,breakable;


// Meteoroid orbital elements
	float a_met, e_met, inc_met, vg_met;

// float weight;
	float cosphi;
	float ratio;
	float diameter;
	float old_diameter;
	float weight;
	int parID;
	float time_int;
	float current_integration_time = 9e99;// Make this number high so it always higher than the first line 
	float deltaV;
	float impactor_weight;
	float particle_cross_section;
	float particle_mass;
	float t_record;
	float tau_coll;
	double test_zody_mass = 0.0, test_zody_area = 0.0;

	Vec3 velvec1,point;
	char debug_char;
	string file_to_load;
	const char * cver_datafile;
	std::vector<int> res_index;
	OctreePoint *results;
	int cntr=0;
	int imp_index;
	int i1_max;
	float pericenter;
	float Qbind;
	float impactor_energy;
	float factor1;
	float mean_anomaly_met;
	float true_anomaly_met;
	Vec3 vec1,vec2;
	string line;
	string getme;


	diameter=100.0;
	t_record= 100.0*86400.0*365.25; //100 years in seconds



	i1_max = 30;

	int omega_max = 30;

	e_met = 0.0;
	inc_met = 0.0;

for( int a_iter=0; a_iter<1; a_iter++) {

	FILE * pFile_temp;
	pFile_temp = fopen(  (std::string("Collisional_Map_diameter_")+std::to_string(a_iter)).c_str(),"w+");


	int iteration_resolution = 100;

	for (int ecc_int=0; ecc_int<iteration_resolution; ecc_int++) {
		diameter = pow(10.0,double(ecc_int*4.0)/double(iteration_resolution));

		particle_cross_section = pow(diameter,2.0) * M_PI * 0.25 *1e-12; // Particle cross-section in meters
		factor1 = particle_cross_section * auday2ms * t_record/volume_unit;
		particle_mass = pow(diameter,3.0)/6.0*M_PI*2000.0*1e-18 ; // Particle mass in kg
		Qbind = As * pow((diameter*0.5*1e-6),Bs) * particle_mass;


		for (int inc_int=1; inc_int<iteration_resolution; inc_int++) {
			a_met = float(inc_int)/float(iteration_resolution)*5.0;

			double average_collisional_lifetime = 0.0;

			for (int i1=0; i1<i1_max; i1++){



				mean_anomaly_met = float(i1)/float(i1_max)*2.0*M_PI;

				true_anomaly_met = solve_kepler_equation_for_TA(mean_anomaly_met, e_met);


				for(int capom_int=0; capom_int<omega_max; capom_int++) {
					float capom_met = float(capom_int)/float(omega_max)*2.0*M_PI;
					for(int omega_int=0; omega_int<omega_max; omega_int++) {

						float omega_met = float(omega_int)/float(omega_max)*2.0*M_PI;
// i1_max = 100;



						convert_orbel_to_xyz(a_met,e_met,inc_met,capom_met,omega_met,true_anomaly_met,point,velvec1);

if(point.maxComponent()>=octree_limit || point.minComponent()<=-octree_limit) { continue;} // removes the points beyond the octree dimension 
// std::cout << point.x << "\t" << point.y << "\t" << point.z << "\n";
res_index = peterito->getIndexes(point); //
for (int jj=0; jj<3;jj++){
	if(res_index[jj]>=400 || res_index[jj]<0) cout << "index:" << jj<< "\t is "<<res_index[jj]<< "\t" << point.x << "\t" << point.y << "\t" << point.z << "\n"; 
}


results = &peterito->data[res_index[0]][res_index[1]][res_index[2]];


float collisional_lifetime = 0;
float collisional_total_tau = 0;

if(results->getRecNum()>0){

	rec_iter = results->getRecNum();
	for (int k=0;k<rec_iter ;++k){

		vec1 = velvec1*(1.0/auday2ms);
		vec2 = results->getVelocity(k);
// cout << velvec1.x << "\t" << velvec1.y << "\t" << velvec1.z  << "\t" << auday2ms << std::endl;
// cout << velvec1.x*(1.0/auday2ms) << "\t" << velvec1.y*(1.0/auday2ms) << "\t" << velvec1.z*(1.0/auday2ms)  << "\t - WTF ---  " << float(auday2ms)/float(auday2ms) << std::endl;
// cout << vec1.x << "\t" << vec1.y << "\t" << vec1.z  << std::endl;
// cout << vec2.x << "\t" << vec2.y << "\t" << vec2.z  << std::endl;
// cin >> debug_char;
		getDeltaV(vec1,vec2, deltaV);
// deltaV = deltaV * auday2ms;
		deltaV = deltaV;

// CHANGE THIS
		imp_index = results->getIndex(k);
		impactor_energy = deltaV*deltaV*auday2ms*auday2ms*0.5*SFD_array[imp_index][2];
// if(diameter>600) {printf("Collision test: Qbind: %13.5e, Impactor Energy: %13.5e, Collision: %d,DeltaV: %13.5e, Target Diameter: %.5f, Impactor Diameter: %.5f\n", Qbind, impactor_energy, Qbind > impactor_energy, deltaV, diameter, SFD_array[imp_index][0]); fflush(stdout); cin >> pink;}
if(Qbind > impactor_energy) {continue;} // Binding Energy is larger than the kinetic energy of the impactor => no effect

impactor_weight = results->getWeight(k) * SFD_array[imp_index][1];
tau_coll = impactor_weight * deltaV * factor1;

collisional_total_tau += tau_coll/t_record;
// printf("Collision test: Qbind: %13.5e, Impactor Energy: %13.5e, Collision: %d,DeltaV: %13.5e, Target Diameter: %.5f, Impactor Diameter: %.5f\n", Qbind, impactor_energy, Qbind > impactor_energy, deltaV, diameter, SFD_array[imp_index][0]); fflush(stdout);
}


}
else {
// printf("DOES THIS EVER HAPPEN?\n");fflush(stdout);
// weight=particle_weight[parID];
// addOctreePoint(octree_prev, point,particle_weight[parID],dia_index,velvec1);
	cntr+=1;
}

average_collisional_lifetime += collisional_total_tau*86400.0*365.25;
}	

}
}   
// cout << "Average Collisional Lifetime: " << average_collisional_lifetime/double(i1_max*i1_max*i1_max) <<" Ecc: " << e_met << " Inc: " << inc_met <<  std::endl;
cout << "Average Collisional Lifetime: " << 1.0/(average_collisional_lifetime/double(i1_max*i1_max*i1_max)) << " A_met: " << a_met << " Diameter: " << diameter <<  std::endl;
fprintf(pFile_temp,"Average Collisional Lifetime: %13.5e  Ecc: %13.5e  Inc: %13.5e\n", 1.0/(average_collisional_lifetime/double(i1_max*i1_max*i1_max)) , a_met, diameter);
}
}


fclose(pFile_temp);
}
return;



}


int main(int argc, char **argv) {

    std::feclearexcept(FE_OVERFLOW);
    std::feclearexcept(FE_UNDERFLOW);
	// TESTING AREA

    int load_save_switch = 0;
    int reweight_switch = 1;

	// END OF TESTING AREA


	std::size_t found;

	double As_modifier = 1;
	alpha = 4.1;
	beta = 3.7;
	Dmid = 60.0;
	As = 1e6*1e-7*1e3; // J kg-1 - Using Krivov et al. (2006)

	Bs = -0.24;
	pericenter_index = -1.3;

    std::ifstream file("input.in");
    if(file.is_open()) {
    std::string str; 
    int input_cntr = 0;
    while (std::getline(file, str))
    {
    	found = str.find("#");
    	// cout << found << "\n";
    	if(found != string::npos)str.erase(str.begin()+found, str.end());
    	if(input_cntr == 0) alpha = stof(str);
    	if(input_cntr == 1) beta = stof(str);
    	if(input_cntr == 2) Dmid = stof(str);
    	if(input_cntr == 3) As_modifier = stof(str);
    	if(input_cntr == 4) Bs = stof(str);
    	if(input_cntr == 5) intended_cloud_area = stof(str);
    	if(input_cntr == 6) pericenter_index = stof(str);
    	if(input_cntr == 7) work_directory = trim(str);
    	if(input_cntr == 8) output_switch = stoi(str);
    	if(input_cntr == 9) load_save_switch = stoi(str);
    	if(input_cntr == 10) reweight_switch = stoi(str);
    	
        // cout << str << input_cntr << "\n";
        input_cntr ++;
    }

    } else {cout << "File not found! \t " << "input.in" << "\n";}

    initial_weight = intended_cloud_area;
	load_radar_constraints();

	As = As * As_modifier;
	
	// pericenter_index = 0.0;


	Dmax = 3000.0;
	N0 = 1.0;
	// float As = 2e5; // erg g-1

	particle_density = 2000; // kg m^-3

	cout << "Fitting parameters: alpha = " << alpha << ", beta = " << beta << ", Dmid = " << Dmid << ",As = " << As << ", Bs = " << Bs << ", Intended Cloud Area: " << intended_cloud_area << "\n";
	cout << "Fitting parameters: pericenter_index = " << pericenter_index << ", Output_Switch = " << output_switch << "\n";
	// Parameters

	peterito = new Grid3D(Vec3(10.0,10.0,10.0),Vec3(0.05,0.05,0.05));
	double finish_critetion = 1e-2;




    // return 0;

	// Variables
	double radial_difference = 1.0;
	double T;
	int iter_count = 0;
	float tree_difference = 1.0;
	// std::vector<OctreePoint*> results;
	
	std::string line;

	pFile = fopen (cradial_file,"w+");
	pFile1 = fopen (cecl_file,"w+");
	pFile2 = fopen ("Outputs_SFD_profile.out","w+");
	pFile3 = fopen ("Outputs_Radar_Profiles_CMOR.out","w+");
	pFile4 = fopen ("Outputs_Radar_Profiles_AMOR.out","w+");
	// pFile666 = fopen ("TEST_probabilities","w+");
   	
	double start = stopwatch();
	
	pre_init();
	calculate_SFD(alpha,beta,Dmid,Dmax,N0);




	for(int i=0;i<(int)SFD_array.size();i++){
		printf("%d, Diameter: %.f, SFD Number: %13.5e, Intergration particle count: %13.5e\n", i, SFD_array[i][0],SFD_array[i][1],Integration_Setup_array[i][1]);
	}

	
	Initialize_Outputs();
    Initiate_Output_Files();


    if(load_save_switch != 1) {	init(); } // Initialize the seed model 

	if(load_save_switch == 1) {
	initial_weight = peterito->load_grid3D();
	cout << "GRID SUCCESSFULLY LOADED!!!!!" << std::endl;

	// calculate_collisional_lifetime_map();
	// // calculate_collisional_lifetime_plot();
	// return 1;


	}

	T = stopwatch() - start; printf("Time-elapsed Init: %.5f - total number of particles loaded: %d\n",T,total_number_of_particles_integration); start = stopwatch();fflush(stdout); // Timer

	float zody_area ;
	float zody_mass ;
	float Earth_mass;

	// While loop for iterations - finishes when the finish_criterion is reached
	while(radial_difference>finish_critetion) {
	// while(tree_difference>finish_critetion) {

		mass_accreted_at_Earth = 0.0;
		tree_difference = 0.0;
		std::fill(ecl_profile.begin(),ecl_profile.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		std::fill(ecl_profile_mass.begin(),ecl_profile_mass.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		std::fill(SFD_profile.begin(),SFD_profile.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		std::fill(SFD_profile_Kessler.begin(),SFD_profile_Kessler.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		std::fill(r_prof_brightness.begin(),r_prof_brightness.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there

		// cout << "T1 " << My_Measurements.CMOR_AH_a_x.size() << "\t" << My_Measurements.CMOR_AH_a_y.size() << "\t" << My_Measurements.CMOR_AH_ecc_x.size() << "\t" << My_Measurements.CMOR_AH_ecc_y.size() << "\n";

		initialize_radar_measurements();

		// cout << My_Constraints.CMOR_AH_a_x.size() << "\t" << My_Constraints.CMOR_AH_a_y.size() << "\t" << My_Constraints.CMOR_AH_ecc_x.size() << "\t" << My_Constraints.CMOR_AH_ecc_y.size() << "\n";
		// cout << My_Measurements.CMOR_AH_a_x.size() << "\t" << My_Measurements.CMOR_AH_a_y.size() << "\t" << My_Measurements.CMOR_AH_ecc_x.size() << "\t" << My_Measurements.CMOR_AH_ecc_y.size() << "\n";

		cout << "Before iteration - area profile sum: " << compute_sum(ecl_profile) << ", mass profile sum: " << compute_sum(ecl_profile_mass) << std::endl;
		iteration(radial_difference);
		cout << "After iteration - area profile sum: " << compute_sum(ecl_profile) << ", mass profile sum: " << compute_sum(ecl_profile_mass) << std::endl;
		normalize_radar_measurements();
		printf("Mass accreted at Earth: %13.5e \n", mass_accreted_at_Earth);

		

	// 	for (int i=0;i<My_Measurements.CMOR_AH_a_y.size(); ++i){
	// 	printf("RADAR CHECK --- %.5f %13.5e %.5f %13.5e %.5f %13.5e %.5f %13.5e %d \n", 
	// 		My_Constraints.CMOR_AH_a_x[i], My_Measurements.CMOR_AH_a_y[i],
	// 		My_Constraints.CMOR_AH_ecc_x[i], My_Measurements.CMOR_AH_ecc_y[i],
	// 		My_Constraints.CMOR_AH_inc_x[i], My_Measurements.CMOR_AH_inc_y[i],
	// 		My_Constraints.CMOR_AH_vg_x[i], My_Measurements.CMOR_AH_vg_y[i],
	// 		iter_count );
	// }


		// pfileTEST = fopen ("TEST_probability","w+");
		// for(int TEST_I=0; TEST_I < TEST_probabilities.size(); TEST_I++){
		// fprintf(pfileTEST,"%15.7e\n",TEST_probabilities[TEST_I]);
		// 	}
		// return 1;
		iter_count += 1;
		// radial_difference = getRadialDifference(r_prof_init,r_prof_iter); 
		radial_difference = peterito->get_the_iteration_difference();
		cout << "Radial Difference New " << radial_difference << std::endl;
		radial_difference = getRadialDifference(r_prof_init,r_prof_iter); 
		cout << "Radial Difference Oold " << radial_difference << std::endl;

		printf("\n ---------------------- \n");

		zody_area = calculate_cloud_area();
		zody_mass = calculate_cloud_mass();

		T = stopwatch() - start; printf("Time-elapsed Iter %d: %.5f, Radial difference: %.5e, Zody area: %.5e, Zody mass %.5e \n",iter_count,T,radial_difference, zody_area, zody_mass); start = stopwatch();printf("---------------------- \n"); fflush(stdout);
		checkRadialProfile(pFile, r_prof_iter, iter_count);fprintf(pFile,"\n");
		checkEclipticProfile(pFile1, ecl_profile, ecl_profile_mass, iter_count);fprintf(pFile1,"\n");
		checkSFDProfile(pFile2, SFD_profile, iter_count);fprintf(pFile2,"\n");
		write_radar_profile(pFile3, iter_count,"CMOR");fprintf(pFile3,"\n");
		write_radar_profile(pFile4, iter_count,"AMOR");fprintf(pFile4,"\n");
		if(output_switch == 1) {Print_Outputs_Big(iter_count);}
		Print_Outputs_Small(iter_count);
		Zero_Outputs();
		std::copy(r_prof_iter.begin(), r_prof_iter.end(), r_prof_init.begin()); std::fill(r_prof_iter.begin(), r_prof_iter.end(), 0); // Copy the r_prof_iter to r_prof_init and put zeros to r_prof_iters
	}
	fclose (pFile);
	fclose (pFile1);
	fclose (pFile2);
	fclose (pFile3);
	fclose (pFile4);



	pFile4 = fopen ("Outputs_Mass_Accreted","w+");
	std::vector<double> mass_vector;
	calculate_Earth_Mass_Flux(mass_vector,iter_count);
	fclose(pFile4);


	// float zody_area = compute_sum(ecl_profile);
	// float zody_mass = compute_sum(ecl_profile_mass);

	zody_area = calculate_cloud_area();
	zody_mass = calculate_cloud_mass();
	Earth_mass = compute_average(mass_vector);

	float ratio1 = zody_area/intended_cloud_area;
	float ratio2 = zody_mass/4e16;
	float ratio3 = Earth_mass/4e4;
	float reweight_factor = 9e9;


	float estimate1;
	float estimate2;
	float estimate3;
	std::vector<float> estimate_vec;


	// // A more sophisticated way of reweighting - for now we reweight just to match the Zody Area
	// estimate1 = abs(log( ratio1/ratio1)) +abs(log( ratio2/ratio1 )) +abs(log( ratio3/ratio1));
	// estimate2 = abs(log( ratio1/ratio2)) +abs(log( ratio2/ratio2 )) +abs(log( ratio3/ratio2));
	// estimate3 = abs(log( ratio1/ratio3)) +abs(log( ratio2/ratio2 )) +abs(log( ratio3/ratio2));
	// estimate_vec = {estimate1,estimate2,estimate3};
	// if(estimate1 == *min_element(estimate_vec.begin(),estimate_vec.end())) reweight_factor = 1.0/ratio1;
	// if(estimate2 == *min_element(estimate_vec.begin(),estimate_vec.end())) reweight_factor = 1.0/ratio2;
	// if(estimate3 == *min_element(estimate_vec.begin(),estimate_vec.end())) reweight_factor = 1.0/ratio3;




	// if(ratio1 < 1.0) {reweight_factor = abs(ratio1-1.0)}
	reweight_factor = 1.0/ratio1;
	// reweight_factor =  exp( ( log(ratio1)+log(ratio2)+log(ratio3) ) / (-3.0)) ; // OLD REWEIGHT
	// float quality_factor = std::max(zody_area/2e17,pow(zody_area/2e17,-1))*std::max(zody_mass/4e16,pow(zody_mass/4e16,-1))*std::max(Earth_mass/4e4,pow(Earth_mass/4e4,-1));
	// float quality_factor = pow(log( ratio1),2) + pow(log( ratio2),2) + pow(log( ratio3),2); // OLD QUALITY
	float quality_factor = abs(log( ratio1)) +abs(log( ratio2 )) +abs(log( ratio3));

	cout << "Area: " << ratio1 << ", Zody Mass: " << ratio2 << ", Earth Mass: " << ratio3 << ", Quality factor: "<< quality_factor << ", Reweight Factor: "<< reweight_factor << "\n";



	if(reweight_switch == 1) 
	{ 
	while(pow(log(reweight_factor),2.0)>0.009085) { // reweighing factor is less than 10%

	// peterito->zero_ALL_Weights();	
	// init(); // Initialize the seed model
	initial_weight = reweight_factor * initial_weight;

	radial_difference = 1.0;
	iter_count += 100;

	printf("\n ---------------------- \n");
	cout << "Recalculating the initial weight: " << initial_weight << ", reweighting factor: " << reweight_factor<<  "\n";
	peterito->reweight_ALL_Weights(reweight_factor);
	printf("\n ---------------------- \n");

	pFile = fopen (cradial_file,"a+");
	pFile1 = fopen (cecl_file,"a+");
	pFile2 = fopen ("Outputs_SFD_profile.out","a+");
	pFile3 = fopen ("Outputs_Radar_Profiles_CMOR.out","a+");
	pFile4 = fopen ("Outputs_Radar_Profiles_AMOR.out","a+");

	while(radial_difference>finish_critetion) {
	// while(tree_difference>finish_critetion) {
		tree_difference = 0.0;
		mass_accreted_at_Earth = 0.0;
		std::fill(ecl_profile.begin(),ecl_profile.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		std::fill(ecl_profile_mass.begin(),ecl_profile_mass.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		std::fill(SFD_profile.begin(),SFD_profile.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		std::fill(SFD_profile_Kessler.begin(),SFD_profile_Kessler.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		std::fill(r_prof_brightness.begin(),r_prof_brightness.end(),0); // Fills the ecliptic latitude profile with zeros before we start putting numbers there
		initialize_radar_measurements();

		iteration(radial_difference);
		normalize_radar_measurements();

		printf("Mass accreted at Earth: %13.5e \n", mass_accreted_at_Earth);
		iter_count += 1;
		// radial_difference = getRadialDifference(r_prof_init,r_prof_iter); 
		radial_difference = peterito->get_the_iteration_difference();
		cout << "Radial Difference New " << radial_difference << std::endl;
		radial_difference = getRadialDifference(r_prof_init,r_prof_iter); 
		cout << "Radial Difference Oold " << radial_difference << std::endl;
		printf("\n ---------------------- \n");
		zody_area = calculate_cloud_area();
		zody_mass = calculate_cloud_mass();

		T = stopwatch() - start; printf("Time-elapsed Iter %d: %.5f, Radial difference: %.5e, Zody area: %.5e, Zody mass %.5e \n",iter_count,T,radial_difference, zody_area, zody_mass); start = stopwatch();printf("---------------------- \n"); fflush(stdout);
		checkRadialProfile(pFile, r_prof_iter, iter_count);fprintf(pFile,"\n");
		checkEclipticProfile(pFile1, ecl_profile, ecl_profile_mass, iter_count);fprintf(pFile1,"\n");
		checkSFDProfile(pFile2, SFD_profile, iter_count);fprintf(pFile2,"\n");
		write_radar_profile(pFile3, iter_count,"CMOR");fprintf(pFile3,"\n");
		write_radar_profile(pFile4, iter_count,"AMOR");fprintf(pFile4,"\n");
		std::copy(r_prof_iter.begin(), r_prof_iter.end(), r_prof_init.begin()); std::fill(r_prof_iter.begin(), r_prof_iter.end(), 0); // Copy the r_prof_iter to r_prof_init and put zeros to r_prof_iters
		if(output_switch == 1) {Print_Outputs_Big(iter_count);}
		Print_Outputs_Small(iter_count);
		Zero_Outputs();
	}
	fclose (pFile);
	fclose (pFile1);
	fclose (pFile2);
	fclose (pFile3);
	fclose (pFile4);

	pFile4 = fopen ("Outputs_Mass_Accreted","a+");
	mass_vector.clear(); 
	calculate_Earth_Mass_Flux(mass_vector,iter_count);
	fclose(pFile4);


	zody_area = calculate_cloud_area();
	zody_mass = calculate_cloud_mass();
	Earth_mass = compute_average(mass_vector);

	ratio1 = zody_area/intended_cloud_area;
	ratio2 = zody_mass/4e16;
	ratio3 = Earth_mass/4e4;
	// float quality_factor = std::max(zody_area/2e17,pow(zody_area/2e17,-1))*std::max(zody_mass/4e16,pow(zody_mass/4e16,-1))*std::max(Earth_mass/4e4,pow(Earth_mass/4e4,-1));
	// quality_factor = pow(log( ratio1),2) + pow(log( ratio2),2) + pow(log( ratio3),2);
	quality_factor = abs(log( ratio1)) +abs(log( ratio2 )) +abs(log( ratio3));

	cout << "Area: " << ratio1 << ", Zody Mass: " << ratio2 << ", Earth Mass: " << ratio3 << ", Quality factor: "<< quality_factor << "\n";
	// reweight_factor =  exp( ( log(ratio1)+log(ratio2)+log(ratio3) ) / (-3.0)) ;

	// estimate1 = abs(log( ratio1/ratio1)) +abs(log( ratio2/ratio1 )) +abs(log( ratio3/ratio1));
	// estimate2 = abs(log( ratio1/ratio2)) +abs(log( ratio2/ratio2 )) +abs(log( ratio3/ratio2));
	// estimate3 = abs(log( ratio1/ratio3)) +abs(log( ratio2/ratio2 )) +abs(log( ratio3/ratio2));
	// estimate_vec = {estimate1,estimate2,estimate3};
	// if(estimate1 == *min_element(estimate_vec.begin(),estimate_vec.end())) reweight_factor = 1.0/ratio1;
	// if(estimate2 == *min_element(estimate_vec.begin(),estimate_vec.end())) reweight_factor = 1.0/ratio2;
	// if(estimate3 == *min_element(estimate_vec.begin(),estimate_vec.end())) reweight_factor = 1.0/ratio3;
	reweight_factor = 1.0/ratio1;

}
}
	pFile = fopen("Fitting_Parameters.results", "a+");
	fprintf(pFile,"alpha: %7.3f\t beta: %7.3f\t Dmid: %8.3f\t As_mod: %7.3f\t Bs: %7.3f\t Initial_Weight: %13.5e\t Reweight_factor: %8.3f\t Quality_Factor: %7.3f, Ratios: %7.3f\t%7.3f\t%7.3f\n",alpha,beta,Dmid,As_modifier,Bs,initial_weight,reweight_factor, quality_factor,ratio1,ratio2,ratio3);
	fclose(pFile);


    std::cout << "Calculation check: Overflow flag: " << (bool)std::fetestexcept(FE_OVERFLOW) << "Underflow flag: " << (bool)std::fetestexcept(FE_UNDERFLOW) << std::endl;

    //
     if(load_save_switch == 2) {peterito->save_grid3D(initial_weight);}

return 0;
}



// NOTES
// tau_coll=n_par * sigma_par * deltaV * t_record
// w_new = w*exp(-tau_coll)
// t_pr = 400*(M_sun/M_star)*pow(r_orig/au,2)/beta
// beta = 5.7e-5/rho/particle_radius*L_star/L_sun*M_sun/M_star
