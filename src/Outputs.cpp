#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Grid3D.h"
#include "Outputs.h"
#include "Settings.h"
#include "Records.h"



// Define namespaces
using namespace petrpokorny;
using namespace std;

Outputs My_Outputs;

vector<double> r_prof_init; // Vector that keeps the radial profile for the 1st iteration
vector<double> r_prof_iter; // Vector that keeps the radial profile for the 2nd iteration
vector<double> r_prof_brightness; // Vector that keeps the radial profile for the 2nd iteration
vector<double> ecl_profile; // Vector that keeps the ecliptic latitude profile
vector<double> ecl_profile_mass; // Vector that keeps the ecliptic latitude profile

void Initiate_Output_Files()
{
	FILE *pFile_temp;
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Face_On_Number", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Face_On_Brightness", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Face_On_Impact_Velocity", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Edge_On_Number", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Edge_On_Brightness", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Edge_On_Impact_Velocity", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Mass_Accreted_Kessler", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Face_On_MassLoss", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Edge_On_MassLoss", "w+");
	::fclose(pFile_temp);
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Radial_MassLoss", "w+");
	::fclose(pFile_temp);
}

void Initialize_1D_Array(std::vector<double> &init_array, int xdim)
{
	init_array.resize(xdim);
}

void Initialize_2D_Array(std::vector<std::vector<double>> &init_array, int xdim, int ydim)
{
	init_array.resize(xdim);
	for (int i = 0; i < xdim; i++)
	{
		init_array[i].resize(ydim);
	}
}

void Initialize_Outputs(Grid3D *grid3Ddata)
{
    Initialize_1D_Array(r_prof_init,10000);		 // Vector that keeps the radial profile for the 1st iteration
	Initialize_1D_Array(r_prof_iter, 10000);		 // Vector that keeps the radial profile for the 2nd iteration
	Initialize_1D_Array(r_prof_brightness, 10000); // Vector that keeps the radial profile for the 2nd iteration
	Initialize_1D_Array(ecl_profile, 180);		 // Vector that keeps the ecliptic latitude profile
	Initialize_1D_Array(ecl_profile_mass, 180);	 // Vector that keeps the ecliptic latitude profile

	// TO DO - This works for now but could be better
	int xdim = grid3Ddata->data.size();
	int ydim = grid3Ddata->data[0].size();

	
	Initialize_2D_Array(My_Outputs.disk_profile_face_on_number, xdim, ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_face_on_brightness, xdim, ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_face_on_impact_velocity, xdim, ydim);

	Initialize_2D_Array(My_Outputs.disk_profile_edge_on_number, xdim, ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_edge_on_brightness, xdim, ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_edge_on_impact_velocity, xdim, ydim);

	Initialize_2D_Array(My_Outputs.disk_profile_face_on_massloss, xdim, ydim);
	Initialize_2D_Array(My_Outputs.disk_profile_edge_on_massloss, xdim, ydim);
	Initialize_1D_Array(My_Outputs.disk_profile_radial_massloss, 10000);

	::cout << "Check the Output Size:" << std::endl;
    ::cout << "1D R_hel pre-iter Profile size: " << r_prof_init.size() << std::endl;
    ::cout << "1D R_hel post-iter Profile size: " << r_prof_iter.size() << std::endl;
    ::cout << "1D R_hel brightness Profile size: " << r_prof_brightness.size() << std::endl;
    ::cout << "1D R_hel ecliptic latitude Profile size: " << ecl_profile.size() << std::endl;
    ::cout << "1D R_hel ecliptic latitude mass Profile size: " << ecl_profile_mass.size() << std::endl;
    ::cout << "------------------" << std::endl;
	::cout << "2D Profile 1 size: " << My_Outputs.disk_profile_face_on_number.size() << "\t" << My_Outputs.disk_profile_face_on_number[0].size() << std::endl;
	::cout << "2D Profile 2 size: " << My_Outputs.disk_profile_face_on_brightness.size() << "\t" << My_Outputs.disk_profile_face_on_brightness[0].size() << std::endl;
	::cout << "2D Profile 3 size: " << My_Outputs.disk_profile_face_on_impact_velocity.size() << "\t" << My_Outputs.disk_profile_face_on_impact_velocity[0].size() << std::endl;
	::cout << "2D Profile 4 size: " << My_Outputs.disk_profile_edge_on_number.size() << "\t" << My_Outputs.disk_profile_edge_on_number[0].size() << std::endl;
	::cout << "2D Profile 5 size: " << My_Outputs.disk_profile_edge_on_brightness.size() << "\t" << My_Outputs.disk_profile_edge_on_brightness[0].size() << std::endl;
	::cout << "2D Profile 6 size: " << My_Outputs.disk_profile_edge_on_impact_velocity.size() << "\t" << My_Outputs.disk_profile_edge_on_impact_velocity[0].size() << std::endl;
}

void Print_Outputs_Big(Grid3D *grid3Ddata, int file_counter, double refactor)
{
	int x = grid3Ddata->data.size();
	int y = grid3Ddata->data[0].size();
	int z = grid3Ddata->data[0][0].size();

	FILE *pFile_temp;

	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Face_On_Number", "a+");
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			if (My_Outputs.disk_profile_face_on_number[i][j] > 0.0)
			{
				::fprintf(pFile_temp, "%.5f %.5f %13.5e %d \n", settings_dr[0] *0.5 + i * settings_dr[0] - settings_grid_half_size, settings_dr[1] *0.5 + j * settings_dr[0] - settings_grid_half_size, My_Outputs.disk_profile_face_on_number[i][j] * refactor, file_counter);
			}
		}
	};
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);

	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Face_On_Brightness", "a+");
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			if (My_Outputs.disk_profile_face_on_brightness[i][j] > 0.0)
			{
				::fprintf(pFile_temp, "%.5f %.5f %13.5e %d \n", settings_dr[0] *0.5 + i * settings_dr[0] - settings_grid_half_size, settings_dr[1] *0.5 + j * settings_dr[0] - settings_grid_half_size, My_Outputs.disk_profile_face_on_brightness[i][j] * refactor, file_counter);
			}
		}
	};
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);

	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Face_On_Impact_Velocity", "a+");
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			if (My_Outputs.disk_profile_face_on_number[i][j] > 0.0)
			{
				::fprintf(pFile_temp, "%.5f %.5f %13.5e %d \n", settings_dr[0] *0.5 + i * settings_dr[0] - settings_grid_half_size, settings_dr[1] *0.5 + j * settings_dr[0] - settings_grid_half_size, My_Outputs.disk_profile_face_on_impact_velocity[i][j] * refactor, file_counter);
			}
		}
	};
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);

	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Edge_On_Number", "a+");
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < z; ++j)
		{
			if (My_Outputs.disk_profile_edge_on_number[i][j] > 0.0)
			{
				::fprintf(pFile_temp, "%.5f %.5f %13.5e %d \n", settings_dr[0] *0.5 + i * settings_dr[0] - settings_grid_half_size, settings_dr[1] *0.5 + j * settings_dr[0] - settings_grid_half_size, My_Outputs.disk_profile_edge_on_number[i][j] * refactor, file_counter);
			}
		};
	}
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);

	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Edge_On_Brightness", "a+");
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < z; ++j)
		{
			if (My_Outputs.disk_profile_edge_on_brightness[i][j] > 0.0)
			{
				::fprintf(pFile_temp, "%.5f %.5f %13.5e %d \n", settings_dr[0] *0.5 + i * settings_dr[0] - settings_grid_half_size, settings_dr[1] *0.5 + j * settings_dr[0] - settings_grid_half_size, My_Outputs.disk_profile_edge_on_brightness[i][j] * refactor, file_counter);
			}
		}
	};
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);

	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Edge_On_Impact_Velocity", "a+");
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < z; ++j)
		{
			if (My_Outputs.disk_profile_edge_on_number[i][j] > 0.0)
			{
				::fprintf(pFile_temp, "%.5f %.5f %13.5e %d \n", settings_dr[0] *0.5 + i * settings_dr[0] - settings_grid_half_size, settings_dr[1] *0.5 + j * settings_dr[0] - settings_grid_half_size, My_Outputs.disk_profile_edge_on_impact_velocity[i][j] * refactor, file_counter);
			}
		}
	};
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);

	///
	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Face_On_MassLoss", "a+");
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			if (My_Outputs.disk_profile_face_on_massloss[i][j] > 0.0)
			{
				::fprintf(pFile_temp, "%.5f %.5f %13.5e %d \n", settings_dr[0] *0.5 + i * settings_dr[0] - settings_grid_half_size, settings_dr[1] *0.5 + j * settings_dr[0] - settings_grid_half_size, My_Outputs.disk_profile_face_on_massloss[i][j] * refactor, file_counter);
			}
		}
	};
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);

	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Edge_On_MassLoss", "a+");
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < z; ++j)
		{
			if (My_Outputs.disk_profile_edge_on_massloss[i][j] > 0.0)
			{
				::fprintf(pFile_temp, "%.5f %.5f %13.5e %d \n", settings_dr[0] *0.5 + i * settings_dr[0] - settings_grid_half_size, settings_dr[1] *0.5 + j * settings_dr[0] - settings_grid_half_size, My_Outputs.disk_profile_edge_on_massloss[i][j] * refactor, file_counter);
			}
		}
	};
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);

	pFile_temp = fopen("./Outputs/Outputs_Disk_Profile_Radial_MassLoss", "a+");
	// for (int i=0;i<(int)My_Outputs.disk_profile_radial_massloss.size(); ++i){if(My_Outputs.disk_profile_radial_massloss[i] >0.0) {::fprintf(pFile_temp,"%.5f %13.5e %d\n", i*dr.x*2.0+dr.x,My_Outputs.disk_profile_radial_massloss[i]*refactor, file_counter );}	} ;  ::fprintf(pFile_temp,"\n");	::fclose(pFile_temp);
	for (int i = 0; i < (int)My_Outputs.disk_profile_radial_massloss.size(); ++i)
	{
		if (My_Outputs.disk_profile_radial_massloss[i] > 0.0)
		{
			::fprintf(pFile_temp, "%.5f %13.5e %d\n", (i + 0.5) * (settings_grid_half_size) / (float)My_Outputs.disk_profile_radial_massloss.size(), My_Outputs.disk_profile_radial_massloss[i] * refactor, file_counter);
		}
	};
	::fprintf(pFile_temp, "\n");
	::fclose(pFile_temp);
}

void Print_Outputs_Small(Grid3D *grid3Ddata, int file_counter, double refactor)
{

	// int x = grid3Ddata->data.size();
	// int y = grid3Ddata->data[0].size();
	// int z = grid3Ddata->data[0][0].size();

	FILE *pFile_temp;

	pFile_temp = fopen("./Outputs/Outputs_Mass_Accreted_Kessler", "a+");
	::fprintf(pFile_temp, "%13.5e %d \n", My_Outputs.mass_accreted_at_Earth * refactor, file_counter);
	::fclose(pFile_temp);
}

void Zero_Outputs(Grid3D *grid3Ddata)
{

	int x = grid3Ddata->data.size();
	int y = grid3Ddata->data[0].size();
	int z = grid3Ddata->data[0][0].size();

	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			My_Outputs.disk_profile_face_on_number[i][j] = 0.0;
		}
	}
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			My_Outputs.disk_profile_face_on_brightness[i][j] = 0.0;
		}
	}
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			My_Outputs.disk_profile_face_on_impact_velocity[i][j] = 0.0;
		}
	}
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < z; ++j)
		{
			My_Outputs.disk_profile_edge_on_number[i][j] = 0.0;
		}
	}
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < z; ++j)
		{
			My_Outputs.disk_profile_edge_on_brightness[i][j] = 0.0;
		}
	}
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < z; ++j)
		{
			My_Outputs.disk_profile_edge_on_impact_velocity[i][j] = 0.0;
		}
	}

	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			My_Outputs.disk_profile_face_on_massloss[i][j] = 0.0;
		}
	}
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < z; ++j)
		{
			My_Outputs.disk_profile_edge_on_massloss[i][j] = 0.0;
		}
	}
	for (int i = 0; i < (int)My_Outputs.disk_profile_radial_massloss.size(); ++i)
	{
		My_Outputs.disk_profile_radial_massloss[i] = 0.0;
	}
}