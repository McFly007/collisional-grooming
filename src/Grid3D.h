#ifndef Grid3D_H
#define Grid3D_H

#include <cstddef>
#include <vector>
#include "GridPoint.h"

using namespace std;

// This file contains all function related to the Grid3D class

namespace petrpokorny
{

	class Grid3D
	{

		Vec3 xyzsize; // Physical size of the grid - it has three components (x,y,z)
		Vec3 dr;	  // The size of one bin - it has three components (x,y,z)

		int xdim; // Number of bins along the x-axis
		int ydim; // Number of bins along the y-axis
		int zdim; // Number of bins along the z-axis

		// Function for rounding doubles that returns an integer
		int round(double d)
		{
			return d < 0 ? ceil(d - 0.5) : floor(d + 0.5);
		}

	public:
		vector<vector<vector<GridPoint>>> data; // The 3D representation of the grid

		Grid3D(const Vec3 &xyzsize, const Vec3 &dr)
			: xyzsize(xyzsize), dr(dr)
		{

			xdim = round(xyzsize.x / dr.x) * 2;
			ydim = round(xyzsize.y / dr.y) * 2;
			zdim = round(xyzsize.z / dr.z) * 2;

			data.resize(xdim);

			// Give the data the correct size
			for (int i = 0; i < (int)data.size(); i++)
			{
				data[i].resize(ydim);
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					data[i][j].resize(zdim);
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						data[i][j][k].setPosition(Vec3(-xyzsize.x + 0.5 * dr.x + dr.x * i, -xyzsize.y + 0.5 * dr.y + dr.y * j, -xyzsize.z + 0.5 * dr.z + dr.z * k)); // Set the grid positions - otherwise the points are empty
					}
				}
			}

			::cout << "--- 3D Grid Successfully created ---"
				   << "\n";
			::cout << "Data x-dim size: " << data.size() << "\n";
			::cout << "Data y-dim size: " << data[0].size() << "\n";
			::cout << "Data z-dim size: " << data[0][0].size() << "\n";
		}

		// Get grid dimensions
		int getXDim() { return xdim; }
		int getYDim() { return ydim; }
		int getZDim() { return zdim; }

		// Redundant - to be removed
		// void getPointsInsideGrid(const Vec3 point, GridPoint *results)
		// {
		// 	int xindex = round(point.x / dr.x) + xdim * 0.5;
		// 	int yindex = round(point.y / dr.y) + ydim * 0.5;
		// 	int zindex = round(point.z / dr.z) + zdim * 0.5;

		// 	cout << data[xindex][yindex][zindex].getRecNum();
		// 	results = &data[xindex][yindex][zindex];
		// 	return;
		// }

		// Returns Vector3 of indices of a given point - used to lookup values in the 3D grid
		std::vector<int> getIndexes(const Vec3 point)
		{
			int xindex = floor(point.x / dr.x) + xdim * 0.5;
			int yindex = floor(point.y / dr.y) + ydim * 0.5;
			int zindex = floor(point.z / dr.z) + zdim * 0.5;

			std::vector<int> out;
			out.push_back(xindex);
			out.push_back(yindex);
			out.push_back(zindex);
			return out;
		}

		// Checks if all weights (pre- and post-iteration) in the 3D grid are not negative
		// Negative weights should never happen and the code is then stopped
		void check_all_Weights()
		{
			int rec_iter;
			for (int i = 0; i < (int)data.size(); i++)
			{
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						rec_iter = data[i][j][k].getRecNum();
						for (int kk = 0; kk < rec_iter; ++kk)
						{
							if (data[i][j][k].getWeight(kk) < 0)
							{
								printf("Aborting program: Problem with checking weights in the 3D grid: Pre-iter weight is %13.5e and should be non-negative\n\
								Debug info: x = %d, y = %d, z = %d, record number = %d", data[i][j][k].getWeight(kk), i,j,k,kk);
								fflush(stdout);
								abort();
							}
							if (data[i][j][k].getWeight_Iter(kk) < 0)
							{
								printf("Aborting program: Problem with checking weights in the 3D grid: Post-iter weight is %13.5e and should be non-negative\n\
								Debug info: x = %d, y = %d, z = %d, record number = %d", data[i][j][k].getWeight_Iter(kk), i,j,k,kk);
								fflush(stdout);
								abort();
							}
						}
					}
				}
			}
		}

		// Sets a weighted average of pre- and post-iteration weights to the pre-iter weights
		// Used to blend the pre- and post-iteration states of the 3D grid to make the convergence smoother
		// Formula: pre-weight = pre-weight* D_weight + post-weight * (1 - D_weight)
		void averagePreIterWeights(const float D_weight)
		{
			int rec_iter;
			for (int i = 0; i < (int)data.size(); i++)
			{
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						rec_iter = data[i][j][k].getRecNum();

						for (int kk = 0; kk < rec_iter; ++kk)
						{
							float new_weight = data[i][j][k].getWeight(kk) * (D_weight) + data[i][j][k].getWeight_Iter(kk) * (1.0 - D_weight);
							if (new_weight < 0)
							{
								printf("Aborting program: Problem with setting weights between iterations: w1: %13.5e, w2: %13.5e, wtot: %13.5e \n\
								Debug info: x = %d, y = %d, z = %d, record number = %d", \
								data[i][j][k].getWeight_Iter(kk),data[i][j][k].getWeight(kk), new_weight, i, j, k, kk);
								fflush(stdout);
								abort();
							}
							
							data[i][j][k].setWeight(new_weight, kk);
						}

					}
				}
			}
		}

		// Assigns zeros to all post-iteration record weights
		void zero_Iter_Weights()
		{
			int rec_iter;

			for (int i = 0; i < (int)data.size(); i++)
			{
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						rec_iter = data[i][j][k].getRecNum();
						for (int kk = 0; kk < rec_iter; ++kk)
						{
							data[i][j][k].setWeight_Iter(0.0, kk);
						}
					}
				}
			}
		}

		// Assigns zeros to all post-iteration record weights
		void zero_ALL_Weights()
		{
			int rec_iter;

			for (int i = 0; i < (int)data.size(); i++)
			{
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						rec_iter = data[i][j][k].getRecNum();
						for (int kk = 0; kk < rec_iter; ++kk)
						{
							data[i][j][k].setWeight_Iter(0.0, kk);
							data[i][j][k].setWeight(0.0, kk);
						}
					}
				}
			}
		}

		// Multiplies all weights by a factor of "reweight_factor"
		void multiplyAllWeights(float reweight_factor)
		{
			int rec_iter;
			for (int i = 0; i < (int)data.size(); i++)
			{
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						rec_iter = data[i][j][k].getRecNum();
						for (int kk = 0; kk < rec_iter; ++kk)
						{
							data[i][j][k].setWeight_Iter(data[i][j][k].getWeight_Iter(kk) * reweight_factor, kk);
							data[i][j][k].setWeight(data[i][j][k].getWeight(kk) * reweight_factor, kk);
						}
					}
				}
			}
		}

		// Calculates the maximum logarithmic difference of all weights in the 3D grid - used to determine if the particle cloud is converging to a stable solution
		float get_the_iteration_difference()
		{
			int rec_iter;
			double log_diff;
			double max_log_diff = 0.0;
			for (int i = 0; i < (int)data.size(); i++)
			{
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						rec_iter = data[i][j][k].getRecNum();
						for (int kk = 0; kk < rec_iter; ++kk)
						{

							if (abs(i - 0.5 * (int)data.size()) > 2)
							{
								if (abs(j - 0.5 * (int)data.size()) > 2)
								{
									if (abs(k - 0.5 * (int)data.size()) > 2)
									{
										// If data is not empty - it should never happen but we still check
										if (data[i][j][k].getWeight(kk) > 0.0 && data[i][j][k].getWeight_Iter(kk) > 0.0)
										{
											log_diff = abs(log10(data[i][j][k].getWeight(kk) / data[i][j][k].getWeight_Iter(kk)));
										}

										// Looking for the maximum logarithmic differnce
										if (log_diff > max_log_diff)
										{
											max_log_diff = log_diff;
											// For logging if necessary
											// std::cout << "Log: Maximum Difference Found: x = " << i << " y = " << j << " z = " << k << " Size index: " << data[i][j][k].getIndex(kk) << " Weights: " << data[i][j][k].getWeight(kk) << "\t" << data[i][j][k].getWeight_Iter(kk) << std::endl;
										}
									}
								}
							}
						}
					}
				}
			}
			return max_log_diff;
		}


		// Read a record from the saved 3D Grid binary file
		void read_grid_record(ifstream &file, float &weight, float &weight_iter, int &SFD_index, Vec3 &velocity)
		{
			// char * memblock_id = new char [4]; // each int has 4 bytes
			char *memblock = new char[24]; // each float has 4 bytes
			float *read_values;
			int *read_values_int;

			file.read(memblock, 8);
			read_values = (float *)memblock; // reinterpret as doubles

			weight = read_values[0];
			weight_iter = read_values[1];

			file.read(memblock, 4);
			read_values_int = (int *)memblock; // reinterpret as doubles
			SFD_index = read_values_int[0];

			file.read(memblock, 12);
			read_values = (float *)memblock; // reinterpret as double

			velocity = Vec3(read_values[0], read_values[1], read_values[2]);
			delete memblock;
		}

		// Saves grid record to the stream in &file - used to save 3D Grid into a binary file
		void save_grid_record(ofstream &file, float weight, float weight_iter, int SFD_index, Vec3 velocity)
		{
			file.write((char *)&weight, sizeof(weight));
			file.write((char *)&weight_iter, sizeof(weight_iter));
			file.write((char *)&SFD_index, sizeof(SFD_index));
			file.write((char *)&velocity, sizeof(velocity));
		}

		// Saves the entire 3D grid file into data_filename
		// First we save the initial weight of the simulation and then the entire grid in a binary form
		void save_grid3D(string data_filename, double initial_weight)
		{
			ofstream file(data_filename, ios::out | ios::binary);
			if (!file)
			{
				cout << "Function: save_grid3d - Cannot open file: " << data_filename << endl;
				abort();
			}

			file.write((char *)&initial_weight, sizeof(initial_weight));
			for (int i = 0; i < (int)data.size(); i++)
			{
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						// now we are cycling through every cell in the 3D Grid
						int rec_iter = data[i][j][k].getRecNum();
						Vec3 position = data[i][j][k].getPosition();
						file.write((char *)&rec_iter, sizeof(rec_iter));
						file.write((char *)&position, sizeof(position));

						for (int kk = 0; kk < rec_iter; ++kk)
						{
							float weight = data[i][j][k].getWeight(kk);
							float weight_iter = data[i][j][k].getWeight_Iter(kk);
							int SFD_index = data[i][j][k].getIndex(kk);
							Vec3 velocity = data[i][j][k].getVelocity(kk);
							save_grid_record(file, weight, weight_iter, SFD_index, velocity);
						}
					}
				}
			}
			file.close();
		}

		// Loads the entire 3D grid file from data_filename
		double load_grid3D(string data_filename)
		{
			char pink;
			int record_counter = 0;
			double initial_weight;
			ifstream file(data_filename, ios::in | ios::binary);
			
			if (!file)
			{
				cout << "Function: load_grid3D - Cannot open file: " << data_filename << endl;
				abort();
			}

			file.read((char *)&initial_weight, sizeof(initial_weight));
			for (int i = 0; i < (int)data.size(); i++)
			{
				for (int j = 0; j < (int)data[0].size(); j++)
				{
					for (int k = 0; k < (int)data[0][0].size(); k++)
					{
						int rec_iter;
						Vec3 position;

						file.read((char *)&rec_iter, sizeof(rec_iter));
						file.read((char *)&position, sizeof(position));

						data[i][j][k].setPosition(position);

						for (int kk = 0; kk < rec_iter; ++kk)
						{
							float weight;
							float weight_iter;
							int SFD_index;
							Vec3 velocity;
							read_grid_record(file, weight, weight_iter, SFD_index, velocity);

							data[i][j][k].insertWeight(weight);
							data[i][j][k].insertWeight_Iter(weight_iter);
							data[i][j][k].insertIndex(SFD_index);
							data[i][j][k].insertVelocity(velocity);
							data[i][j][k].addRecord();
						}
					}
				}
			}
			file.close();
			return initial_weight;
		}
	};

}
#endif
