#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Grid3D.h"
#include "Outputs.h"
#include "Settings.h"
#include "Records.h"
#include "FreeParameters.h"

int find_diameter_index(float diameter)
{
	// Finds the diameter index in the SFD_array
	for (int i = 0; i < (int)SFD_array.size(); i++)
	{
		if (float(SFD_array[i][0]) == diameter)
		{
			return i;
		}
	}
	::cout << "Error finding dimameters, something is terribly wrong, this should not happen at all! --- function int find_diameter_index(float diameter)\n";
	::exit(EXIT_FAILURE);
	return 999;
}