#ifndef GridPoint_H
#define GridPoint_H

#include "Vec3.h"
using namespace std;

// Simple point data type to insert into the Grid3D.
// Currently we are storing 6 floats + 1 int per record = 28 bytes per record + for each GridPoint we store 13 bytes
// 
// Diameter and Mass can be exchanged by just an index (interger) that we take a different array (SFD_array) -> instead of two floats per record we will have one integer


class GridPoint {
	Vec3 position; // Position of the point in the grid
	std::vector<float> weight; // a float that tracks the weight BEFORE the algorithm iteration is done
	std::vector<float> weight_iter; // a float that tracks the weight AFTER the algorithm iteration is done
	std::vector<int> SFD_index; // an integer that tracks the record position in the size-frequency distribution array
	std::vector<Vec3> velocity; // three floats that provide the particle velocity in the Cartesian coordinate system
	int rec_num=0; // internal counter of number of records

public:
	GridPoint() { }
	inline const Vec3& getPosition() const { return position; }
	inline float getWeight(int i) { return weight[i]; }
	inline float getWeight_Iter(int i) { return weight_iter[i]; }
	inline int getRecNum() {return rec_num;}
    inline const Vec3& getVelocity(int i) const { return velocity[i]; }
	inline int getIndex(int i) {return SFD_index[i];}

	inline void setPosition(const Vec3& p) { position = p; }

	inline void insertWeight(float w) {weight.push_back(w);}
	inline void insertWeight_Iter(float w) {weight_iter.push_back(w);}
	inline void insertIndex(int p) {SFD_index.push_back(p);}
	inline void insertVelocity(const Vec3& p) { velocity.push_back(p); }
	inline void addRecord(){rec_num+=1;}
	inline void setRecord(int p) {rec_num = p;}

	inline void setWeight(float w, int i) { weight[i] = w; }
	inline void setWeight_Iter(float w, int i) { weight_iter[i] = w; }
	inline void setVelocity(const Vec3& p,int i) { velocity[i] = p; }
	inline void setIndex(int w, int i) { SFD_index[i] = w; }
	inline void printReport(int k) {printf("x = %f, y = %f, z = %f, vx = %f, vy = %f, vz = %f, weight: %f, weight_iter: %f, SFD_index: %d\n", position.x,position.y,position.z,velocity[k].x,velocity[k].y,velocity[k].z,weight[k],weight_iter[k],SFD_index[k]);}


};

#endif
