#ifndef OctreePoint_H
#define OctreePoint_H

#include "Vec3.h"

// Simple point data type to insert into the tree.
// Have something with more interesting behavior inherit
// from this in order to store other attributes in the tree.
class OctreePoint {
	Vec3 position;

	//added
	// float weight;
	// float diameter;
	// Vec3 velocity;

	// //

	// Diameter and Mass can be exchanged by just an index (interger) that will take it from a different array -> instead of two floats per record we will have one integer

	// Currently we are storing 7 floats per record = 28 bytes per record + for each OctreePoint we store 13 bytes

	// Alternative
	std::vector<float> weight;
	std::vector<float> weight_iter;
	// std::vector<float> diameter;
	// std::vector<float> mass;
	std::vector<int> SFD_index;
	std::vector<Vec3> velocity;
	int rec_num=0;

public:
	OctreePoint() { }
	// OctreePoint(const Vec3& position) : position(position) { }
	inline const Vec3& getPosition() const { return position; }
	inline float getWeight(int i) { return weight[i]; }
	inline float getWeight_Iter(int i) { return weight_iter[i]; }
	inline int getRecNum() {return rec_num;}
    inline const Vec3& getVelocity(int i) const { return velocity[i]; }
	// inline float getDiameter(int i) { return diameter[i]; }
	// inline float getMass(int i) { return mass[i]; }
	inline int getIndex(int i) {return SFD_index[i];}



	inline void setPosition(const Vec3& p) { position = p; }

	inline void insertWeight(float w) {weight.push_back(w);}
	inline void insertWeight_Iter(float w) {weight_iter.push_back(w);}
	// inline void insertDiameter(float p) { diameter.push_back(p); mass.push_back(pow(p,3.0)/6.0*M_PI*1e-18*2000.0); }
	inline void insertIndex(int p) {SFD_index.push_back(p);}
	inline void insertVelocity(const Vec3& p) { velocity.push_back(p); }
	inline void addRecord(){rec_num+=1;}
	inline void setRecord(int p) {rec_num = p;}

	inline void setWeight(float w, int i) { weight[i] = w; }
	inline void setWeight_Iter(float w, int i) { weight_iter[i] = w; }
	inline void setVelocity(const Vec3& p,int i) { velocity[i] = p; }
	inline void setIndex(int w, int i) { SFD_index[i] = w; }
	inline void printReport(int k) {printf("x = %f, y = %f, z = %f, vx = %f, vy = %f, vz = %f, weight: %f, weight_iter: %f, SFD_index: %d\n", position.x,position.y,position.z,velocity[k].x,velocity[k].y,velocity[k].z,weight[k],weight_iter[k],SFD_index[k]);}
	// inline void setDiameter(float w,int i) { diameter[i] = w; mass[i] = (pow(w,3.0)/6.0*M_PI*1e-18*2000.0);}


};

#endif
