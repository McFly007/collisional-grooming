#ifndef Octree_H
#define Octree_H

#include <cstddef>
#include <vector>
#include "OctreePoint.h"

namespace brandonpelfrey {


	// std::vector<float> cmax_N(3);
	// std::vector<float> bmax_N(3);
	// std::vector<float> cmin_N(3);
	// std::vector<float> bmin_N(3);
	/**!
	 *
	 */
	class Octree {
		// Physical position/size. This implicitly defines the bounding 
		// box of this node
		Vec3 origin;         //! The physical center of this node
		Vec3 halfDimension;  //! Half the width/height/depth of this node

		// The tree has up to eight children and can additionally store
		// a point, though in many applications only, the leaves will store data.
		Octree *children[8]; //! Pointers to child octants
		OctreePoint *data;   //! Data point to be stored at a node

		/*
				Children follow a predictable pattern to make accesses simple.
				Here, - means less than 'origin' in that dimension, + means greater than.
				child:	0 1 2 3 4 5 6 7
				x:      - - - - + + + +
				y:      - - + + - - + +
				z:      - + - + - + - +
		 */

		public:
		Octree(const Vec3& origin, const Vec3& halfDimension) 
			: origin(origin), halfDimension(halfDimension), data(NULL) {
				// Initially, there are no children
				for(int i=0; i<8; ++i) 
					children[i] = NULL;
			}

		Octree(const Octree& copy)
			: origin(copy.origin), halfDimension(copy.halfDimension), data(copy.data) {

			}

		~Octree() {
			// Recursively destroy octants
			for(int i=0; i<8; ++i) 
				delete children[i];
		}

		// Determine which octant of the tree would contain 'point'
		int getOctantContainingPoint(const Vec3& point) const {
			int oct = 0;
			if(point.x >= origin.x) oct |= 4;
			if(point.y >= origin.y) oct |= 2;
			if(point.z >= origin.z) oct |= 1;
			return oct;
		}

		bool isLeafNode() const {
			// This is correct, but overkill. See below.
			/*
				 for(int i=0; i<8; ++i)
				 if(children[i] != NULL) 
				 return false;
				 return true;
			 */

			// We are a leaf iff we have no children. Since we either have none, or 
			// all eight, it is sufficient to just check the first.
			return children[0] == NULL;
		}

		void insert(OctreePoint* point) {
			// If this node doesn't have a data point yet assigned 
			// and it is a leaf, then we're done!
			if(isLeafNode()) {
				if(data==NULL) {
					data = point;
					return;
				} else {
					// We're at a leaf, but there's already something here
					// We will split this node so that it has 8 child octants
					// and then insert the old data that was here, along with 
					// this new data point

					// Save this data point that was here for a later re-insert
					OctreePoint *oldPoint = data;
					data = NULL;

					// Split the current node and create new empty trees for each
					// child octant.
					for(int i=0; i<8; ++i) {
						// Compute new bounding box for this child
						Vec3 newOrigin = origin;
						newOrigin.x += halfDimension.x * (i&4 ? .5f : -.5f);
						newOrigin.y += halfDimension.y * (i&2 ? .5f : -.5f);
						newOrigin.z += halfDimension.z * (i&1 ? .5f : -.5f);
						children[i] = new Octree(newOrigin, halfDimension*.5f);
					}

					// Re-insert the old point, and insert this new point
					// (We wouldn't need to insert from the root, because we already
					// know it's guaranteed to be in this section of the tree)
					children[getOctantContainingPoint(oldPoint->getPosition())]->insert(oldPoint);
					children[getOctantContainingPoint(point->getPosition())]->insert(point);
				}
			} else {
				// We are at an interior node. Insert recursively into the 
				// appropriate child octant
				int octant = getOctantContainingPoint(point->getPosition());
				children[octant]->insert(point);
			}
		}

		void initialize_Second_Weights() {
			int rec_iter;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					for (int k=0;k<rec_iter ;++k){
						data->setWeight_Iter(0.0,k);    	
					}

				}
			}

		else{
			for(int i=0; i<8; ++i) { 
					// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
						// children[i]->origin.y,
						// children[i]->origin.z);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);
				children[i]->initialize_Second_Weights();
			}

		}
	}


		void zero_Weights() {
			int rec_iter;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					for (int k=0;k<rec_iter ;++k){
						data->setWeight_Iter(0.0,k);    	
					}

				}
			}

		else{
			for(int i=0; i<8; ++i) { 
					// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
						// children[i]->origin.y,
						// children[i]->origin.z);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);
				children[i]->zero_Weights();
			}

		}
	}


	int get_tree_size(int tree_size) {
			int rec_iter;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					// for (int k=0;k<rec_iter ;++k){
							tree_size += rec_iter;   
					// }
					return tree_size;    

				}
			}

		else{
			for(int i=0; i<8; ++i) { 
					// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
						// children[i]->origin.y,
						// children[i]->origin.z);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);
				tree_size = children[i]->get_tree_size(tree_size);
			}

		}
		return tree_size;
	}

	void get_all_Weights(std::vector<float>& results_weights) {
			int rec_iter;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					for (int k=0;k<rec_iter ;++k){
						results_weights.push_back(data->getWeight(k));    	
					}

				}
			}

		else{
			for(int i=0; i<8; ++i) { 
					// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
						// children[i]->origin.y,
						// children[i]->origin.z);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);
				children[i]->get_all_Weights(results_weights);
			}

		}
	}


	void get_all_Diameter_Weights(std::vector<float>& results_weights,std::vector<float>& diameter_vector) {
			int rec_iter;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					for (int k=0;k<rec_iter ;++k){
						for(int i=0;i<(int)diameter_vector.size();i++){
							// if(float(diameter_vector[i]) ==  data->getDiameter(k)){
							if(i ==  data->getIndex(k)){
								// results_weights[i] += data->getWeight(k);
								results_weights[i] += 1.0;
							break;
						}
						} 


						// results_weights.push_back(data->getWeight(k));    	
					}

				}
			}

		else{
			for(int i=0; i<8; ++i) { 
					// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
						// children[i]->origin.y,
						// children[i]->origin.z);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);
				children[i]->get_all_Diameter_Weights(results_weights, diameter_vector);
			}

		}
	}


	void get_total_mass(float &total_mass, const std::vector<float>& diameter_vector, const std::vector < std::vector<float> > & SFD_array) {
			int rec_iter;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					for (int k=0;k<rec_iter ;++k){



						for(int i=0;i<(int)diameter_vector.size();i++){
							if(i ==  data->getIndex(k)){
								total_mass += SFD_array[i][1]*SFD_array[i][2]*data->getWeight(k);
								// results_weights[i] += 1.0;
							break;}
						// total_mass += pow(data->getDiameter(k),3)*M_PI/6.0*2000.0*1e-18*
						} 


						// results_weights.push_back(data->getWeight(k));    	
					}

				}
			}

		else{
			for(int i=0; i<8; ++i) { 
					// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
						// children[i]->origin.y,
						// children[i]->origin.z);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);
				children[i]->get_total_mass(total_mass,diameter_vector,SFD_array);
			}

		}
	}



	void tree_max_diff(float &max_diff) {
			int rec_iter;
			float difference;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					for (int k=0;k<rec_iter ;++k){
							difference = std::max(data->getWeight(k)/data->getWeight_Iter(k),data->getWeight_Iter(k)/data->getWeight(k))-1.0;
							// if(difference>1e30) printf("Weight1: %13.5e, Weight2: %13.5e\n", data->getWeight(k),data->getWeight_Iter(k));
							if(max_diff<difference && difference<1e30) { max_diff=difference;} 	
					}
				}
			}
		else{
			for(int i=0; i<8; ++i) { 
				children[i]->tree_max_diff(max_diff);
			}

		}
	}


	void check_all_Weights() {
			int rec_iter;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					for (int k=0;k<rec_iter ;++k){
						if(data->getWeight(k)<0) {printf("HOUSTON: Weight1 is %13.5e\n",data->getWeight(k)); fflush(stdout);  }  	
						if(data->getWeight_Iter(k)<0) {printf("HOUSTON: Weight2 is %13.5e\n",data->getWeight_Iter(k)); fflush(stdout);  }  	
					}

				}
			}

		else{
			for(int i=0; i<8; ++i) { 
					// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
						// children[i]->origin.y,
						// children[i]->origin.z);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);
				children[i]->check_all_Weights();
			}

		}
	}


	void set_all_Weights(const float weight) {
			int rec_iter;
			if(isLeafNode()) {
				if(data!=NULL) {
					rec_iter = data->getRecNum();
					// printf("RECITER: %d\n",rec_iter);fflush(stdout);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);fflush(stdout);
					for (int k=0;k<rec_iter ;++k){
						// data->setweight(data->getweight(k)*0.5+setting_weights[j]*0.5,k);
						if(data->getWeight(k)*(weight)+data->getWeight_Iter(k)*(1.0-weight) < 0) { printf("HOUSTON WE HAVE A PROBLEM: w1: %13.5e, w2: %13.5e, wtot: %13.5e \n", data->getWeight(k),data->getWeight_Iter(k), data->getWeight(k)*(weight)+data->getWeight_Iter(k)*(1.0-weight)); fflush(stdout);}
						data->setWeight(data->getWeight(k)*(weight)+data->getWeight_Iter(k)*(1.0-weight),k);
						
						// j+=1;
						// printf("Debug %d\n",j);
					}

				}
			}

		else{
			for(int i=0; i<8; ++i) { 
					// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
						// children[i]->origin.y,
						// children[i]->origin.z);
					// printf("LEAF HERE %.5f\n",data->getPosition().x);fflush(stdout);
				children[i]->set_all_Weights(weight);
			}

		}
	}

	// 	void set_all_weights() {
	// 		int rec_iter;
	// 		if(isLeafNode()) {
	// 			if(data!=NULL) {
	// 				rec_iter = data->getRecNum();
	// 				// printf("RECITER: %d\n",rec_iter);fflush(stdout);
	// 				// printf("LEAF HERE %.5f\n",data->getPosition().x);fflush(stdout);
	// 				for (int k=0;k<rec_iter ;++k){
	// 					// data->setweight(data->getweight(k)*0.5+setting_weights[j]*0.5,k);
	// 					data->setweight(data->getweight(k)*0.15+data->getweight_Iter(k)*0.85,k);
	// 					// j+=1;
	// 					// printf("Debug %d\n",j);
	// 				}

	// 			}
	// 		}

	// 	else{
	// 		for(int i=0; i<8; ++i) { 
	// 				// printf("%.5f, %.5f, %.5f\n", children[i]->origin.x,
	// 					// children[i]->origin.y,
	// 					// children[i]->origin.z);
	// 				// printf("LEAF HERE %.5f\n",data->getPosition().x);fflush(stdout);
	// 			children[i]->set_all_weights();
	// 		}

	// 	}
	// }

		//

		// This is a really simple routine for querying the tree for points
		// within a bounding box defined by min/max points (bmin, bmax)
		// All results are pushed into 'results'
		void getPointsInsideBox(const Vec3& bmin, const Vec3& bmax, std::vector<OctreePoint*>& results) {
			// If we're at a leaf node, just see if the current data point is inside
			// the query bounding box
			if(isLeafNode()) {
				if(data!=NULL) {
					const Vec3& p = data->getPosition();
					if(p.x>bmax.x || p.y>bmax.y || p.z>bmax.z) return;
					if(p.x<bmin.x || p.y<bmin.y || p.z<bmin.z) return;
					results.push_back(data);
				}
			} else {
				// We're at an interior node of the tree. We will check to see if
				// the query bounding box lies outside the octants of this node.
				for(int i=0; i<8; ++i) {
					// Compute the min/max corners of this child octant
					Vec3 cmax = children[i]->origin + children[i]->halfDimension;
					// Vec3 cmin = children[i]->origin - children[i]->halfDimension;

					// If the query rectangle is outside the child's bounding box, 
					// then continue
					if(cmax.x<bmin.x || cmax.y<bmin.y || cmax.z<bmin.z) continue;

					Vec3 cmin = children[i]->origin - children[i]->halfDimension;

					if(cmin.x>bmax.x || cmin.y>bmax.y || cmin.z>bmax.z) continue;

					// At this point, we've determined that this child is intersecting 
					// the query bounding box
					children[i]->getPointsInsideBox(bmin,bmax,results);
				} 
			}
		}


		// void getPointsInsideBox(const Vec3& bmin, const Vec3& bmax, std::vector<OctreePoint*>& results) {
		// 	// If we're at a leaf node, just see if the current data point is inside
		// 	// the query bounding box
		// 	if(isLeafNode()) {
		// 		if(data!=NULL) {
		// 			const Vec3& p = data->getPosition();
		// 			if(p.x>bmax.x || p.y>bmax.y || p.z>bmax.z) return;
		// 			if(p.x<bmin.x || p.y<bmin.y || p.z<bmin.z) return;
		// 			results.push_back(data);
		// 		}
		// 	} else {
		// 		// We're at an interior node of the tree. We will check to see if
		// 		// the query bounding box lies outside the octants of this node.
		// 		for(int i=0; i<8; ++i) {

		// 			float cmaxx=children[i]->origin.x + children[i]->halfDimension.x;if(cmaxx<bmin.x) continue;
		// 			float cmaxy=children[i]->origin.y + children[i]->halfDimension.y;if(cmaxy<bmin.y) continue;
		// 			float cmaxz=children[i]->origin.z + children[i]->halfDimension.z;if(cmaxz<bmin.z) continue;

		// 			float cminx=children[i]->origin.x - children[i]->halfDimension.x;if(cminx>bmax.x) continue;
		// 			float cminy=children[i]->origin.y - children[i]->halfDimension.y;if(cminy>bmax.y) continue;
		// 			float cminz=children[i]->origin.z - children[i]->halfDimension.z;if(cminz>bmax.z) continue;

		// 			// // Compute the min/max corners of this child octant
		// 			// Vec3 cmax = children[i]->origin + children[i]->halfDimension;
		// 			// // Vec3 cmin = children[i]->origin - children[i]->halfDimension;

		// 			// // If the query rectangle is outside the child's bounding box, 
		// 			// // then continue
		// 			// if(cmax.x<bmin.x || cmax.y<bmin.y || cmax.z<bmin.z) continue;

		// 			// Vec3 cmin = children[i]->origin - children[i]->halfDimension;
					
		// 			// if(cmin.x>bmax.x || cmin.y>bmax.y || cmin.z>bmax.z) continue;

		// 			// At this point, we've determined that this child is intersecting 
		// 			// the query bounding box
		// 			children[i]->getPointsInsideBox(bmin,bmax,results);
		// 		} 
		// 	}
		// }


		// void getPointsInsideBox(const std::vector<float>& bmin_N, const std::vector<float>& bmax_N, std::vector<OctreePoint*>& results) {
		// 	// If we're at a leaf node, just see if the current data point is inside
		// 	// the query bounding box
		// 	if(isLeafNode()) {
		// 		if(data!=NULL) {
		// 			const Vec3& p = data->getPosition();
		// 			if(p.x>bmax_N[0] || p.y>bmax_N[1] || p.z>bmax_N[2]) return;
		// 			if(p.x<bmin_N[0] || p.y<bmin_N[1] || p.z<bmin_N[2]) return;
		// 			results.push_back(data);
		// 		}
		// 	} else {
		// 		// We're at an interior node of the tree. We will check to see if
		// 		// the query bounding box lies outside the octants of this node.
		// 		for(int i=0; i<8; ++i) {
		// 			// Compute the min/max corners of this child octant
		// 			cmax_N[0] = children[i]->origin.x + children[i]->halfDimension.x;
		// 			cmax_N[1] = children[i]->origin.y + children[i]->halfDimension.y;
		// 			cmax_N[2] = children[i]->origin.z + children[i]->halfDimension.z;
		// 			cmin_N[0] = children[i]->origin.x- children[i]->halfDimension.x;
		// 			cmin_N[1] = children[i]->origin.y- children[i]->halfDimension.y;
		// 			cmin_N[2] = children[i]->origin.z- children[i]->halfDimension.z;

		// 			// If the query rectangle is outside the child's bounding box, 
		// 			// then continue
		// 			// if(cmax.x<bmin.x || cmax.y<bmin.y || cmax.z<bmin.z) continue;
		// 			// if(cmin.x>bmax.x || cmin.y>bmax.y || cmin.z>bmax.z) continue;
		// 			if(cmax_N[0]<bmin_N[0] || cmax_N[1]<bmin_N[1] || cmax_N[2]<bmin_N[2]) continue;
		// 			if(cmin_N[0]>bmax_N[0] || cmin_N[1]>bmax_N[1] || cmin_N[2]>bmax_N[2]) continue;

		// 			// At this point, we've determined that this child is intersecting 
		// 			// the query bounding box
		// 			children[i]->getPointsInsideBox(bmin_N,bmax_N,results);
		// 		} 
		// 	}
		// }

	};

}
#endif
