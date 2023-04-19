#ifndef Grid3D_H
#define Grid3D_H

#include <cstddef>
#include <vector>
#include "GridPoint.h"

using namespace std;

namespace petrpokorny {


	// std::vector<float> cmax_N(3);
	// std::vector<float> bmax_N(3);
	// std::vector<float> cmin_N(3);
	// std::vector<float> bmin_N(3);
	/**!
	 *
	 */
	class Grid3D {
		// Physical position/size. This implicitly defines the bounding 
		// box of this node
		Vec3 xyzsize;         //! The physical center of this node
		Vec3 dr;  //! Half the width/height/depth of this node
		

		
		int xdim;
		int ydim;
		int zdim;
		/*
				Children follow a predictable pattern to make accesses simple.
				Here, - means less than 'origin' in that dimension, + means greater than.
				child:	0 1 2 3 4 5 6 7
				x:      - - - - + + + +
				y:      - - + + - - + +
				z:      - + - + - + - +
		 */

		double round(double d)
		{
    	return d < 0 ? ceil(d - 0.5) : floor(d + 0.5);
		}

		public:
			vector<vector<vector<GridPoint>>> data;   //! Data point to be stored at a node
		// 	vector<vector<vector<GridPoint>>>  data;   //! Data point to be stored at a node
		
		Grid3D(const Vec3& xyzsize, const Vec3& dr) 
			: xyzsize(xyzsize), dr(dr){

				 xdim = round(xyzsize.x/dr.x)*2;
				 ydim = round(xyzsize.y/dr.y)*2;
				 zdim = round(xyzsize.z/dr.z)*2;
	
		// 		data(xdim, std::vector<std::vector<GridPoint> >(ydim,  std::vector<GridPoint>(zdim)));



  // vector < vector < vector<int> > > tube;
  // for(int i = 0; i < 2; i++)
  // {
  //   vector < vector < int > > w;
  //   tube.push_back( w );
  //   for(int j = 0; j < 4; j++)
  //   {
  //     vector <int> v;
  //     tube[i].push_back( v );
  //     for(int k = 0; k < 15; k++)
  //     {
  //       tube[i][j].push_back( rand());
  //     }
  //   }
  // }



				data.resize(xdim);

				for(int i=0; i<(int)data.size(); i++){
    				data[i].resize(ydim);
					for(int j=0; j<(int)data[0].size(); j++){
						data[i][j].resize(zdim);
						for(int k=0; k<(int)data[0][0].size(); k++){
							data[i][j][k].setPosition(Vec3(-xyzsize.x+0.5*dr.x+dr.x*i,-xyzsize.y+0.5*dr.y+dr.y*j,-xyzsize.z+0.5*dr.z+dr.z*k)); // Set the grid positions - otherwise the points are empty 
							// cout << i << "\t" << j << "\t" << k <<"\t" << data[i][j][k].getPosition().x <<"\t" << data[i][j][k].getPosition().y <<"\t" << data[i][j][k].getPosition().z << "\t" << data[i][j][k].getRecNum() << "\n" ;
						}	
					}
				}
				// cout << "Data size: " <<  data.size() << "\n";
				cout << "Data xdim size: " <<  data.size() << "\n";
				cout << "Data ydim size: " <<  data[0].size() << "\n";
				cout << "Data zdim size: " <<  data[0][0].size() << "\n";

			}

		int getXDim() { return xdim;}
		int getYDim() { return ydim;}
		int getZDim() { return zdim;}


	void getPointsInsideGrid(const Vec3 point, GridPoint *results) {
		int xindex = round(point.x/dr.x)+xdim*0.5;
		int yindex = round(point.y/dr.y)+ydim*0.5;
		int zindex = round(point.z/dr.z)+zdim*0.5;
		// cout << xindex << "\t" << yindex << "\t" << zindex << "\t" << point.x << "\t" << point.y << "\t" << point.z << "\n";
		// char pinkie;
		// cin >> pinkie;
		cout << data[xindex][yindex][zindex].getRecNum();
		results = &data[xindex][yindex][zindex];
		return;
	}


	std::vector<int> getIndexes(const Vec3 point) {
		int xindex = floor(point.x/dr.x)+xdim*0.5;
		int yindex = floor(point.y/dr.y)+ydim*0.5;
		int zindex = floor(point.z/dr.z)+zdim*0.5;

		std::vector<int> out;
		// cout << xindex << "\t" << yindex << "\t" << zindex << "\t" << point.x << "\t" << point.y << "\t" << point.z << "\n";
		// char pinkie;
		// cin >> pinkie;
		// cout << data[xindex][yindex][zindex].getRecNum();
		out.push_back(xindex);
		out.push_back(yindex);
		out.push_back(zindex);
		return out;
	}


void check_all_Weights() {
			int rec_iter;
			char pinkie;
			for(int i=0; i<(int)data.size(); i++){
					for(int j=0; j<(int)data[0].size(); j++){
						for(int k=0; k<(int)data[0][0].size(); k++){
							rec_iter = data[i][j][k].getRecNum();
							for (int kk=0;kk<rec_iter ;++kk){
						if(data[i][j][k].getWeight(kk)<0) {printf("HOUSTON: Weight1 is %13.5e\n",data[i][j][k].getWeight(kk)); fflush(stdout); cin >> pinkie; }  	
						if(data[i][j][k].getWeight_Iter(kk)<0) {printf("HOUSTON: Weight2 is %13.5e\n",data[i][j][k].getWeight_Iter(kk)); fflush(stdout); cin >> pinkie; }  	
							// cout << i << "\t" << j << "\t" << k <<"\t" << data[i][j][k].getPosition().x <<"\t" << data[i][j][k].getPosition().y <<"\t" << data[i][j][k].getPosition().z << "\t" << data[i][j][k].getRecNum() << "\n" ;
						}	
					}
				}
				}

	}


	void set_all_Weights(const float weight) {
			int rec_iter;
			char pinkie;
			for(int i=0; i<(int)data.size(); i++){
					for(int j=0; j<(int)data[0].size(); j++){
						for(int k=0; k<(int)data[0][0].size(); k++){
							rec_iter = data[i][j][k].getRecNum();
							for (int kk=0;kk<rec_iter ;++kk){
						if(data[i][j][k].getWeight(kk)*(weight)+data[i][j][k].getWeight_Iter(kk)*(1.0-weight) < 0) { printf("HOUSTON WE HAVE A PROBLEM: w1: %13.5e, w2: %13.5e, wtot: %13.5e \n", data[i][j][k].getWeight(kk),data[i][j][k].getWeight_Iter(kk), data[i][j][k].getWeight(kk)*(weight)+data[i][j][k].getWeight_Iter(kk)*(1.0-weight)); fflush(stdout); cin >> pinkie;}
						data[i][j][k].setWeight(data[i][j][k].getWeight(kk)*(weight)+data[i][j][k].getWeight_Iter(kk)*(1.0-weight),kk);
							// cout << i << "\t" << j << "\t" << k <<"\t" << data[i][j][k].getPosition().x <<"\t" << data[i][j][k].getPosition().y <<"\t" << data[i][j][k].getPosition().z << "\t" << data[i][j][k].getRecNum() << "\n" ;
						}	
					}
				}
				}

	}


	void zero_Weights() {
			int rec_iter;

			for(int i=0; i<(int)data.size(); i++){
					for(int j=0; j<(int)data[0].size(); j++){
						for(int k=0; k<(int)data[0][0].size(); k++){
							rec_iter = data[i][j][k].getRecNum();
							for (int kk=0;kk<rec_iter ;++kk){
						data[i][j][k].setWeight_Iter(0.0,kk); 
							// cout << i << "\t" << j << "\t" << k <<"\t" << data[i][j][k].getPosition().x <<"\t" << data[i][j][k].getPosition().y <<"\t" << data[i][j][k].getPosition().z << "\t" << data[i][j][k].getRecNum() << "\n" ;
						}	
					}
				}
				}

	}

	void zero_ALL_Weights() {
			int rec_iter;

			for(int i=0; i<(int)data.size(); i++){
					for(int j=0; j<(int)data[0].size(); j++){
						for(int k=0; k<(int)data[0][0].size(); k++){
							rec_iter = data[i][j][k].getRecNum();
							for (int kk=0;kk<rec_iter ;++kk){
						data[i][j][k].setWeight_Iter(0.0,kk); 
						data[i][j][k].setWeight(0.0,kk); 
							// cout << i << "\t" << j << "\t" << k <<"\t" << data[i][j][k].getPosition().x <<"\t" << data[i][j][k].getPosition().y <<"\t" << data[i][j][k].getPosition().z << "\t" << data[i][j][k].getRecNum() << "\n" ;
						}	
					}
				}
				}

	}

	void reweight_ALL_Weights(float reweight_factor) {
			int rec_iter;

			for(int i=0; i<(int)data.size(); i++){
					for(int j=0; j<(int)data[0].size(); j++){
						for(int k=0; k<(int)data[0][0].size(); k++){
							rec_iter = data[i][j][k].getRecNum();
							for (int kk=0;kk<rec_iter ;++kk){
						data[i][j][k].setWeight_Iter( data[i][j][k].getWeight_Iter(kk)*reweight_factor ,kk); 
						data[i][j][k].setWeight( data[i][j][k].getWeight(kk)*reweight_factor,kk); 
							// cout << i << "\t" << j << "\t" << k <<"\t" << data[i][j][k].getPosition().x <<"\t" << data[i][j][k].getPosition().y <<"\t" << data[i][j][k].getPosition().z << "\t" << data[i][j][k].getRecNum() << "\n" ;
						}	
					}
				}
				}

	}


	float get_the_iteration_difference() {
			int rec_iter;
			double log_diff;
			double max_log_diff=0.0;
			char pink;
			for(int i=0; i<(int)data.size(); i++){
					for(int j=0; j<(int)data[0].size(); j++){
						for(int k=0; k<(int)data[0][0].size(); k++){
							rec_iter = data[i][j][k].getRecNum();
							for (int kk=0;kk<rec_iter ;++kk){

								if(  abs(i-0.5*(int)data.size()) > 2 ){
								if(  abs(j-0.5*(int)data.size()) > 2 ){
								if(  abs(k-0.5*(int)data.size()) > 2 ){
								if(data[i][j][k].getWeight(kk) > 0.0 && data[i][j][k].getWeight_Iter(kk) >0.0)  { log_diff = abs(log10(data[i][j][k].getWeight(kk)/data[i][j][k].getWeight_Iter(kk)));}
								if(log_diff>max_log_diff){max_log_diff = log_diff;  cout << "Difference: " << i << " " << j << " " << k << " " << data[i][j][k].getIndex(kk) << " " << data[i][j][k].getWeight(kk) << "\t" << data[i][j][k].getWeight_Iter(kk) << std::endl;}
							}
						}
					}
								
								// cin >> pink;
							// cout << i << "\t" << j << "\t" << k <<"\t" << data[i][j][k].getPosition().x <<"\t" << data[i][j][k].getPosition().y <<"\t" << data[i][j][k].getPosition().z << "\t" << data[i][j][k].getRecNum() << "\n" ;
						}	
					}
				}
				}
	return max_log_diff;
	}



void read_grid_record(ifstream &file, float &weight, float &weight_iter, int &SFD_index, Vec3 &velocity){
	// char * memblock_id = new char [4]; // each int has 4 bytes
	char * memblock = new char [24]; // each float has 4 bytes
	float* read_values; 
	int* read_values_int; 

	// file.read (memblock_id, 4);
	// int* int_val = (int*)memblock_id;
	// parID = int_val[0]-1;

	file.read (memblock, 8);
	read_values = (float*)memblock;//reinterpret as doubles

	weight = read_values[0];
	weight_iter = read_values[1];


	file.read (memblock, 4);
	read_values_int = (int*)memblock;//reinterpret as doubles
	SFD_index = read_values_int[0];

	file.read (memblock, 12);
	read_values = (float*)memblock;//reinterpret as double

	velocity = Vec3(read_values[0],read_values[1],read_values[2]);
	delete memblock;
}

void save_grid_record(ofstream &file, float weight, float weight_iter, int SFD_index, Vec3 velocity){
    // file myFile ("test_grid3D.bin", ios::out | ios::binary);
    file.write ((char*)&weight, sizeof (weight));
    file.write ((char*)&weight_iter, sizeof (weight_iter));
    file.write ((char*)&SFD_index, sizeof (SFD_index));
    file.write ((char*)&velocity, sizeof (velocity));
}


	void save_grid3D(double initial_weight) {
		// FILE * pFile_temp;
		char pink;
		ofstream file("test_grid3D.bin", ios::out | ios::binary);
		   // ofstream wf("student.dat", ios::out | ios::binary);
   		if(!file) {
      		cout << "Function: save_grid3d - Cannot open file: test_grid3D.bin!" << endl;
      	return;
   	}
		// int rec_iter;
		// Vec3 position;

   		file.write((char*)&initial_weight, sizeof (initial_weight));
		for(int i=0; i<(int)data.size(); i++){
					for(int j=0; j<(int)data[0].size(); j++){
						for(int k=0; k<(int)data[0][0].size(); k++){
							// now we are cycling through every cell in the 3D Grid
							int rec_iter = data[i][j][k].getRecNum();
							Vec3 position = data[i][j][k].getPosition();
							file.write((char*)&rec_iter, sizeof (rec_iter));
							file.write((char*)&position, sizeof (position));

							for (int kk=0;kk<rec_iter ;++kk){
								float weight = data[i][j][k].getWeight(kk);
								float weight_iter = data[i][j][k].getWeight_Iter(kk);
								int SFD_index = data[i][j][k].getIndex(kk);
								Vec3 velocity = data[i][j][k].getVelocity(kk);
								save_grid_record(file, weight, weight_iter, SFD_index, velocity);
						// 		if(SFD_index>0){
						// 	save_grid_record(file, weight, weight_iter, SFD_index, velocity);
						// 	cout << "Saving_Sire: " << weight << "\t" << weight_iter << "\t" << SFD_index << "\t" << velocity.x << "\t" << velocity.y << "\t" << velocity.z << std::endl;
						// 	cin >> pink;
						// 	file.close();

						// 	ifstream file_in("test_grid3D.bin", ios::in | ios::binary);
						// 	read_grid_record(file_in, weight, weight_iter, SFD_index, velocity);
						// 	cout << "Saving_Sire Check: " << weight << "\t" << weight_iter << "\t" << SFD_index << "\t" << velocity.x << "\t" << velocity.y << "\t" << velocity.z << std::endl;
						// 	cin >> pink;
						// }
							}
						}
					}
				}
			file.close();
			}



		// getPosition()
		// getWeight(int i)
		// getWeight_Iter(int i)
		// getIndex(int i) 
		// getVelocity(int i)


	double load_grid3D() {
		char pink;
		int record_counter = 0;
		double initial_weight;
	ifstream file("test_grid3D.bin", ios::in | ios::binary);
		   // ofstream wf("student.dat", ios::out | ios::binary);
   		if(!file) {
      		cout << "Function: save_grid3d - Cannot open file: test_grid3D.bin!" << endl;
      	return 0.0;
   	}

   	file.read((char*) &initial_weight,sizeof (initial_weight));
		for(int i=0; i<(int)data.size(); i++){
			for(int j=0; j<(int)data[0].size(); j++){
				for(int k=0; k<(int)data[0][0].size(); k++){
					int rec_iter;
					Vec3 position;


					file.read((char*) &rec_iter,sizeof (rec_iter));
					file.read((char*) &position,sizeof (position));

				// if(rec_iter>0){
				// cout << "Loading_Sire: position: " << i << "\t" << j << "\t" << k << "\t" << position.x << "\t" << position.y << "\t" << position.z << "\t" << rec_iter << std::endl;
				// cin >> pink;
				// }

					data[i][j][k].setPosition(position);
					// data[i][j][k].setRecord(rec_iter);

					for (int kk=0;kk<rec_iter ;++kk){
						float weight;
						float weight_iter;
						int SFD_index;
						Vec3 velocity;
						read_grid_record(file, weight, weight_iter, SFD_index, velocity);
						// if(SFD_index>0){
						// cout << "Loading_Sire: " << weight << "\t" << weight_iter << "\t" << SFD_index << "\t" << velocity.x << "\t" << velocity.y << "\t" << velocity.z << std::endl;
						// cin >> pink;
						// }
						
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

