#include "LUTMaker.h"

#include <chrono>

#include "Misc_Functions.h"

// Constructor
LUTMaker::LUTMaker(){
	v_phi = linspace(0, 90, SNSTEPS_PHI);
	v_d = linspace(0, fdiag/2, SNSTEPS_D);
	v_z = linspace(0, 1000, SNSTEPS_Z);
}

// Destructor
LUTMaker::~LUTMaker(){
}


// Compute the length in the pad for all (phi, d) pairs
void LUTMaker::ComputeLengthMap(){

    auto start = std::chrono::high_resolution_clock::now();
	double time = 0;

	for (int i = 0; i < SNSTEPS_PHI; i++){
		float phi_rad = v_phi[i]/180*M_PI;
		if(v_phi[i] == 0) phi_rad = 1e-5/180*M_PI;
		if(v_phi[i] == 90) phi_rad = (90-1e-5)/180*M_PI;

		for (int j = 0; j < SNSTEPS_D; j++){
			float d = v_d[j];
			std::vector<float> v_x, v_y;

			// Compute intersections with lines defining the pad borders
			float y_xmin = GetY(0, 		  phi_rad, d);
			float y_xmax = GetY(SXWIDTH,  phi_rad, d);
			float x_ymin = GetX(0, 		  phi_rad, d);
			float x_ymax = GetX(SYHEIGHT, phi_rad, d);

			if(0 <= y_xmin and y_xmin < SYHEIGHT){
				v_x.push_back(0);
				v_y.push_back(y_xmin);
			}
			if(0 <= y_xmax and y_xmax < SYHEIGHT){
				v_x.push_back(SXWIDTH);
				v_y.push_back(y_xmax);
			}
			if(0 <= x_ymin and x_ymin < SXWIDTH){
				v_x.push_back(x_ymin);
				v_y.push_back(0);
			}
			if(0 <= x_ymax and x_ymax < SXWIDTH){
				v_x.push_back(x_ymax);
				v_y.push_back(x_ymax);
			}

			if(v_x.size() != 2 or v_y.size() != 2) continue;
			arr_length[i][j] = std::sqrt(std::pow(v_x[1]-v_x[0], 2) + std::pow(v_y[1]-v_y[0], 2));	
			std::cout << v_x[0] << " " << v_y[0] << " " << v_x[1] << " " << v_y[1] << " " << v_phi[i] << " " << v_d[j] << " " << arr_length[i][j] << std::endl;
			v_x.clear();
			v_y.clear();
		}
	}

    // Calculate elapsed time in milliseconds
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    time += duration.count();
    std::cout << "Length map done in: " << duration.count() << " ms" << std::endl;
}