#include "GeomFunctions.h"

#include <vector>

float ComputeLength(const float &phi, const float &d){
	float phi_rad = phi/180*M_PI;
	std::vector<float> v_x, v_y;

	// Compute intersections with lines defining the pad borders
	float y_xmin = GetY(0, 		  phi_rad, d);
	float y_xmax = GetY(11.28,  phi_rad, d);
	float x_ymin = GetX(0, 		  phi_rad, d);
	float x_ymax = GetX(10.19, phi_rad, d);
	if(0 <= y_xmin and y_xmin < 10.19){
		v_x.push_back(0);
		v_y.push_back(y_xmin);
	}
	if(0 <= y_xmax and y_xmax < 10.19){
		v_x.push_back(11.28);
		v_y.push_back(y_xmax);
	}
	if(0 <= x_ymin and x_ymin < 11.28){
		v_x.push_back(x_ymin);
		v_y.push_back(0);
	}
	if(0 <= x_ymax and x_ymax < 11.28){
		v_x.push_back(x_ymax);
		v_y.push_back(10.19);
	}
	return std::sqrt(std::pow(v_x[1]-v_x[0], 2) + std::pow(v_y[1]-v_y[0], 2));	
}