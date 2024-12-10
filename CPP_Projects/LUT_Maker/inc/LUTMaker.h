#ifndef LUT_MAKER_H
#define LUT_MAKER_H

#include <cmath>
#include <vector>

class LUTMaker
{
	public:

	LUTMaker();
	virtual ~LUTMaker();

	float GetX(const float &y, const float &phi_rad, const float &d) { return (y - (d-sin(phi_rad)*SXWIDTH/2+cos(phi_rad)*SYHEIGHT)/cos(phi_rad))/tan(phi_rad);}
	float GetY(const float &x, const float &phi_rad, const float &d) { return tan(phi_rad)*x + (d-sin(phi_rad)*SXWIDTH/2+cos(phi_rad)*SYHEIGHT/2)/cos(phi_rad);}

	void ComputeLengthMap();

	// Getters
	std::vector<std::vector<float>> GetLengthMap() { return arr_length; }
	std::vector<double> GetVecPhi() { return v_phi; }
	std::vector<double> GetVecD() { return v_d; }
	std::vector<double> GetVecZ() { return v_z; }
	float GetDiag() { return fdiag; }

	private:

	// Pad dimensions
	static constexpr float SXWIDTH =	11.28; // mm
	static constexpr float SYHEIGHT =	10.19; // mm
	const float fdiag = 				std::sqrt(std::pow(SXWIDTH, 2) + std::pow(SYHEIGHT, 2));

	// Number of discrete steps in each dimension of the Look Up Table
	static const int SNSTEPS_PHI      = 250;
	static const int SNSTEPS_D        = 250;
	static const int SNSTEPS_Z        = 101;
	static const int SNSTEPS_RC       = 2;
	static const int SNSTEPS_TRANS    = 2;

	std::vector<std::vector<float>> arr_length = std::vector<std::vector<float>>(SNSTEPS_PHI, std::vector<float>(SNSTEPS_D, 0));
	std::vector<double> v_phi, v_d, v_z;

};

#endif