#include <iomanip>

#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "Util.h"
#include "Misc_Functions.h"
#include "SetStyle.h"

#include "TrackModel.h"
#include "LUTMaker.h"

int main()
{

	LUTMaker *p_lutmaker = new LUTMaker();
	p_lutmaker->ComputeLengthMap();
	// p_lutmaker->DrawLengthMap();

	p_lutmaker->MakeLUT();

	// p_lutmaker->LoadLUT("Output_LUT/LUT_test.root");
	// // p_lutmaker->DrawLUT(120, 0, 310/pow(10, 	3.5));

	// float diag = p_lutmaker->GetDiag();
	// p_lutmaker->DrawDiffInterpolExact(0, 90, 0, diag/2);

	// TrackModel trackmodel;
	// trackmodel.ComputeRealCharge(671, 10.3876, 0.134317, 78.8042, 143.559, 640.9, 310/pow(10, 3.5));
	// std::cout << "Amplitude: " << 311 << std::endl;
	// std::cout << "Length: " << 8.25727 << std::endl;
	// std::cout << "Impact: " << 5.02252 << std::endl;
	// std::cout << "Phi: " << 79.9085 << std::endl;
	// std::cout << "RC: " << 121.395 << std::endl;
	// std::cout << "Drift: " << 671.9 << std::endl;
	// std::cout << "Real charge: " << trackmodel.GetRealCharge() << std::endl;

	delete p_lutmaker;

	return 0;
}
// 726.653 671, 10.3876, 0.134317, 78.8042, 143.559, 640.9,
// 736.87 - 726.653 = 10.2166(relative: 1.40598%)