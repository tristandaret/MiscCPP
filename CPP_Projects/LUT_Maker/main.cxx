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
	// p_lutmaker->ComputeLengthMap();
	// p_lutmaker->DrawLengthMap();

	// p_lutmaker->ComputeScaleFactor();

	p_lutmaker->LoadLUT("Output_LUT/LUT_test.root");
	p_lutmaker->DrawLUT(120, 0, 310/pow(10, 3.5));

	delete p_lutmaker;

	return 0;
}