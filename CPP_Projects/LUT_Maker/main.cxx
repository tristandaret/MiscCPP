#include <iomanip>

#include "Util.h"

#include "Misc_Functions.h"

#include "TrackModel.h"
#include "LUTMaker.h"

int main()
{
	LUTMaker *p_lutmaker = new LUTMaker();
	p_lutmaker->ComputeLengthMap();

	std::vector<std::vector<float>> arr_length = p_lutmaker->GetLengthMap();
	std::vector<double> v_phi = p_lutmaker->GetVecPhi();
	std::vector<double> v_d = p_lutmaker->GetVecD();
	std::vector<double> v_z = p_lutmaker->GetVecZ();
	float diag = p_lutmaker->GetDiag();

	TH2F h2("h2", "Length map", v_phi.size(), 0, 90, v_d.size(), 0, diag/2);

	for (int i = 0; i < (int)v_phi.size(); ++i) {
		for (int j = 0; j < (int)v_d.size(); ++j) {
			h2.SetBinContent(i+1, j+1, arr_length[i][j]);
		}
	}

	TCanvas c1("c1", "Canvas", 1800, 1350);
	gStyle->SetOptStat(0);
	h2.Draw("COLZ");
	c1.SaveAs("length_map.pdf");

	delete p_lutmaker;
	
	return 0;
}