#include "LUTMaker.h"

#include <chrono>

#include "Misc_Functions.h"
#include "SetStyle.h"

#include "TrackModel.h"

// Constructor
LUTMaker::LUTMaker(){
	v_phi =	linspace(0, 90, SNSTEPS_PHI);
	v_d =	linspace(0, fdiag/2, SNSTEPS_D);
	v_z =	linspace(0, 1000, SNSTEPS_Z);
	v_RC =	{112, 158};
	v_Dt =	{310/pow(10, 3.5), 350/pow(10, 3.5)};
}

// Destructor
LUTMaker::~LUTMaker(){
}


// Compute the length in the pad for all (phi, d) pairs
void LUTMaker::ComputeLengthMap(){

    auto start = std::chrono::high_resolution_clock::now();
	double time = 0;

	for (int i = 0; i < SNSTEPS_PHI; i++){
		float phi = v_phi[i];
		float phi_rad = phi/180*M_PI;
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
				v_y.push_back(SYHEIGHT);
			}
			if(v_x.size() == 2 and v_y.size() == 2){
				arr_length[i][j] = std::sqrt(std::pow(v_x[1]-v_x[0], 2) + std::pow(v_y[1]-v_y[0], 2));	
			}
		}
	}

    // Calculate elapsed time in milliseconds
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    time += duration.count();
    std::cout << "Length map done in: " << duration.count() << " ms" << std::endl;
}

// Draw length map
void LUTMaker::DrawLengthMap(){

	TCanvas c1("c1", "c1", 1800, 1350);
	TStyle* ptstyle =				SetMyStyle();
	gROOT->							SetStyle(ptstyle->GetName());
	gStyle->						SetPadLeftMargin(0.04);
	gStyle->						SetPadRightMargin(0.11);
	gStyle->						SetPadBottomMargin(0.06);
	gStyle->						SetPadTopMargin(0.08);
	gPad->							UseCurrentStyle();
	gStyle->SetOptStat(0);

	TH2F h2("h2", "Length map", v_phi.size(), 0, 90, v_d.size(), 0, fdiag/2);
	for (int i = 0; i < (int)v_phi.size(); ++i) {
		for (int j = 0; j < (int)v_d.size(); ++j) {
			h2.SetBinContent(i+1, j+1, arr_length[i][j]);
		}
	}
	h2.Draw("COLZ");
	c1.SaveAs("Output_PDF/length_map.pdf");

	delete ptstyle;
}
// Compute scale factor
void LUTMaker::ComputeScaleFactor(){

	TFile tfileLUT("Output_LUT/LUT_test.root", "RECREATE");
	TTree ttreeLUT("LUT", "LUT");
	float length, d, phi, RC, Dt, drift, scalefactor;
	ttreeLUT.Branch("length", &length);
	ttreeLUT.Branch("d", &d);
	ttreeLUT.Branch("phi", &phi);
	ttreeLUT.Branch("RC", &RC);
	ttreeLUT.Branch("Dt", &Dt);
	ttreeLUT.Branch("drift", &drift);
	ttreeLUT.Branch("scalefactor", &scalefactor);

    auto start = std::chrono::high_resolution_clock::now();
	double time = 0;

	TrackModel trackmodel;

	Dt = v_Dt[0];
	RC = v_RC[0];
	drift = v_z[0];

	for(int iphi = 0; iphi < SNSTEPS_PHI; iphi++){
		phi = v_phi[iphi];
		std::cout << "Computing phi: " << phi << "°" << std::endl;
		for(int id = 0; id < SNSTEPS_D; id++){
			d = v_d[id];
			length = arr_length[iphi][id];
			trackmodel.ComputeAmplitudeLoss(length, v_d[id], v_phi[iphi], RC, drift, Dt);
			scalefactor = 1/trackmodel.GetAmplitudeLoss();
			if (std::isnan(scalefactor) || std::isinf(scalefactor) || scalefactor <= 0) continue;
			ttreeLUT.Fill();
		}
	}

	tfileLUT.Write();
	tfileLUT.Close();						

    // Calculate elapsed time in milliseconds
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    time += duration.count();
    std::cout << "LUT done in: " << duration.count() << " s" << std::endl;
}


// Load LUT
void LUTMaker::LoadLUT(std::string LUTpath){
	// Initialize the LUT with zeros
	for (int i = 0; i < SNSTEPS_PHI; i++){
		for (int j = 0; j < SNSTEPS_D; j++){
			LUTValues[i][j] = 0;
		}
	}

	// Load the LUT from the file
	TFile tfile(LUTpath.c_str(), "READ");
	TTree &ptree = *(TTree*)tfile.Get("LUT");
	float length, d, phi, RC, Dt, drift, scalefactor;
	ptree.SetBranchAddress("length", &length);
	ptree.SetBranchAddress("d", &d);
	ptree.SetBranchAddress("phi", &phi);
	ptree.SetBranchAddress("RC", &RC);
	ptree.SetBranchAddress("Dt", &Dt);
	ptree.SetBranchAddress("drift", &drift);
	ptree.SetBranchAddress("scalefactor", &scalefactor);


	TCanvas c1("c1", "c1", 1800, 1350);
	TStyle* ptstyle =				SetMyStyle();
	gROOT->							SetStyle(ptstyle->GetName());
	gStyle->						SetPadLeftMargin(0.04);
	gStyle->						SetPadRightMargin(0.11);
	gStyle->						SetPadBottomMargin(0.06);
	gStyle->						SetPadTopMargin(0.08);
	gPad->							UseCurrentStyle();
	gStyle->SetOptStat(0);
	TH2F h2_scalefactor("h2_scalefactor", "Scale factor", SNSTEPS_PHI, 0, 90.01, SNSTEPS_D, 0, fdiag/2);

	for (int i = 0; i < ptree.GetEntries(); i++) {
		ptree.GetEntry(i);
		int iphi = (int)std::round(phi/sSTEP_PHI);
		int id = (int)std::round(d/sSTEP_D);
		LUTValues[iphi][id] = scalefactor;
			h2_scalefactor.SetBinContent(h2_scalefactor.FindBin(phi, d), scalefactor);
	}
	h2_scalefactor.Draw("COLZ");
	c1.SaveAs("Output_PDF/scalefactor_map1.pdf");
	tfile.Close();
}


// Get factor from LUT
float LUTMaker::GetFactorFromLUT(const float &phi, const float &d){
	float iphi =			phi/sSTEP_PHI;
	float iphi_min =		std::min(std::floor(phi/sSTEP_PHI), (float)SNSTEPS_PHI-1);
	float iphi_max =		std::max(std::ceil(phi/sSTEP_PHI), (float)0);
	float id =				d/sSTEP_D;
	float id_min =			std::min(std::floor(d/sSTEP_D), (float)SNSTEPS_D-1);
	float id_max =			std::max(std::ceil(d/sSTEP_D), (float)0);

	// weights
	float w_phi, w_d;
	if(iphi_min == iphi_max)	w_phi = 1;
	else						w_phi = 1 - (iphi - iphi_min)/(iphi_max - iphi_min);
	if(id_min == id_max)		w_d = 1;
	else						w_d = 1 - (id - id_min)/(id_max - id_min);

	// Interpolation
	float factor = 0;
	factor += w_phi     * w_d     * LUTValues[(int)iphi_min][(int)id_min];
	factor += w_phi     * (1-w_d) * LUTValues[(int)iphi_min][(int)id_max];
	factor += (1-w_phi) * w_d     * LUTValues[(int)iphi_max][(int)id_min];
	factor += (1-w_phi) * (1-w_d) * LUTValues[(int)iphi_max][(int)id_max];

	std::cout << phi << "(" << iphi << ", " << iphi_min << ", " << iphi_max << ") " << w_phi << " " << d << "(" << id << ", " << id_min << ", " << id_max << ") " << w_d << " " << factor << std::endl;

	return factor;
}


// Draw LUT
void LUTMaker::DrawLUT(const float &RC, const float &drift, const float &Dt){

	TCanvas c1("c1", "c1", 1800, 1350);
	TStyle* ptstyle =				SetMyStyle();
	gROOT->							SetStyle(ptstyle->GetName());
	gStyle->						SetPadLeftMargin(0.04);
	gStyle->						SetPadRightMargin(0.11);
	gStyle->						SetPadBottomMargin(0.06);
	gStyle->						SetPadTopMargin(0.08);
	gPad->							UseCurrentStyle();
	gStyle->SetOptStat(0);

	TH2F h2_scalefactor("h2_scalefactor", "Scale factor", SNSTEPS_PHI, 0, 90.01, SNSTEPS_D, 0, fdiag/2);
	for (int i = 0; i < SNSTEPS_PHI; i++) {
		for (int j = 0; j < SNSTEPS_D; j++) {
			float phi = v_phi[i];
			float d = v_d[j];
			float scalefactor = GetFactorFromLUT(phi, d);
			h2_scalefactor.SetBinContent(h2_scalefactor.FindBin(phi, d), scalefactor);
		}
	}
	h2_scalefactor.Draw("COLZ");
	c1.SaveAs("Output_PDF/scalefactor_map2.pdf");

	delete ptstyle;
}