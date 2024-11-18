#include "Includes.h"
#include "Util.h"
#include "SetStyle.h"

int main(){

	// Get parameter files
	std::string param_file =	"params.dat";
    std::ifstream infile(param_file);
	if (!infile) {
		std::cerr << "Unable to open parameter file: " << param_file << std::endl;
		return 1;
	}
    std::string line;

	// Create vectors of input files info
	std::vector<std::string>	v_files;
	std::vector<TFile*>			v_tfiles;
	std::vector<TTree*>			v_trees;
	std::vector<int>			v_nentries;
    
	// Get input files
    while (std::getline(infile, line)) {
		if (line.empty() || line[0] == '#') {
			continue;
		}
        if (line.rfind("FILE: ", 0) == 0) {
            v_files.			push_back(line.substr(5));
        }
    }

	// Open input files
	for (auto file : v_files) {
		TFile *fpInputFile =	TFile::Open(file.c_str());
		TTree *fpInputTree =	(TTree*)fpInputFile->Get("outTree");
		v_tfiles.				push_back(fpInputFile);
		v_trees.				push_back(fpInputTree);
		v_nentries.				push_back(fpInputTree->GetEntries());
	}

	// Set branches variables
	int hat, ep, eram, eram_id, col, row, cathode;
	double resolution, mom, drift, x, y, z;

	// Set branches addresses
	for (auto tree : v_trees) {
		tree->SetBranchAddress("hat", &hat);
		tree->SetBranchAddress("ep", &ep);
		tree->SetBranchAddress("eram", &eram);
		tree->SetBranchAddress("eram_id", &eram_id);
		tree->SetBranchAddress("col", &col);
		tree->SetBranchAddress("row", &row);
		tree->SetBranchAddress("cathode", &cathode);
		tree->SetBranchAddress("resolution", &resolution);
		tree->SetBranchAddress("mom", &mom);
		tree->SetBranchAddress("drift", &drift);
		tree->SetBranchAddress("x", &x);
		tree->SetBranchAddress("y", &y);
		tree->SetBranchAddress("z", &z);
	}

	// Set style (must be done before creating histograms)
	TStyle *style =				SetMyStyle();
	gROOT->						SetStyle(style->GetName());
	TCanvas *pcanvas =			new TCanvas("pcanvas", "pcanvas", 1800, 1350);
	gPad->						UseCurrentStyle();
	std::string output =		"output_pdf/SpatialResolution_MC_vs_Data.pdf";
	std::vector<int> colors =	{kRed+1, kBlue+1, kGreen+1};
	std::vector<std::string> legEntries = {"  DOG1 1022", "  Old sim", "  New sim"};
	std::vector<std::string> legdiffentries = {" Old - data", " New - data"};
	float norm, normdiff;

	// Set histograms
	std::vector<TH1F> 			v_h1f_botreso;
	std::vector<TH1F> 			v_h1f_topreso;
	std::vector<TH1F> 			v_h1f_botresodiffold;
	std::vector<TH1F> 			v_h1f_botresodiffnew;
	std::vector<TH1F> 			v_h1f_topresodiffold;
	std::vector<TH1F> 			v_h1f_topresodiffnew;
	for(int i=0; i<(int)v_files.size(); i++){
		v_h1f_botreso.			push_back(TH1F(Form("h1f_botreso_%d", i), "; Resolution (mm); Count (normalized)", 100, -4, 4));
		v_h1f_topreso.			push_back(TH1F(Form("h1f_topreso_%d", i), "; Resolution (mm); Count (normalized)", 100, -4, 4));
		v_h1f_botresodiffold.	push_back(TH1F(Form("h1f_botresodiffold_%d", i), "; Resolution (mm); Count (normalized)", 100, -4, 4));
		v_h1f_botresodiffnew.	push_back(TH1F(Form("h1f_botresodiffnew_%d", i), "; Resolution (mm); Count (normalized)", 100, -4, 4));
		v_h1f_topresodiffold.	push_back(TH1F(Form("h1f_topresodiffold_%d", i), "; Resolution (mm); Count (normalized)", 100, -4, 4));
		v_h1f_topresodiffnew.	push_back(TH1F(Form("h1f_topresodiffnew_%d", i), "; Resolution (mm); Count (normalized)", 100, -4, 4));
	}

	// Fill histograms
	for(int i=0; i<(int)v_files.size(); i++){
		std::cout << "Processing file: " << v_files[i] << std::endl;
		for(int j=0; j<v_nentries[i]; j++){
			v_trees[i]->		GetEntry(j);
			if(x < 700) continue;
			if(hat==0) v_h1f_botreso[i].	Fill(resolution);
			if(hat==1) v_h1f_topreso[i].	Fill(resolution);
		}
	}

	// Normalization of the histograms
	for(int i=0; i<(int)v_files.size(); i++){
		v_h1f_botreso[i].		Scale(1/v_h1f_botreso[i].Integral());
		v_h1f_topreso[i].		Scale(1/v_h1f_topreso[i].Integral());
	}

	// Histogram bin to bin difference
	for(int i=1; i<(int)v_files.size(); i++){
		v_h1f_botresodiffold[i].	Add(&v_h1f_botreso[i], &v_h1f_botreso[0], 1, -1);
		v_h1f_botresodiffnew[i].	Add(&v_h1f_botreso[i], &v_h1f_botreso[0], 1, -1);
		v_h1f_topresodiffold[i].	Add(&v_h1f_topreso[i], &v_h1f_topreso[0], 1, -1);
		v_h1f_topresodiffnew[i].	Add(&v_h1f_topreso[i], &v_h1f_topreso[0], 1, -1);
	}

	// Set legend
	TLegend *legend =			new TLegend(0.65,0.7,0.85,0.91);
	legend->					SetTextSize(0.06);
	legend->					SetFillStyle(0);
	legend->					SetTextColor(kBlue-1);
	legend->					SetNColumns(2);
	legend->					SetHeader("tHAT bHAT");

	TLegend *legdiff =			new TLegend(0.75,0.65,0.9,0.91);
	legdiff->					SetTextSize(0.1);
	legdiff->					SetFillStyle(0);
	legdiff->					SetTextColor(kBlue-1);
	legdiff->					SetNColumns(2);
	legdiff->					SetHeader("tHAT  bHAT");

	// Draw settings
	pcanvas->					cd();
	gStyle->					SetOptStat(0);
	gStyle->					SetOptFit(0);

	// Draw histogram
	TPad *pad1 =				new TPad("pad1", "pad1", 0, 0.3, 1, 1, 0, 0, 0);
	pad1->						SetBottomMargin(0);
	pad1->						SetTopMargin(0.03);
	pad1->						Draw();
	TPad *pad2 =				new TPad("pad2", "pad2", 0, 0, 1, 0.3, 0, 0, 0);
	pad2->						SetTopMargin(0);
	pad2->						SetBottomMargin(0.3);
	pad2->						Draw();
	norm = 0, normdiff = 0;
	for(int i=0; i<(int)v_files.size(); i++){
		pad1->					cd();
		v_h1f_botreso[i].GetMaximum() > norm ? norm = v_h1f_botreso[i].GetMaximum() : norm;
		v_h1f_topreso[i].GetMaximum() > norm ? norm = v_h1f_topreso[i].GetMaximum() : norm;
		Graphic_setup(&v_h1f_topreso[i], 0.5, 1, colors[i]-1, 2, colors[i]-1, colors[i]-1, 0.1);
		Graphic_setup(&v_h1f_botreso[i], 0.5, 1, colors[i]+1, 2, colors[i]+1, colors[i]+1, 0.1);
		v_h1f_topreso[0].		SetMaximum(norm*1.1);
		v_h1f_topreso[0].		GetXaxis()->SetLabelSize(0);
		v_h1f_topreso[0].		GetXaxis()->SetTitle("");
		v_h1f_topreso[i].		Draw(i==0 ? "HIST" : "HIST same");
		v_h1f_botreso[i].		Draw("HIST same");
		legend->				AddEntry(&v_h1f_topreso[i], "         ", "f");
		legend->				AddEntry(&v_h1f_botreso[i], legEntries[i].c_str(), "f");

		if(i==0) continue;
		pad2->					cd();
		v_h1f_botresodiffold[i].GetMaximum() > normdiff ? normdiff = v_h1f_botresodiffold[i].GetMaximum() : normdiff;
		v_h1f_botresodiffnew[i].GetMaximum() > normdiff ? normdiff = v_h1f_botresodiffnew[i].GetMaximum() : normdiff;
		v_h1f_topresodiffold[i].GetMaximum() > normdiff ? normdiff = v_h1f_topresodiffold[i].GetMaximum() : normdiff;
		v_h1f_topresodiffnew[i].GetMaximum() > normdiff ? normdiff = v_h1f_topresodiffnew[i].GetMaximum() : normdiff;
		v_h1f_botresodiffold[1].SetMaximum(normdiff*1.1);
		v_h1f_botresodiffold[1].GetXaxis()->SetTitleSize(0.12);
		v_h1f_botresodiffold[1].GetXaxis()->SetLabelSize(0.1);
		v_h1f_botresodiffold[1].GetYaxis()->SetTitleSize(0.12);
		v_h1f_botresodiffold[1].GetYaxis()->SetLabelSize(0.1);
		Graphic_setup(&v_h1f_botresodiffold[i], 0.5, 1, colors[i]+1, 1, colors[i]+1, 0, 0);
		Graphic_setup(&v_h1f_botresodiffnew[i], 0.5, 1, colors[i]+1, 1, colors[i]+1, 0, 0);
		Graphic_setup(&v_h1f_topresodiffold[i], 0.5, 1, colors[i]-1, 1, colors[i]-1, 0, 0);
		Graphic_setup(&v_h1f_topresodiffnew[i], 0.5, 1, colors[i]-1, 1, colors[i]-1, 0, 0);
		v_h1f_botresodiffold[i].Draw(i==0 ? "HIST" : "HIST same");
		v_h1f_botresodiffnew[i].Draw("HIST same");
		v_h1f_topresodiffold[i].Draw("HIST same");
		v_h1f_topresodiffnew[i].Draw("HIST same");
		legdiff->				AddEntry(&v_h1f_topresodiffold[i], "          ", "f");
		legdiff->				AddEntry(&v_h1f_botresodiffnew[i], legdiffentries[i-1].c_str(), "f");
	}
	legdiff->					Draw();
	TLine *baseline = new TLine(-4, 0, 4, 0);
	baseline->SetLineColor(kGray+2);
	baseline->SetLineStyle(2);
	baseline->Draw();
	pad1->						cd();
	legend->					Draw();
	pcanvas->					SaveAs(output.c_str());

	// Close input files
	for (auto file : v_tfiles) {
		file->					Close();
	}

	// Deleting pointers
	delete style;
	delete pcanvas;

	return 0;
}