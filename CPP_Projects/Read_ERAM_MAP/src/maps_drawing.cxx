#include "maps_drawing.h"
#include "LUTs.h"
#include "Misc_Functions.h"
#include "SetStyle.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TROOT.h>

Drawing::Drawing(){
	pERAMMaps = new Reconstruction::ERAMMaps();

	TStyle* ptstyle = SetMyStyle();
	gROOT->SetStyle(ptstyle->GetName());
    pTCanvas           = new TCanvas("pTCanvas", "pTCanvas", 1800, 1400) ;
}



Drawing::~Drawing(){
	delete pERAMMaps;
	delete pTCanvas;
}




void Drawing::maps_pdf(){
	std::cout << "Running maps_pdf" << std::endl;

  	std::vector<TH2F*> v_h2_gain;
  	std::vector<TH2F*> v_h2_rc;
  	for(int i=0;i<34;i++){
		v_h2_gain.push_back(new TH2F(Form("h2_gain_%d",i),"",36, -0.5, 35.5, 32, -0.5, 31.5));
		v_h2_rc.push_back(new TH2F(Form("h2_rc_%d",i),"",36, -0.5, 35.5, 32, -0.5, 31.5));
  	}
  	std::vector<float> v_low_rc;
  	std::vector<float> v_low_gain;
  	float low_rc    = 1e6;
  	float low_gain  = 1e6;
  	std::vector<float> v_high_rc;
  	std::vector<float> v_high_gain;
  	float high_rc     = 0;
  	float high_gain   = 0;
  	std::vector<float> v_avg_rc;
  	std::vector<float> v_avg_gain;
  	float avg_rc      = 0;
  	float avg_gain    = 0;
  	float n_pads      = 0;
  	for(int i=0;i<34;i++){
		float n_pads_eram = 0;
		float eram_gain = 0;
		float eram_rc   = 0;
		low_rc          = 1e6;
		low_gain        = 1e6;
		high_rc         = 0;
		high_gain       = 0;
		int id          = pERAMMaps->ID(i);
		for(int j=0;j<36;j++){
			for(int k=0;k<32;k++){
			float gain  = pERAMMaps->Gain(i,j,k);
			float rc    = pERAMMaps->RC(i,j,k);
			if(gain != 0) v_h2_gain[i]->Fill(j,k,gain);
			if(rc != 0)   v_h2_rc[i]->Fill(j,k,rc);
			if(gain == 0 or rc == 0) std::cout << "ERAM#" << id << " (" << j << "," << k << ") has a hole" << std::endl;
			// lowest values
			if(rc < low_rc and rc != 0) low_rc = rc;
			if(gain < low_gain and gain != 0) low_gain = gain;
			// highest values
			if(rc > high_rc) high_rc = rc;
			if(gain > high_gain) high_gain = gain;
			// average values
			avg_gain += gain;
			eram_gain += gain;
			eram_rc += rc;
			avg_rc += rc;
			n_pads++;
			n_pads_eram++;
			}
		}
		v_h2_gain[i]->    	SetTitle(Form("ERAM#%i (%i) Gain = %0.f",id, i, eram_gain/n_pads_eram));
		v_h2_rc[i]->      	SetTitle(Form("ERAM#%i (%i) RC = %0.f ns/mm^{2}",id, i, eram_rc/n_pads_eram));
		v_low_rc.			push_back(low_rc);
		v_low_gain.			push_back(low_gain);
		v_high_rc.			push_back(high_rc);
		v_high_gain.		push_back(high_gain);
		v_avg_rc.			push_back(eram_rc/n_pads_eram);
		v_avg_gain.			push_back(eram_gain/n_pads_eram);
  	}
  	low_rc = 				*std::min_element(v_low_rc.begin(), v_low_rc.end());
  	low_gain = 				*std::min_element(v_low_gain.begin(), v_low_gain.end());
  	high_rc = 				*std::max_element(v_high_rc.begin(), v_high_rc.end());
  	high_gain = 			*std::max_element(v_high_gain.begin(), v_high_gain.end());
  	avg_rc /= 				n_pads;
  	avg_gain /= 			n_pads;
  	// std::cout << "Average RC: " << avg_rc << std::endl;
  	// std::cout << "Average Gain: " << avg_gain << std::endl;

  	// for(int i=0;i<(int)v_avg_gain.size();i++) std::cout <<"Gain ERAM#" << std::setw(2) << pERAMMaps->ID(i) << "(" << std::setw(2) << i << "): " << std::setprecision(4) << v_avg_gain[i] << std::endl;
  	// for(int i=0;i<(int)v_avg_rc.size();i++)   std::cout <<"RC   ERAM#" << std::setw(2) << pERAMMaps->ID(i) << "(" << std::setw(2) << i << "): " << std::setprecision(3) << v_avg_rc[i] << std::endl;
	
	std::vector<int> CERN22 = {26, 12, 10, 7, 17, 19, 13, 32, 33};

  	/* Display */
  	OutputFile      = "Output/maps_filled_test.pdf";
  	OutputFile_Beg  = OutputFile + "(" ;
  	OutputFile_End  = OutputFile + ")" ;
	gPad->						UseCurrentStyle();
  	gStyle->                    SetOptStat(0);
  	gStyle->                    SetPalette(kRainBow);

	pTCanvas->					Divide(4,4);
	for(int i=0;i<16;i++){
		pTCanvas->cd(i+1);
		v_h2_gain[i]->			SetMinimum(low_gain);
		v_h2_gain[i]->			SetMaximum(high_gain);
		v_h2_gain[i]->			Draw("colz");
  	}
  	pTCanvas->					SaveAs(OutputFile_Beg.c_str());
  	pTCanvas->					Clear();

	pTCanvas->					Divide(4,4);
  	for(int i=16;i<32;i++){
		pTCanvas->cd(i-15);
		v_h2_gain[i]->			SetMinimum(low_gain);
		v_h2_gain[i]->          SetMaximum(high_gain);
		v_h2_gain[i]->          Draw("colz");
  	}
  	pTCanvas->					SaveAs(OutputFile.c_str());
  	pTCanvas->					Clear();

	pTCanvas->					Divide(4,4);
	for(int i=0;i<9;i++){
		pTCanvas->				cd(i+1);
		v_h2_gain[CERN22[i]]->	SetMinimum(low_gain);
		v_h2_gain[CERN22[i]]->	SetMaximum(high_gain);
		v_h2_gain[CERN22[i]]->	Draw("colz");
	}
	pTCanvas->					SaveAs(OutputFile.c_str());
	pTCanvas->					Clear();

	pTCanvas->					Divide(4,4);
  	for(int i=0;i<16;i++){
  	  	pTCanvas->				cd(i+1);
		v_h2_rc[i]->            SetMinimum(low_rc);
		v_h2_rc[i]->            SetMaximum(high_rc);
		v_h2_rc[i]->            Draw("colz");
  	}
  	pTCanvas->					SaveAs(OutputFile.c_str());
  	pTCanvas->					Clear();

	pTCanvas->					Divide(4,4);
  	for(int i=16;i<32;i++){
		pTCanvas->				cd(i-15);
		v_h2_rc[i]->            SetMinimum(low_rc);
		v_h2_rc[i]->            SetMaximum(high_rc);
		v_h2_rc[i]->            Draw("colz");
  	}
  	pTCanvas->					SaveAs(OutputFile.c_str());
  	pTCanvas->					Clear();

	pTCanvas->					Divide(4,4);
	for(int i=0;i<9;i++){
		pTCanvas->cd(i+1);
		v_h2_rc[CERN22[i]]->SetMinimum(low_rc);
		v_h2_rc[CERN22[i]]->SetMaximum(high_rc);
		v_h2_rc[CERN22[i]]->Draw("colz");
	}
	pTCanvas->SaveAs(OutputFile.c_str());
	pTCanvas->Clear();

  	for(int i=0;i<34;i++){
  	  pTCanvas->cd();
  	  v_h2_gain[i]->SetMinimum(v_low_gain[i]);
  	  v_h2_gain[i]->SetMaximum(v_high_gain[i]);
  	  v_h2_gain[i]->Draw("colz");
  	  pTCanvas->SaveAs(OutputFile.c_str());
  	  pTCanvas->Clear();
  	}

  	for(int i=0;i<34;i++){
  	  pTCanvas->cd();
  	  v_h2_rc[i]->SetMinimum(v_low_rc[i]);
  	  v_h2_rc[i]->SetMaximum(v_high_rc[i]);
  	  v_h2_rc[i]->Draw("colz");
  	  if(i == 33) pTCanvas->SaveAs(OutputFile_End.c_str());
  	  else pTCanvas->SaveAs(OutputFile.c_str());
  	  pTCanvas->Clear();
  	}
}







void Drawing::maps_csv(){
	std::cout << "Running maps_csv" << std::endl;

  	std::vector<float> v_high_rc;
  	std::vector<float> v_avg_rc;
  	float high_rc = 		0;
  	float avg_rc =			0;
  	float n_pads = 			0;
  	for(int i=0;i<34;i++){
		float n_pads_eram = 0;
		float eram_rc = 	0;
		high_rc = 			0;
		for(int j=0;j<36;j++){
			for(int k=0;k<32;k++){
			float rc    = pERAMMaps->RC(i,j,k);
			if(rc > high_rc) high_rc = rc;
			eram_rc += rc;
			avg_rc += rc;
			n_pads++;
			n_pads_eram++;
			}
		}
		v_high_rc.			push_back(high_rc);
		v_avg_rc.			push_back(eram_rc/n_pads_eram);
  	}
  	high_rc = 				*std::max_element(v_high_rc.begin(), v_high_rc.end());
  	avg_rc /= 				n_pads;

	std::cout << "Average RC: " << avg_rc << std::endl;
	std::cout << "Highest RC: " << high_rc << std::endl;

	for(int i=0;i<34;i++){
		std::cout <<"RC   ERAM#" << std::setw(2) << pERAMMaps->ID(i) << "(" << std::setw(2) << i << "): " << std::setprecision(3) << v_avg_rc[i] << std::endl;
		std::ofstream csvFile(Form("Output/CSV/ERAM_%i_RC_grayscale.csv", pERAMMaps->ID(i)));
		for(int j=0;j<32;j++){
			for(int k=0;k<36;k++){
				float rc = pERAMMaps->RC(i,k,31-j);
				csvFile << rc/high_rc*255;
				if(k < 35) csvFile << ",";
			}
			csvFile << std::endl;
		}
		csvFile.close();
	}
}





void Drawing::Distributions(){
	int col = 36;
	int row = 32;
	int nerams = 32;
	std::cout << "Running Distributions" << std::endl;

	TH2F *ph2_GRC =			new TH2F("ph2_GRC", "Gain vs RC;RC value (ns/mm^{2});Gain value (ADC count)", 50, 50, 270, 50, 700, 2700);
	std::vector<TH2F*> v_ph2_GRC;
	std::vector<float> v_min_rc;
	std::vector<float> v_max_rc;
	for(int i=0;i<nerams;i++)v_ph2_GRC.push_back(new TH2F(Form("ph2_GRC_%d",i), Form("ERAM#%i Gain vs RC;RC value (ns/mm^{2});Gain value (ADC count)", pERAMMaps->ID(i)), 50, 50, 270, 50, 700, 2700));
	TH2F *ph2_gainmap = 	new TH2F("ph2_gainmap", "Average gain map;;;ns/mm^{2}", col, -0.5, 35.5, row, -0.5, 31.5);

	TH2F *ph2_rcmap = 		new TH2F("ph2_rcmap", "Average RC map;;;ns/mm^{2}", col, -0.5, 35.5, row, -0.5, 31.5);
	TH2F *ph2_rcmapbot = 	new TH2F("ph2_rcmapbot", "Average RC map (bHAT);;;ns/mm^{2}", col, -0.5, 35.5, row, -0.5, 31.5);
	TH2F *ph2_rcmaptop = 	new TH2F("ph2_rcmaptop", "Average RC map (tHAT);;;ns/mm^{2}", col, -0.5, 35.5, row, -0.5, 31.5);
	TH2F *ph2_rcdiffTPC =	new TH2F("ph2_rcdiffTPC", "Difference of tHAT vs bHAT average RC map;;;ns/mm^{2}", col, -0.5, 35.5, row, -0.5, 31.5);

	TH2F *ph2_rcnorm = 		new TH2F("ph2_rcnorm", "Normalized average RC map", col, -0.5, 35.5, row, -0.5, 31.5);
	TH2F *ph2_rcnormbot = 	new TH2F("ph2_rcnormbot", "Normalized average RC map (bHAT)", col, -0.5, 35.5, row, -0.5, 31.5);
	TH2F *ph2_rcnormtop = 	new TH2F("ph2_rcnormtop", "Normalized average RC map (tHAT)", col, -0.5, 35.5, row, -0.5, 31.5);
	TH2F *ph2_rcnormdiffTPC = new TH2F("ph2_rcnormdiffTPC", "Absolute difference of tHAT vs bHAT normalized average RC map (in %)", col, -0.5, 35.5, row, -0.5, 31.5);
	
	float csv_rcnorm[col][row];
	float csv_rcnormbot[col][row];
	float csv_rcnormtop[col][row];
	for(int j=0;j<col;j++) for(int k=0;k<row;k++){
		csv_rcnorm[j][k] = 0;
		csv_rcnormbot[j][k] = 0;
		csv_rcnormtop[j][k] = 0;
	}
    std::ofstream CSVRCNorm("Output/CSV_RC_normed/RC_normed.csv");
	std::ofstream CSVRCNormBot("Output/CSV_RC_normed/RC_normed_bottom.csv");
	std::ofstream CSVRCNormTop("Output/CSV_RC_normed/RC_normed_top.csv");


	TH2F *ph2_rcnormstd = 	new TH2F("ph2_rcnormstd", "Std of normalized average RC difference map (in %)", col, -0.5, 35.5, row, -0.5, 31.5);
	TH1F *ph1_gain =		new TH1F("ph1_gain", "Gain distribution;Gain value (ADC count);Count", 200, 700, 2700);
	TH1F* ph1_rc = 			new TH1F("ph1_rc", "RC distribution;RC value (ns/mm^{2});Count", 200, 50, 275);
	TH1F* ph1_rchigh = 		new TH1F("ph1_rchigh", "High RC distribution;RC value (ns/mm^{2});Count", 200, 50, 275);
	TH1F* ph1_rclow = 		new TH1F("ph1_rclow", "Low RC distribution;RC value (ns/mm^{2});Count", 200, 50, 275);
	TH1F *ph1_rcbot = 		new TH1F("ph1_rcbot", "RC distribution (bHAT);RC value (ns/mm^{2});Count", 200, 50, 275);
	TH1F *ph1_rctop = 		new TH1F("ph1_rctop", "RC distribution (tHAT);RC value (ns/mm^{2});Count", 200, 50, 275);
	TH1F* ph1_rcdiff = 		new TH1F("ph1_rcdiff", "RC difference between neighbours;RC difference (ns/mm^{2});Count", 61, -30, 30);
	TH1F* ph1_rcdiffcut = 	new TH1F("ph1_rcdiffcut", "RC difference between neighbours;RC difference (ns/mm^{2});Count", 61, -30, 30);
	TH1F* ph1_rcdiffhori =	new TH1F("ph1_rcdiffhori", "ph1_rcdiffhori", 61, -30, 30);
	TH1F* ph1_rcdiffvert = 	new TH1F("ph1_rcdiffvert", "ph1_rc", 61, -30, 30);

	float rc = 0;
	float gain = 0;
	// Loop to get min and max values of RC in each ERAM
	for(int i=0;i<nerams;i++){
		float min_rc = 1e6;
		float max_rc = 0;
		for(int j=0;j<col;j++){
			for(int k=0;k<row;k++){
				rc = pERAMMaps->RC(i,j,k);
				if(rc < min_rc) min_rc = rc;
				if(rc > max_rc) max_rc = rc;
			}
		}
		v_min_rc.push_back(min_rc);
		v_max_rc.push_back(max_rc);
	}

	// Normalized mean rc value of each pad
	float mean_rc[col][row];
	for(int j=0;j<col;j++) for(int k=0;k<row;k++) mean_rc[j][k] = 0;
	for(int j=0;j<col;j++){
		for(int k=0;k<row;k++){
			for(int i=0;i<nerams;i++){
				mean_rc[j][k] += (pERAMMaps->RC(i,j,k)-v_min_rc[i])/(v_max_rc[i]-v_min_rc[i])/(nerams-1);
			}
		}
	}

	// Mean RC per ERAM
	std::vector<float> v_eram_rc;
	std::vector<int> v_ID_low, v_ID_high, v_channel_low, v_channel_high, v_RCs;
	float eram_rc, eram_npads;
	for(int i=0;i<nerams;i++){
		eram_rc = 0, eram_npads = 0;
		for(int j=0;j<col;j++){
			for(int k=0;k<row;k++){
				rc = 				pERAMMaps->RC(i,j,k);
				if(rc != 0) eram_rc += rc, eram_npads++;
			}
		}
		v_eram_rc.push_back(eram_rc/eram_npads);
		if(eram_rc/eram_npads < 130){
			v_ID_low.push_back(pERAMMaps->ID(i));
			v_channel_low.push_back(i);
			v_RCs.push_back(0);
		}
		else{
			v_ID_high.push_back(pERAMMaps->ID(i));
			v_channel_high.push_back(i);
			v_RCs.push_back(1);
		}
	}
	std::cout << "Low RC ERAMs:" << std::endl;
	std::cout << "Channel: ";
	for(int i=0;i<(int)v_channel_low.size();i++) std::cout << v_channel_low[i] << " ";
	std::cout << std::endl;
	std::cout << "ID:      ";
	for(int i=0;i<(int)v_ID_low.size();i++) std::cout << v_ID_low[i] << " ";
	std::cout << std::endl;
	std::cout << "High RC ERAMs:" << std::endl;
	std::cout << "Channel: ";
	for(int i=0;i<(int)v_channel_high.size();i++) std::cout << v_channel_high[i] << " ";
	std::cout << std::endl;
	std::cout << "ID:      ";
	for(int i=0;i<(int)v_ID_high.size();i++) std::cout << v_ID_high[i] << " ";
	std::cout << std::endl;
	std::cout << "RC:      ";
	for(int i=0;i<(int)v_RCs.size();i++) std::cout << v_RCs[i] << ", ";
	std::cout << std::endl;

	// Maps
	for(int i=0;i<nerams;i++){
		eram_rc = 0, eram_npads = 0;
		for(int j=0;j<col;j++){
			for(int k=0;k<row;k++){
				gain = 				pERAMMaps->Gain(i,j,k);
				rc = 				pERAMMaps->RC(i,j,k);
				ph2_gainmap->		Fill(j,k,gain/(nerams-1));
				ph2_rcmap->			Fill(j,k,rc/(nerams-1));
				ph2_rcnorm->		Fill(j,k,(rc-v_min_rc[i])/(v_max_rc[i]-v_min_rc[i])/(nerams-1));
				csv_rcnorm[j][31-k] += (rc-v_min_rc[i])/(v_max_rc[i]-v_min_rc[i])/(nerams-1);
				if(i<nerams/2){
					ph1_rcbot->		Fill(rc);
					ph2_rcmapbot->	Fill(j,k,rc/(nerams/2-1));
					ph2_rcnormbot->	Fill(j,k,(rc-v_min_rc[i])/(v_max_rc[i]-v_min_rc[i])/(nerams/2-1));
					csv_rcnormbot[j][31-k] += (rc-v_min_rc[i])/(v_max_rc[i]-v_min_rc[i])/(nerams/2-1);
				}
				else{
					ph1_rctop->		Fill(rc);
					ph2_rcmaptop->	Fill(j,k,rc/(nerams/2));
					ph2_rcnormtop->	Fill(j,k,(rc-v_min_rc[i])/(v_max_rc[i]-v_min_rc[i])/(nerams/2));
					csv_rcnormtop[j][31-k] += (rc-v_min_rc[i])/(v_max_rc[i]-v_min_rc[i])/(nerams/2);
				}
				ph2_GRC->			Fill(rc, gain);
				v_ph2_GRC[i]->		Fill(rc, gain);
				ph1_rc->			Fill(rc);
				v_eram_rc[i] < 130 ? ph1_rclow->Fill(rc) : ph1_rchigh->Fill(rc);
				ph1_gain->			Fill(gain);
				if(j>0) ph1_rcdiffhori->Fill(rc-pERAMMaps->RC(i,j-1,k));
				if(k>0) ph1_rcdiffvert->Fill(rc-pERAMMaps->RC(i,j,k-1));
				// Std of RC difference vertical and horizontal with cutoff for differences > 20
				if(fabs(rc-pERAMMaps->RC(i,j-1,k)) < 20 ) ph1_rcdiffcut->Fill(rc-pERAMMaps->RC(i,j-1,k));
				if(fabs(rc-pERAMMaps->RC(i,j,k-1)) < 20 ) ph1_rcdiffcut->Fill(rc-pERAMMaps->RC(i,j,k-1));
			}
		}
	}
	ph1_rcdiff->Add(ph1_rcdiffhori);
	ph1_rcdiff->Add(ph1_rcdiffvert);
	// Fit ph1_rc with two Gaussian functions
	TF1 *gausn1 = new TF1("gaus1", "gausn", 50, 130);
	TF1 *gausn2 = new TF1("gaus2", "gausn", 130, 220);
	ph1_rc->Fit(gausn1, "RQ");
	ph1_rc->Fit(gausn2, "RQ+");
	TF1* doubleGaus = new TF1("doubleGaus", "gausn(0) + gausn(3)", 50, 230);
	doubleGaus->SetParameters(gausn1->GetParameter(0), gausn1->GetParameter(1), gausn1->GetParameter(2), 
							  gausn2->GetParameter(0), gausn2->GetParameter(1), gausn2->GetParameter(2));
	ph1_rc->Fit(doubleGaus, "RQ");

	// Fill CSV files
	for(int k=0;k<row;k++){
		for(int j=0;j<col;j++){
			CSVRCNorm << csv_rcnorm[j][k];
			CSVRCNormBot << csv_rcnormbot[j][k];
			CSVRCNormTop << csv_rcnormtop[j][k];
			if(j < 35){
				CSVRCNorm << ",";
				CSVRCNormBot << ",";
				CSVRCNormTop << ",";
			}
		}
		CSVRCNorm << std::endl;
		CSVRCNormBot << std::endl;
		CSVRCNormTop << std::endl;
	}
	CSVRCNorm.close();
	CSVRCNormBot.close();
	CSVRCNormTop.close();

	// Differences at TPC level
	float rc_top = 0, rc_bot = 0, rc_diff = 0;
	float rc_top_norm = 0, rc_bot_norm = 0, rc_diff_norm = 0;
	for(int j=0;j<col;j++){
		for(int k=0;k<row;k++){
			rc_bot = ph2_rcmapbot->GetBinContent(j+1, k+1);
			rc_top = ph2_rcmaptop->GetBinContent(j+1, k+1);
			rc_diff = rc_top - rc_bot;
			ph2_rcdiffTPC->Fill(j, k, rc_diff);
			rc_bot_norm = ph2_rcnormbot->GetBinContent(j+1, k+1);
			rc_top_norm = ph2_rcnormtop->GetBinContent(j+1, k+1);
			rc_diff_norm = fabs(rc_top_norm - rc_bot_norm)*100;
			ph2_rcnormdiffTPC->Fill(j, k, rc_diff_norm);
		}
	}

	// Standard deviation
	float sigma[col][row];
	for(int j=0;j<col;j++) for(int k=0;k<row;k++) sigma[j][k] = 0;
	for(int j=0;j<col;j++){
		for(int k=0;k<row;k++){
			for(int i=0;i<nerams;i++){
				float rc_normed = (pERAMMaps->RC(i,j,k)-v_min_rc[i])/(v_max_rc[i]-v_min_rc[i]);
				sigma[j][k] += pow(rc_normed-mean_rc[j][k], 2)/(nerams-1);
			}
			ph2_rcnormstd->Fill(j,k,sigma[j][k]*100);
		}
	}

	// Get min and max average values of RC
	std::vector<float> v_min_avg_rc;
	std::vector<float> v_max_avg_rc;
	v_min_avg_rc.push_back(ph2_rcmap->GetMinimum());
	v_max_avg_rc.push_back(ph2_rcmap->GetMaximum());
	v_min_avg_rc.push_back(ph2_rcmapbot->GetMinimum());
	v_max_avg_rc.push_back(ph2_rcmapbot->GetMaximum());
	v_min_avg_rc.push_back(ph2_rcmaptop->GetMinimum());
	v_max_avg_rc.push_back(ph2_rcmaptop->GetMaximum());
	float min_avg_rc = *std::min_element(v_min_avg_rc.begin(), v_min_avg_rc.end());
	float max_avg_rc = *std::max_element(v_max_avg_rc.begin(), v_max_avg_rc.end());


	
	// Drawing
	gStyle->					SetOptStat(0);
	gStyle->					SetOptFit(0);
  	OutputFile = 				"Output/Distributions.pdf";
  	OutputFile_Beg = 			OutputFile + "(" ;
  	OutputFile_End = 			OutputFile + ")" ;

	gPad->						UseCurrentStyle();
	pTCanvas->					cd();
	Graphic_setup(ph1_gain, 0.5, 1, kMagenta+1, 2, kMagenta-2, kMagenta, 0.2);
	ph1_gain->					Draw("hist");
	pTCanvas->					SaveAs(OutputFile_Beg.c_str());
	pTCanvas->					Clear();


	Graphic_setup(ph1_rc, 0.5, 1, kMagenta+1, 2, kMagenta-2, kMagenta, 0.2);
	Graphic_setup(ph1_rctop, 0.5, 1, kRed+1, 2, kRed-2, kRed, 0.2);
	Graphic_setup(ph1_rcbot, 0.5, 1, kBlue+1, 2, kBlue-2, kBlue, 0.2);
	ph1_rc->					Draw("hist");
	ph1_rcbot->					Draw("hist same");
	ph1_rctop->					Draw("hist same");
	doubleGaus->				SetLineColor(kMagenta+3);
	doubleGaus->				Draw("same");
	pTLegend = 					new TLegend(0.77, 0.7, 0.95, 0.89);
	pTLegend->					SetTextSize(0.05);
	pTLegend->					AddEntry(ph1_rc, "Both", "f");
	pTLegend->					AddEntry(ph1_rctop, "tHAT", "f");
	pTLegend->					AddEntry(ph1_rcbot, "bHAT", "f");
	pTLegend->					SetFillStyle(0);
	pTLegend->					SetFillColorAlpha(0, 0);
	pTLegend->					Draw();
	pTCanvas->					SaveAs(OutputFile.c_str());
	pTCanvas->					Clear();

	Graphic_setup(ph1_rchigh, 0.5, 1, kRed+1, 2, kRed-2, kRed, 0.2);
	Graphic_setup(ph1_rclow, 0.5, 1, kBlue+1, 2, kBlue-2, kBlue, 0.2);
	ph1_rclow->					Draw("hist");
	ph1_rchigh->				Draw("hist same");
	ph1_rclow->					SetTitle("RC distributions");
	PrintResolution(ph1_rclow, pTCanvas, 0.65, 0.6, kBlue+1, "Low mean RC ERAMs");
	PrintResolution(ph1_rchigh, pTCanvas, 0.65, 0.25, kRed+1, "High mean RC ERAMs");
	ph1_rclow->					GetFunction("gausn")->SetLineColor(kBlue+3);
	ph1_rclow->					GetFunction("gausn")->Draw("same");
	ph1_rchigh->				GetFunction("gausn")->SetLineColor(kRed+3);
	ph1_rchigh->				GetFunction("gausn")->Draw("same");
	pTCanvas->					SaveAs(OutputFile.c_str());
	pTCanvas->					Clear();	

	Graphic_setup(ph1_rcdiff, 0.5, 1, kMagenta+1, 2, kMagenta-2, kMagenta, 0.2);
	Graphic_setup(ph1_rcdiffhori, 0.5, 1, kRed+1, 2, kRed-2, kRed, 0.2);
	Graphic_setup(ph1_rcdiffvert, 0.5, 1, kBlue+1, 2, kBlue-2, kBlue, 0.2);
	ph1_rcdiff->				Draw("hist");
	ph1_rcdiffhori->			Draw("hist same");
	ph1_rcdiffvert->			Draw("hist same");
	TLatex* pLatexRCDiff = 		new TLatex();
	pLatexRCDiff->				SetTextSize(0.05);
	pLatexRCDiff->				SetTextAlign(32);
	pLatexRCDiff->				SetNDC();
	pLatexRCDiff->				DrawLatex(0.95, 0.65, Form("#mu #pm #sigma = %.2f #pm %.2f", ph1_rcdiffcut->GetMean(), ph1_rcdiffcut->GetRMS()));
	// Add legend
	pTLegend = 					new TLegend(0.67, 0.7, 0.89, 0.89);
	pTLegend->					SetTextSize(0.05);
	pTLegend->					AddEntry(ph1_rcdiff, "RC Difference", "f");
	pTLegend->					AddEntry(ph1_rcdiffhori, "Horizontal", "f");
	pTLegend->					AddEntry(ph1_rcdiffvert, "Vertical", "f");
	pTLegend->					SetFillStyle(0);
	pTLegend->					SetFillColorAlpha(0, 0);
	pTLegend->					Draw();
	pTCanvas->					SaveAs(OutputFile.c_str());
	pTCanvas->					Clear();

	gPad->						SetRightMargin(0.12);
	ph2_GRC->					Draw("colz");
	pTCanvas->					SaveAs(OutputFile.c_str());

	gStyle->					SetPadRightMargin(0.12);
	pTCanvas->					Clear();
	pTCanvas->					Divide(4,4);
	for(int i=0;i<(nerams/2-1);i++){
		pTCanvas->				cd(i+1);
		v_ph2_GRC[i]->			Draw("colz");
	}
	pTCanvas->					SaveAs(OutputFile.c_str());

	pTCanvas->					Clear();
	pTCanvas->					Divide(4,4);
	for(int i=(nerams/2-1);i<nerams;i++){
		pTCanvas->				cd(i-15);
		v_ph2_GRC[i]->			Draw("colz");
	}
	pTCanvas->					SaveAs(OutputFile.c_str());

	pTCanvas->					Clear();
  	gStyle->                    SetPalette(kRainBow);
	ph2_gainmap->				Draw("colz");
	pTCanvas->					SaveAs(OutputFile.c_str());

	pTCanvas->					Clear();
	pTCanvas->					Divide(2,2);
	pTCanvas->					cd(1);
	ph2_rcmap->					SetMinimum(min_avg_rc);
	ph2_rcmap->					SetMaximum(max_avg_rc);
	ph2_rcmap->					Draw("colz");
	pTCanvas->					cd(2);
	ph2_rcmaptop->				SetMinimum(min_avg_rc);
	ph2_rcmaptop->				SetMaximum(max_avg_rc);
	ph2_rcmaptop->				Draw("colz");
	pTCanvas->					cd(3);
	ph2_rcdiffTPC->				Draw("colz");
	pTCanvas->					cd(4);
	ph2_rcmapbot->				SetMinimum(min_avg_rc);
	ph2_rcmapbot->				SetMaximum(max_avg_rc);
	ph2_rcmapbot->				Draw("colz");
	pTCanvas->					SaveAs(OutputFile.c_str());

	pTCanvas->					Clear();
	pTCanvas->					Divide(2,2);
	pTCanvas->					cd(1);
	ph2_rcnorm->				SetMinimum(0);
	ph2_rcnorm->				SetMaximum(1);
	ph2_rcnorm->				Draw("colz");
	pTCanvas->					cd(2);
	ph2_rcnormtop->				SetMinimum(0);
	ph2_rcnormtop->				SetMaximum(1);
	ph2_rcnormtop->				Draw("colz");
	pTCanvas->					cd(3);
	ph2_rcnormdiffTPC->			SetMinimum(0);
	ph2_rcnormdiffTPC->			Draw("colz");
	pTCanvas->					cd(4);
	ph2_rcnormbot->				SetMinimum(0);
	ph2_rcnormbot->				SetMaximum(1);
	ph2_rcnormbot->				Draw("colz");
	pTCanvas->					SaveAs(OutputFile.c_str());

	pTCanvas->					Clear();
	ph2_rcnormstd->				SetMinimum(0);
	ph2_rcnormstd->				Draw("colz");
	pTCanvas->					SaveAs(OutputFile_End.c_str());
}