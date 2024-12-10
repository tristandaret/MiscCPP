#include "Draw.h"
#include "Misc_Functions.h"
#include "TFrame.h"
#include "TPaveText.h"
#include "SetStyle.h"

#include <fstream>

// Constructor
Draw::Draw(){
	fpStyle =						SetMyStyle();
	gROOT->							SetStyle(fpStyle->GetName());
	fpCanvas =						new TCanvas("fpCanvas", "fpCanvas", 1800, 1350);
	gPad->							UseCurrentStyle();
}

// Destructor
Draw::~Draw(){}



void Draw::Run(const std::string &filepath){
	TFile inputFile(filepath.c_str(), "READ");
	Process pr = *(Process*)inputFile.Get("Process");
	fRealpathPDF =             foutputPDFFolder  + "/" + pr.ffileName + ".pdf";
	std::cout << "drawout:" << fRealpathPDF << std::endl;

	// Plot setup ------------------------------------------------------------------------------------------------------------------------------
	int dEdxmax = 					pr.dEdxmax;
	int momrange = 					pr.momrange;
	int nmombins = 					pr.nmombins;
	int nthetabins = 				pr.nthetabins;
	fpCanvas->						cd();
	gStyle->						SetOptStat(0);
	gStyle->						SetOptFit(0);
	fpLegend =						new TLegend(0.6,0.78,0.9,0.91);
	fpLegend->						SetTextSize(0.06);
	fpLegend->						SetFillStyle(0);
	fpLegend->						SetTextColor(kBlue-1);
	TLegendEntry *pentryWF =		fpLegend->AddEntry((TObject*)0, "Waveforms sum", "p");
	pentryWF->						SetMarkerColor(kCyan+2);
	pentryWF->						SetMarkerSize(7);
	pentryWF->						SetMarkerStyle(33);
	TLegendEntry *pentryXP =		fpLegend->AddEntry((TObject*)0, "Crossed pads", "p");
	pentryXP->						SetMarkerColor(kMagenta+2);
	pentryXP->						SetMarkerSize(7);
	pentryXP->						SetMarkerStyle(47);
	gPad->							SetTopMargin(0.05);
	float invX = 					0;

	// Global dE/dx plot -------------------------------------------------------------------------------------------------------------------
	Graphic_setup(pr.fph1f_WF, 0.5, 1, kCyan+1, 2, kCyan-2, kCyan, 0.2);
	Graphic_setup(pr.fph1f_XP, 0.5, 1, kMagenta+2, 2, kMagenta-2, kMagenta, 0.2);
	pr.fph1f_WF->					SetAxisRange(0, 1.1	*std::max({pr.fph1f_WF->GetMaximum(), pr.fph1f_XP->GetMaximum()}),	"Y");
	if(pr.fph1f_WF->GetMean() > dEdxmax/2) invX = 0.4;
	pr.fph1f_WF->					Draw("HIST");
	pr.fph1f_XP->					Draw("HIST sames");
	PrintResolution(pr.fph1f_XP, fpCanvas, 0.65-invX, 0.58, kMagenta+2, "XP");
	PrintResolution(pr.fph1f_WF, fpCanvas, 0.65-invX, 0.25, kCyan+2, "WF");
	fpCanvas->						SaveAs((fRealpathPDF + "(").c_str());

	// dE/dx ERAM by ERAM ----------------------------------------------------------------------------------------------------------------------
	float maxdEdx = 0;
	for (TH1F *hist : pr.vmod_fph1f_WF) if (hist->GetMaximum() > maxdEdx) maxdEdx = hist->GetMaximum();
	for (TH1F *hist : pr.vmod_fph1f_XP) if (hist->GetMaximum() > maxdEdx) maxdEdx = hist->GetMaximum();

	float xMax = 0, yMax = 0, reso=0, dreso = 0;
	fpCanvas->						Clear();
	fpCanvas->						Divide(4,4);
	for(int i=0;i<32;i++){
		pr.vmod_fph1f_WF[i]->		SetAxisRange(0, 1.1	*maxdEdx,	"Y");
		pr.vmod_fph1f_XP[i]->		SetAxisRange(0, 1.1	*maxdEdx,	"Y");
		Graphic_setup(pr.vmod_fph1f_WF[i], 0.5, 1, kCyan+1, 1, kCyan-2, kCyan, 0.2);
		Graphic_setup(pr.vmod_fph1f_XP[i], 0.5, 1, kMagenta+2, 1, kMagenta-2, kMagenta, 0.2);
	}
	for(int i = 0; i < 16; i++){
		fpCanvas->					cd(i+1);
		pr.vmod_fph1f_WF[i]->		Draw("HIST");
		pr.vmod_fph1f_XP[i]->		Draw("HIST same");
		if(pr.vmod_fph1f_WF[i]->GetEntries() < 100) continue;
		xMax = pr.vmod_fph1f_WF[i]->GetXaxis()->GetXmax();
		yMax = pr.vmod_fph1f_WF[i]->GetMaximum();
		pr.vmod_fph1f_WF[i]->GetMean() > dEdxmax/2 ? invX = 0.4 : invX = 0;
		PrintResolution(pr.vmod_fph1f_WF[i], fpCanvas, 0.65-invX, 0.58, kCyan+2, "WF");
		PrintResolution(pr.vmod_fph1f_XP[i], fpCanvas, 0.65-invX, 0.25, kMagenta+2, "XP");
	}
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	for(int i = 16; i < 32; i++){
		fpCanvas->					cd(i+1-16);
		pr.vmod_fph1f_WF[i]->		Draw("HIST");
		pr.vmod_fph1f_XP[i]->		Draw("HIST same");

		if(pr.vmod_fph1f_WF[i]->GetEntries() < 100) continue;
		xMax = pr.vmod_fph1f_WF[i]->GetXaxis()->GetXmax();
		yMax = pr.vmod_fph1f_WF[i]->GetMaximum();
		pr.vmod_fph1f_WF[i]->GetMean() > dEdxmax/2 ? invX = 0.4 : invX = 0;
		PrintResolution(pr.vmod_fph1f_WF[i], fpCanvas, 0.65-invX, 0.58, kCyan+2, "WF");
		PrintResolution(pr.vmod_fph1f_XP[i], fpCanvas, 0.65-invX, 0.25, kMagenta+2, "XP");
	}
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// dE/dx XP vs WF ---------------------------------------------------------------------------------------------------------------------------
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	gStyle->						SetOptStat(111111);
	gStyle->						SetStatX(0.33);
	gStyle->						SetStatY(0.95);
	pr.fph2f_WFXP->					Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// dE/dx per time bin -----------------------------------------------------------------------------------------------------------------------
	// Resolution
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.03);
	Graphic_setup(pr.ptge_dd_reso_WF, 2, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_dd_reso_XP, 2, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_dd_reso_XP->			SetTitle(";Drift time (timebins);Resolution (%)");
	pr.ptge_dd_reso_XP->			GetXaxis()->SetLimits(0, 510);
	pr.ptge_dd_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	pr.ptge_dd_reso_XP->			DrawClone("AP");
	pr.ptge_dd_reso_WF->			DrawClone("P same");
	pr.ptge_dd_reso_XP->			SetMarkerSize(7);
	pr.ptge_dd_reso_WF->			SetMarkerSize(7);
	fpLegend->						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_dd_mean_WF, 0, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_dd_mean_XP, 0, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_dd_mean_XP->			SetTitle(";Drift time (timebins);Mean (ADC counts/cm)");
	pr.ptge_dd_mean_XP->			GetXaxis()->SetLimits(0, 510);
	pr.ptge_dd_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	pr.ptge_dd_mean_XP->			Draw("AP");
	pr.ptge_dd_mean_WF->			Draw("P same");
	fpLegend->						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// 2D distribution
	fpCanvas->						Clear();
	gPad->							SetTopMargin(0.02);
	gPad->							SetRightMargin(0.1);
	gStyle->						SetOptStat(0);
	gStyle->						SetOptFit(0); 
	pr.fph2f_XPdrift->				Draw("colz");
	TGraphErrors *ptge_mom_XPdrift =Convert_TH2_TGE(pr.fph2f_XPdrift);
	TF1 *linearFit =				new TF1("linearFit", "pol1", 50, 250);
	ptge_mom_XPdrift->				Fit(linearFit, "RQ");
	linearFit->						SetLineColor(kRed);
	linearFit->						Draw("same");
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.05);
	latex.SetTextColor(kBlue-1);
	latex.DrawLatex(0.6, 0.3, Form("y = %.2fx + %.2f", linearFit->GetParameter(1), linearFit->GetParameter(0)));
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// 2D distribution module per module
	for(int i = 0; i < 32; i++){
		if(i%16 == 0){
			fpCanvas->				Clear();
			fpCanvas->				Divide(4,4);
		}
		fpCanvas->					cd(i%16+1);
		fpCanvas->cd(i%16+1)->		SetRightMargin(0.1);
		pr.vmod_fph2f_XPtmean[i]->	Draw("colz");
		if(i%16 == 15)fpCanvas->	SaveAs(fRealpathPDF.c_str());
	}

	// dE/dx vs track angle phi -----------------------------------------------------------------------------------------------------------------
	// Resolution
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.03);
	Graphic_setup(pr.ptge_phi_reso_WF, 2, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_phi_reso_XP, 2, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_phi_reso_XP->			SetTitle(";#varphi (#circ);Resolution (%)");
	pr.ptge_phi_reso_XP->			GetXaxis()->SetLimits(-90, 90);
	pr.ptge_phi_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	pr.ptge_phi_reso_XP->			DrawClone("AP");
	pr.ptge_phi_reso_WF->			DrawClone("P same");
	pr.ptge_phi_reso_XP->			SetMarkerSize(7);
	pr.ptge_phi_reso_WF->			SetMarkerSize(7);
	fpLegend->						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_phi_mean_WF, 0, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_phi_mean_XP, 0, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_phi_mean_XP->			SetTitle(";#varphi (#circ);Mean (ADC counts/cm)");
	pr.ptge_phi_mean_XP->			GetXaxis()->SetLimits(-90, 90);
	pr.ptge_phi_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	pr.ptge_phi_mean_XP->			Draw("AP");
	pr.ptge_phi_mean_WF->			Draw("P same");
	fpLegend->						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// 2D distribution
	fpCanvas->						Clear();
	gPad->							SetTopMargin(0.02);
	gPad->							SetRightMargin(0.13);
	gStyle->						SetOptStat(0);
	gStyle->						SetOptFit(0);
	pr.fph2f_XPphi->				SetTitle(";#varphi (#circ);dE/dx with XP (ADC counts/cm)");
	pr.fph2f_XPphi->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	pr.fph2f_lenphi->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	pr.fph2f_chi2ndfphi->			Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.02);
	pr.fph1f_phi->					Draw("HIST");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// dE/dx vs track angle theta -----------------------------------------------------------------------------------------------------------------
	// Resolution
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.03);
	Graphic_setup(pr.ptge_theta_reso_WF, 2, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_theta_reso_XP, 2, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_theta_reso_XP->			SetTitle(";#theta (#circ);Resolution (%)");
	pr.ptge_theta_reso_XP->			GetXaxis()->SetLimits(-90, 90);
	pr.ptge_theta_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	pr.ptge_theta_reso_XP->			DrawClone("AP");
	pr.ptge_theta_reso_WF->			DrawClone("P same");
	pr.ptge_theta_reso_XP->			SetMarkerSize(7);
	pr.ptge_theta_reso_WF->			SetMarkerSize(7);
	fpLegend->						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_theta_mean_WF, 0, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_theta_mean_XP, 0, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_theta_mean_XP->			SetTitle(";#theta (#circ);Mean (ADC counts/cm)");
	pr.ptge_theta_mean_XP->			GetXaxis()->SetLimits(-90, 90);
	pr.ptge_theta_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	pr.ptge_theta_mean_XP->			Draw("AP");
	pr.ptge_theta_mean_WF->			Draw("P same");
	fpLegend->						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// 2D distribution
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	pr.fph2f_XPtheta->				SetTitle(";#theta (#circ);dE/dx with XP (ADC counts/cm)");
	pr.fph2f_XPtheta->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	pr.fph2f_lentheta->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	pr.fph2f_momtheta->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.03);
	pr.fph1f_theta->				Draw("HIST");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());


	// Both angles -----------------------------------------------------------------------------------------------------------------------------
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	pr.fph2f_phitheta->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// dEdx vs momentum --------------------------------------------------------------------------------------------------------------------
	// Resolution
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.03);
	Graphic_setup(pr.ptge_mom_reso_WF, 2, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_mom_reso_XP, 2, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_mom_reso_XP->			SetTitle(";Momentum (MeV/c);Resolution (%)");
	pr.ptge_mom_reso_XP->			GetXaxis()->SetLimits(-momrange, momrange);
	pr.ptge_mom_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	pr.ptge_mom_reso_XP->			DrawClone("AP");
	pr.ptge_mom_reso_WF->			DrawClone("P same");
	pr.ptge_mom_reso_XP->			SetMarkerSize(7);
	pr.ptge_mom_reso_WF->			SetMarkerSize(7);
	fpLegend->						Draw();
	fpCanvas->						SaveAs((fRealpathPDF + "(").c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_mom_mean_WF, 0, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_mom_mean_XP, 0, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_mom_mean_XP->			SetTitle(";Momentum (MeV/c);Mean (ADC counts/cm)");
	pr.ptge_mom_mean_XP->			GetXaxis()->SetLimits(-momrange, momrange);
	pr.ptge_mom_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	pr.ptge_mom_mean_XP->			Draw("AP");
	pr.ptge_mom_mean_WF->			Draw("P same");
	fpLegend->						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// 2D distribution
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	pr.fph2f_XPmom->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	pr.fph2f_WFmom->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	gStyle->						SetOptStat(0);
	gPad->							SetRightMargin(0.03);
	Graphic_setup(pr.fph1i_mom, 	 0.5, 1, kMagenta+2, 2, kMagenta-2, kMagenta, 0.2);
	Graphic_setup(pr.fph1i_mom_tHAT, 0.5, 1, kRed+2, 2, kRed-2, kRed, 0.2);
	Graphic_setup(pr.fph1i_mom_bHAT, 0.5, 1, kBlue+2, 2, kBlue-2, kBlue, 0.2);
	// pr.fph1i_mom->				SetAxisRange(0, 1e4, "Y");
	pr.fph1i_mom->					Draw();
	pr.fph1i_mom_bHAT->				Draw("same");
	pr.fph1i_mom_tHAT->				Draw("same");
	TLegend *plegmom = 				new TLegend(0.7, 0.7, 0.9, 0.9);
	plegmom->						AddEntry(pr.fph1i_mom, "Both HATs", "l");
	plegmom->						AddEntry(pr.fph1i_mom_tHAT, "tHAT", "l");
	plegmom->						AddEntry(pr.fph1i_mom_bHAT, "bHAT", "l");
	plegmom->						SetTextSize(0.06);
	plegmom->						SetTextColor(kBlue-1);
	plegmom->						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());
	delete plegmom;

	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	pr.fph2f_chi2mom->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Pulls  -----------------------------------------------------------------------------------------------------------------------
	// Pulls vs theta ---------------------------------------------------------------------------------------------------------------
	// Standard deviation
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.03);
	Graphic_setup(pr.ptge_theta_std_pullmu, 2, 21, kBlue+1, 2, kBlue+1);
	Graphic_setup(pr.ptge_theta_std_pullelec, 2, 20, kOrange+7, 2, kOrange+7);
	pr.ptge_theta_std_pullmu->		SetTitle(";#theta (#circ);Pulls standard deviation");
	float minpull =					0; 
	float maxpull = 				(pr.ptge_theta_std_pullmu->GetHistogram()->GetMaximum() > pr.ptge_theta_std_pullelec->GetHistogram()->GetMaximum()) ? pr.ptge_theta_std_pullmu->GetHistogram()->GetMaximum() : pr.ptge_theta_std_pullelec->GetHistogram()->GetMaximum();
	float diffpull =				maxpull - minpull;
	pr.ptge_theta_std_pullmu->		GetXaxis()->SetLimits(-90, 90);
	pr.ptge_theta_std_pullmu->		GetYaxis()->SetRangeUser(0, maxpull+0.1*diffpull);
	pr.ptge_theta_std_pullmu->		DrawClone("AP");
	pr.ptge_theta_std_pullelec->	DrawClone("P same");
	pr.ptge_theta_std_pullmu->		SetMarkerSize(7);
	pr.ptge_theta_std_pullelec->	SetMarkerSize(7);
	TLegend legpull(0.45, 0.75, 0.65, 0.95);
	legpull.SetTextAlign(22);
	legpull.						SetNColumns(2);
	legpull.						AddEntry(pr.ptge_theta_std_pullmu, " #mu^{#pm}", "p");
	legpull.						AddEntry(pr.ptge_theta_std_pullelec, " e^{#pm}", "p");
	legpull.						SetTextSize(0.06);
	legpull.						SetTextColor(kBlue-1);
	legpull.						SetFillStyle(0);
	legpull.						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_theta_mean_pullmu, 2, 21, kBlue+1, 2, kBlue+1);
	Graphic_setup(pr.ptge_theta_mean_pullelec, 2, 20, kOrange+7, 2, kOrange+7);
	pr.ptge_theta_mean_pullmu->		SetTitle(";#theta (#circ);Pulls mean");
	minpull = 						(pr.ptge_theta_mean_pullmu->GetHistogram()->GetMinimum() < pr.ptge_theta_mean_pullelec->GetHistogram()->GetMinimum()) ? pr.ptge_theta_mean_pullmu->GetHistogram()->GetMinimum() : pr.ptge_theta_mean_pullelec->GetHistogram()->GetMinimum();
	maxpull = 						(pr.ptge_theta_mean_pullmu->GetHistogram()->GetMaximum() > pr.ptge_theta_mean_pullelec->GetHistogram()->GetMaximum()) ? pr.ptge_theta_mean_pullmu->GetHistogram()->GetMaximum() : pr.ptge_theta_mean_pullelec->GetHistogram()->GetMaximum();
	pr.ptge_theta_mean_pullmu->		GetXaxis()->SetLimits(-90, 90);
	pr.ptge_theta_mean_pullmu->		GetYaxis()->SetRangeUser(minpull-0.1*diffpull, maxpull+0.1*diffpull);
	pr.ptge_theta_mean_pullmu->		Draw("AP");
	pr.ptge_theta_mean_pullelec->	Draw("P same");
	legpull.						SetX1NDC(0.15);
	legpull.						SetX2NDC(0.35);
	legpull.						SetY1NDC(0.1);
	legpull.						SetY2NDC(0.3);
	legpull.						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Pulls vs momentum ---------------------------------------------------------------------------------------------------------------
	// Standard deviation
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.03);
	Graphic_setup(pr.ptge_mom_std_pullmu, 2, 21, kBlue+1, 2, kBlue+1);
	Graphic_setup(pr.ptge_mom_std_pullelec, 2, 20, kOrange+7, 2, kOrange+7);
	pr.ptge_mom_std_pullmu->		SetTitle(";Momentum (MeV/c);Pulls standard deviation");
	minpull =						0;
	maxpull = 						(pr.ptge_mom_std_pullmu->GetHistogram()->GetMaximum() > pr.ptge_mom_std_pullelec->GetHistogram()->GetMaximum()) ? pr.ptge_mom_std_pullmu->GetHistogram()->GetMaximum() : pr.ptge_mom_std_pullelec->GetHistogram()->GetMaximum();
	pr.ptge_mom_std_pullmu->		GetXaxis()->SetLimits(-momrange, momrange);
	pr.ptge_mom_std_pullmu->		GetYaxis()->SetRangeUser(0, maxpull+0.1*diffpull);
	pr.ptge_mom_std_pullmu->		DrawClone("AP");
	pr.ptge_mom_std_pullelec->		DrawClone("P same");
	pr.ptge_mom_std_pullmu->		SetMarkerSize(7);
	pr.ptge_mom_std_pullelec->		SetMarkerSize(7);
	legpull.SetX1NDC(0.75); legpull.SetX2NDC(0.95); legpull.SetY1NDC(0.75); legpull.SetY2NDC(0.9);
	legpull.						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_mom_mean_pullmu, 2, 21, kBlue+1, 2, kBlue+1);
	Graphic_setup(pr.ptge_mom_mean_pullelec, 2, 20, kOrange+7, 2, kOrange+7);
	pr.ptge_mom_mean_pullmu->		SetTitle(";Momentum (MeV/c);Pulls mean");
	minpull = 						(pr.ptge_mom_mean_pullmu->GetHistogram()->GetMinimum() < pr.ptge_mom_mean_pullelec->GetHistogram()->GetMinimum()) ? pr.ptge_mom_mean_pullmu->GetHistogram()->GetMinimum() : pr.ptge_mom_mean_pullelec->GetHistogram()->GetMinimum();
	maxpull = 						(pr.ptge_mom_mean_pullmu->GetHistogram()->GetMaximum() > pr.ptge_mom_mean_pullelec->GetHistogram()->GetMaximum()) ? pr.ptge_mom_mean_pullmu->GetHistogram()->GetMaximum() : pr.ptge_mom_mean_pullelec->GetHistogram()->GetMaximum();
	pr.ptge_mom_mean_pullmu->		GetXaxis()->SetLimits(-momrange, momrange);
	pr.ptge_mom_mean_pullmu->		GetYaxis()->SetRangeUser(minpull-0.1*diffpull, maxpull+0.1*diffpull);
	pr.ptge_mom_mean_pullmu->		Draw("AP");
	pr.ptge_mom_mean_pullelec->		Draw("P same");
	legpull.SetX1NDC(0.75); legpull.SetX2NDC(0.95); legpull.SetY1NDC(0.25); legpull.SetY2NDC(0.45);
	legpull.						Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	pr.fph2f_pullelecmu->			Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Pulls distribution 
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.03);
	gPad->							SetTopMargin(0.05);
	Graphic_setup(pr.fph1f_pullmu, 0.5, 1, kCyan+1, 2, kCyan-2, kCyan, 0.2);
	Graphic_setup(pr.fph1f_pullelec, 0.5, 1, kOrange+7, 2, kOrange-2, kOrange, 0.2);
	Graphic_setup(pr.fph1f_pullproton, 0.5, 1, kRed+1, 2, kRed-2, kRed, 0.2);
	TLegend legpull2(0.75, 0.7, 0.9, 0.93);
	legpull2.						AddEntry(pr.fph1f_pullmu, " #mu hyp.", "f");
	legpull2.						AddEntry(pr.fph1f_pullelec, " e hyp.", "f");
	legpull2.						AddEntry(pr.fph1f_pullproton, " p hyp.", "f");
	legpull2.						SetTextSize(0.06);
	legpull2.						SetTextColor(kBlue-1);
	pr.fph1f_pullmu->				SetXTitle("Pull");
	pr.fph1f_pullmu->				Draw("HIST");
	pr.fph1f_pullelec->				Draw("HIST same");
	pr.fph1f_pullproton->			Draw("HIST same");
	TF1 tf1_mu = 					*Fit1Gauss(pr.fph1f_pullmu);
	TF1 tf1_elec = 					*Fit1Gauss(pr.fph1f_pullelec);
	float meanmu = 					tf1_mu.GetParameter(1);
	float dmeanmu = 				tf1_mu.GetParError(1);
	float meanelec = 				tf1_elec.GetParameter(1);
	float dmeanelec = 				tf1_elec.GetParError(1);
	float stdmu = 					tf1_mu.GetParameter(2);
	float dstdmu = 					tf1_mu.GetParError(2);
	float stdelec = 				tf1_elec.GetParameter(2);
	float dstdelec = 				tf1_elec.GetParError(2);
	float separation = 				GetSeparation(&tf1_mu, &tf1_elec);
	float dseparation = 			GetSeparationError(&tf1_mu, &tf1_elec);
	TPaveText pavetextmu(0.6, 0.4, 0.93, 0.55, "NDC");
	pavetextmu.						SetTextSize(0.05);
	pavetextmu.						SetTextColor(kCyan+4);
	pavetextmu.						SetFillColorAlpha(kCyan, 0.2);
	pavetextmu.						SetLineWidth(2);
	pavetextmu.						SetBorderSize(2);
	pavetextmu.						SetLineColor(kCyan-2);
	pavetextmu.						AddText(Form("#mu_{#mu} = %.3f #pm %.3f", meanmu, dmeanmu));
	pavetextmu.						AddText(Form("#sigma_{#mu} = %.3f #pm %.3f", stdmu, dstdmu));

	TPaveText pavetextelec(0.6, 0.2, 0.93, 0.35, "NDC");
	pavetextelec.					SetTextSize(0.05);
	pavetextelec.					SetTextColor(kOrange+3);
	pavetextelec.					SetFillColorAlpha(kOrange, 0.2);
	pavetextelec.					SetLineWidth(2);
	pavetextelec.					SetBorderSize(2);
	pavetextelec.					SetLineColor(kOrange-2);
	pavetextelec.					AddText(Form("#mu_{e} = %.3f #pm %.3f", meanelec, dmeanelec));
	pavetextelec.					AddText(Form("#sigma_{e} = %.3f #pm %.3f", stdelec, dstdelec));

	TLatex latexpull;
	latexpull.						SetNDC();
	latexpull.						SetTextSize(0.06);
	latexpull.						SetTextColor(kBlue-1);
	latexpull.						DrawLatex(0.6, 0.6, Form("S(#mu/e) = %.2f #pm %.2f", separation, dseparation));
	legpull2.						Draw();
	pavetextmu.						Draw();
	pavetextelec.					Draw();
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Geometry plots ----------------------------------------------------------------------------------------------------------------------------

	fpCanvas->						Clear();
	gPad->							SetTopMargin(0.02);
	gPad->							SetLogz();
	pr.fph2f_XZ->					Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());
	gPad->							SetLogz(0);

	// dE/dx vs track length
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	gPad->							SetTopMargin(0.02);
	gStyle->						SetStatX(0.87);
	pr.fph2f_XPlen->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.02);
	pr.fph1f_trklen->				Draw("HIST");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	pr.fph1f_chi2->					Draw("HIST");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// Debugging ------------------------------------------------------------------------------------------------------------------------------

	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	pr.fph2f_momR->					Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	fpCanvas->						Clear();
	pr.fph2f_chi2ndfR->				Draw("colz");
	fpCanvas->						SaveAs(fRealpathPDF.c_str());

	// dE/dx per momentum bin --------------------------------------------------------------------------------------------------------------------
	gPad->							SetTopMargin(0.08);
	gPad->							SetRightMargin(0.02);
	for(int i=0;i<nmombins;i++){
		fpCanvas->					Clear();
		pr.vmom_fph1f_WF[i]->		SetAxisRange(0, 1.1	*std::max({pr.vmom_fph1f_WF[i]->GetMaximum(), pr.vmom_fph1f_XP[i]->GetMaximum()}),	"Y");
		float invX = 				0;
		if(pr.vmom_fph1f_WF[i]->GetMean() > dEdxmax/2) invX = 0.4;
		Graphic_setup(pr.vmom_fph1f_WF[i], 0.5, 1, kCyan+1, 2, kCyan-2, kCyan, 0.2);
		Graphic_setup(pr.vmom_fph1f_XP[i], 0.5, 1, kMagenta+1, 2, kMagenta-2, kMagenta, 0.2);
		pr.vmom_fph1f_WF[i]->				Draw("HIST");
		pr.vmom_fph1f_XP[i]->				Draw("HIST same");
		PrintResolution(pr.vmom_fph1f_XP[i], fpCanvas, 0.65-invX, 0.58, kMagenta+2, "XP");
		PrintResolution(pr.vmom_fph1f_WF[i], fpCanvas, 0.65-invX, 0.25, kCyan+2, "WF");
		fpCanvas->					SaveAs(fRealpathPDF.c_str());
	}

	// dE/dx per theta bin --------------------------------------------------------------------------------------------------------------------
	gPad->							SetTopMargin(0.08);
	gPad->							SetRightMargin(0.02);
	for(int i=0;i<nthetabins;i++){
		fpCanvas->					Clear();
		pr.vtheta_fph1f_WF[i]->		SetAxisRange(0, 1.1	*std::max({pr.vtheta_fph1f_WF[i]->GetMaximum(), pr.vtheta_fph1f_XP[i]->GetMaximum()}),	"Y");
		float invX = 				0;
		if(pr.vtheta_fph1f_WF[i]->GetMean() > dEdxmax/2) invX = 0.4;
		Graphic_setup(pr.vtheta_fph1f_WF[i], 0.5, 1, kCyan+1, 2, kCyan-2, kCyan, 0.2);
		Graphic_setup(pr.vtheta_fph1f_XP[i], 0.5, 1, kMagenta+1, 2, kMagenta-2, kMagenta, 0.2);
		pr.vtheta_fph1f_WF[i]->				Draw("");
		pr.vtheta_fph1f_XP[i]->				Draw("same");
		PrintResolution(pr.vtheta_fph1f_XP[i], fpCanvas, 0.65-invX, 0.58, kMagenta+2, "XP");
		PrintResolution(pr.vtheta_fph1f_WF[i], fpCanvas, 0.65-invX, 0.25, kCyan+2, "WF");
		fpCanvas->					SaveAs(fRealpathPDF.c_str());
	}

	// Drift velocity plots ----------------------------------------------------------------------------------------------------------------------
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.02);
	int bins =						100;
	int tmin =						-999;
	int max =						-999;
	for(int i=0;i<bins;i++) 
	{
		if(pr.fph1i_tminEP0->GetBinContent(i) > max){
			tmin = i;
			max = pr.fph1i_tminEP0->GetBinContent(i);
		}
	}

	int tmax = -999;
	for(int i=300;i<300+bins;i++) 
	{
		if(pr.fph1i_tmaxBotCath->GetBinContent(i) > max){
			tmax = i;
			max = pr.fph1i_tmaxBotCath->GetBinContent(i);
		}
	}
	std::cout << std::fixed << std::setprecision(3);
	std::cout << "tmin, tmax: " << tmin << " " << tmax << std::endl;


	// pr.fph1i_tminEP0->				Fit("gaus", "RQ", "", tmin-2, tmin+2);
	// pr.fph1i_tminEP2->				Fit("gaus", "RQ", "", tmin-2, tmin+2);
	// pr.fph1i_tmaxBotCath->			Fit("gaus", "RQ", "", tmax-4, tmax+4);
	// pr.fph1i_tmaxTopCath->			Fit("gaus", "RQ", "", tmax-4, tmax+4);
	// TF1 *ptf1_tminEP0 =			pr.fph1i_tminEP0->GetFunction("gaus");
	// TF1 *ptf1_tminEP2 =			pr.fph1i_tminEP2->GetFunction("gaus");
	// TF1 *ptf1_tmaxBotCath =		pr.fph1i_tmaxBotCath->GetFunction("gaus");
	// TF1 *ptf1_tmaxTopCath =		pr.fph1i_tmaxTopCath->GetFunction("gaus");

	// float Tmin =				ptf1_tminEP0->GetParameter(1)*40/1000;
	// float dTmin =				ptf1_tminEP0->GetParError(1)*40/1000;
	// float Tmax =				ptf1_tmaxBotCath->GetParameter(1)*40/1000;
	// float dTmax =				ptf1_tmaxBotCath->GetParError(1)*40/1000;
	// float Tdrift =				(Tmax - Tmin);
	// float dTdrift =				sqrt(dTmin*dTmin + dTmax*dTmax);
	// std::cout << "Drift time in bHATPC: " << Tdrift << " +/- " << dTdrift << " ns" << std::endl;
	// float driftVel =			98.5/Tdrift;
	// float ddriftVel =			driftVel * dTdrift/Tdrift;
	// std::cout << "Drift velocity in bHATPC: " << driftVel << " +/- " << ddriftVel << " cm/µs" << std::endl;

	// Tmin =						ptf1_tminEP2->GetParameter(1)*40/1000;
	// dTmin =						ptf1_tminEP2->GetParError(1)*40/1000;
	// Tmax =						ptf1_tmaxTopCath->GetParameter(1)*40/1000;
	// dTmax =						ptf1_tmaxTopCath->GetParError(1)*40/1000;
	// Tdrift =					(Tmax - Tmin);
	// dTdrift =					sqrt(dTmin*dTmin + dTmax*dTmax);
	// std::cout << "Drift time in tHATPC: " << Tdrift << " +/- " << dTdrift << " ns" << std::endl;
	// driftVel =					98.5/Tdrift;
	// ddriftVel =					driftVel * dTdrift/Tdrift;
	// std::cout << "Drift velocity in tHATPC: " << driftVel << " +/- " << ddriftVel << " cm/µs" << std::endl;

	pr.fph1i_tminBotCath->			SetAxisRange(0, 1.1	*std::max({pr.fph1i_tminBotCath->GetMaximum(), pr.fph1i_tminEP0->GetMaximum(), pr.fph1i_tminEP1->GetMaximum(), pr.fph1i_tminTopCath->GetMaximum(), pr.fph1i_tminEP2->GetMaximum(), pr.fph1i_tminEP3->GetMaximum()}),	"Y");
	pr.fph1i_tminBotCath->			SetLineWidth(2);
	pr.fph1i_tminBotCath->			SetLineColor(kGreen+2);
	pr.fph1i_tminEP0->				SetLineWidth(2);
	pr.fph1i_tminEP0->				SetLineColor(kBlue);
	pr.fph1i_tminEP1->				SetLineWidth(2);
	pr.fph1i_tminEP1->				SetLineColor(kRed);
	pr.fph1i_tminTopCath->			SetLineWidth(2);
	pr.fph1i_tminTopCath->			SetLineColor(kGreen-6);
	pr.fph1i_tminEP2->				SetLineWidth(2);
	pr.fph1i_tminEP2->				SetLineColor(kBlue-6);
	pr.fph1i_tminEP3->				SetLineWidth(2);
	pr.fph1i_tminEP3->				SetLineColor(kRed-6);
	pr.fph1i_tminBotCath->			Draw();
	pr.fph1i_tminEP0->				Draw("same");	
	pr.fph1i_tminEP1->				Draw("same");	
	pr.fph1i_tminTopCath->			Draw("same");
	pr.fph1i_tminEP2->				Draw("same");
	pr.fph1i_tminEP3->				Draw("same"); 
	// ptf1_tminEP0->				SetLineColor(kOrange+2);
	// ptf1_tminEP2->				SetLineColor(kViolet+2);
	// ptf1_tminEP0->				Draw("same");
	// ptf1_tminEP2->				Draw("same");
	TLegend *legtmin =			new TLegend(0.6,0.5,0.9,0.9); 
	legtmin->					AddEntry(pr.fph1i_tminBotCath, "T_{min} (bottom cathode) ", "l");	
	legtmin->					AddEntry(pr.fph1i_tminEP0, "T_{min} EP0 ", "l");	
	legtmin->					AddEntry(pr.fph1i_tminEP1, "T_{min} EP1 ", "l"); 
	legtmin->					AddEntry(pr.fph1i_tminTopCath, "T_{min} (top cathode) ", "l");
	legtmin->					AddEntry(pr.fph1i_tminEP2, "T_{min} EP2 ", "l");
	legtmin->					AddEntry(pr.fph1i_tminEP3, "T_{min} EP3 ", "l");
	legtmin->					DrawClone();
	fpCanvas->					Update();
	delete legtmin;
	fpCanvas->					SaveAs(fRealpathPDF.c_str());


	fpCanvas->					Clear();
	gStyle->					SetOptStat(0);
	gStyle->					SetOptFit(0); 

	pr.fph1i_tmaxBotCath->			SetAxisRange(0, 1.1	*std::max({pr.fph1i_tmaxBotCath->GetMaximum(), pr.fph1i_tmaxEP0->GetMaximum(), pr.fph1i_tmaxEP1->GetMaximum(), pr.fph1i_tmaxTopCath->GetMaximum(), pr.fph1i_tmaxEP2->GetMaximum(), pr.fph1i_tmaxEP3->GetMaximum()}),	"Y");
	pr.fph1i_tmaxBotCath->			SetLineWidth(2);
	pr.fph1i_tmaxBotCath->			SetLineColor(kGreen+2);
	pr.fph1i_tmaxEP0->				SetLineWidth(2);
	pr.fph1i_tmaxEP0->				SetLineColor(kBlue);
	pr.fph1i_tmaxEP1->				SetLineWidth(2);
	pr.fph1i_tmaxEP1->				SetLineColor(kRed);
	pr.fph1i_tmaxTopCath->			SetLineWidth(2);
	pr.fph1i_tmaxTopCath->			SetLineColor(kGreen-6);
	pr.fph1i_tmaxEP2->				SetLineWidth(2);
	pr.fph1i_tmaxEP2->				SetLineColor(kBlue-6);
	pr.fph1i_tmaxEP3->				SetLineWidth(2);
	pr.fph1i_tmaxEP3->				SetLineColor(kRed-6);
	pr.fph1i_tmaxBotCath->			Draw();
	pr.fph1i_tmaxEP0->				Draw("same");
	pr.fph1i_tmaxEP1->				Draw("same");
	pr.fph1i_tmaxTopCath->			Draw("same");
	pr.fph1i_tmaxEP2->				Draw("same");
	pr.fph1i_tmaxEP3->				Draw("same");
	// ptf1_tmaxBotCath->			SetLineColor(kOrange+2);
	// ptf1_tmaxTopCath->			SetLineColor(kViolet+2);
	// ptf1_tmaxTopCath->			Draw("same");
	// ptf1_tmaxBotCath->			Draw("same");
	TLegend *legtmax =				new TLegend(0.15,0.5,0.4,0.9);
	legtmax->						AddEntry(pr.fph1i_tmaxBotCath, "T_{max} (bottom cathode) ", "l");
	legtmax->						AddEntry(pr.fph1i_tmaxEP0, "T_{max} EP0 ", "l");
	legtmax->						AddEntry(pr.fph1i_tmaxEP1, "T_{max} EP1 ", "l");
	legtmax->						AddEntry(pr.fph1i_tmaxTopCath, "T_{max} (top cathode) ", "l");
	legtmax->						AddEntry(pr.fph1i_tmaxEP2, "T_{max} EP2 ", "l");
	legtmax->						AddEntry(pr.fph1i_tmaxEP3, "T_{max} EP3 ", "l");
	legtmax->						DrawClone();
	fpCanvas->						Update();
	delete legtmax;	
	fpCanvas->						SaveAs((fRealpathPDF + ")").c_str());
	fpCanvas->						Clear();

	// Reset style
	gPad->							UseCurrentStyle();
}


void Draw::CompareRuns(const std::vector<std::string> &v_filepaths, const std::string &type){

	// Load processes
	std::vector<Process*> v_processes;
	for(int i=0; i < v_filepaths.size(); i++) {
		std::cout << "File " << i << ": " << v_filepaths[i] << std::endl;
		TFile inputFile(v_filepaths[i].c_str(), "READ");
		Process* process = dynamic_cast<Process*>(inputFile.Get("Process"));
		v_processes.push_back(process);
	}

	Process *pr0 = v_processes[0];
	int ncomparisons =				(int)v_processes.size();
	std::string OutputFile;
	if(type == "comments"){
		OutputFile =									foutputComparisonFolder + "/Comparison_" + pr0->frun + pr0->ftag + pr0->fcomment;
		for(int i=1;i<ncomparisons;i++) OutputFile += 	"_VS" + v_processes[i]->fcomment;
		OutputFile +=									pr0->fcutslist + ".pdf";
	}
	else if(type == "cuts"){
		OutputFile =									foutputComparisonFolder + "/" + pr0->fcomment + "/Comparison_" + pr0->frun + pr0->ftag + pr0->fcomment + pr0->fcutslist;
		for(int i=1;i<ncomparisons;i++) OutputFile += 	"_VS" + v_processes[i]->fcutslist;
		OutputFile +=									".pdf";
	}


	// Plot setup ------------------------------------------------------------------------------------------------------------------------------
	int dEdxmax = 						pr0->dEdxmax;
	int momrange = 					pr0->momrange;
	int nmombins = 					pr0->nmombins;
	fpCanvas->						cd();
	gStyle->						SetOptStat(0);
	gStyle->						SetOptFit(0);
	gPad->							SetTopMargin(0.05);
	gPad->							SetRightMargin(0.03);

	// Legend
	fpLegend =						new TLegend(0.6,0.85-0.04*ncomparisons,0.91,0.92);
	fpLegend->						SetTextSize(0.06-0.005*ncomparisons);
	// fpLegend->						SetFillStyle(0);
	fpLegend->						SetTextColor(kBlue-1);

	float ampmax = 0;
	for(int i=0;i<ncomparisons;i++){
		Process *pr_tmp = v_processes[i];
		Graphic_setup(pr_tmp->fph1f_WF, 0.5, 1, colors[i], 2, colors[i], colors[i]-2, 0.2);
		Graphic_setup(pr_tmp->fph1f_XP, 0.5, 1, colors[i], 2, colors[i], colors[i]-2, 0.2);
		if(pr_tmp->fph1f_WF->GetMaximum() > ampmax) ampmax = pr_tmp->fph1f_WF->GetMaximum();
		if(pr_tmp->fph1f_XP->GetMaximum() > ampmax) ampmax = pr_tmp->fph1f_XP->GetMaximum();
		Graphic_setup(pr_tmp->ptge_mom_reso_WF,   2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_mom_reso_XP,   2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_mom_mean_WF,   2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_mom_mean_XP,   2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_dd_reso_WF,    2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_dd_reso_XP,    2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_dd_mean_WF,    2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_dd_mean_XP,    2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_phi_reso_WF,   2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_phi_reso_XP,   2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_phi_mean_WF,   2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_phi_mean_XP,   2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_theta_reso_WF, 2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_theta_reso_XP, 2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_theta_mean_WF, 2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp->ptge_theta_mean_XP, 2, markers[i], colors[i], 2, colors[i]);
	}

	std::vector<std::string> legEntries;
	// for(int i=0;i<ncomparisons;i++) legEntries.push_back(type == "comments" ? v_processes[i]->fcomment : v_processes[i]->fcutslist);
	legEntries.push_back("No selections");
	legEntries.push_back("#chi^{2}/ndf < 5 + L > 25 cm");
	legEntries.push_back("+ 200 < P < 1000 MeV/c");
	legEntries.push_back("+ 250 < P < 500 MeV/c");

	// Global dE/dx plot -------------------------------------------------------------------------------------------------------------------
	pr0->fph1f_XP->	SetAxisRange(0, 1.1	*ampmax,	"Y");
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->fph1f_XP->	Draw(i==0 ? "HIST" : "HIST same");
		float xwidth =	0.23;
		float xright =	0.94;
		float ywidth =	(0.92-0.15)/ncomparisons - 0.01;
		float ytop =	0.92 - i*(ywidth + 0.01);
		std::cout << i << " " << xright-xwidth << " < x < " << xright << " | " << ytop-ywidth << " < y < " << ytop << std::endl;
		PrintResolution(v_processes[i]->fph1f_XP, fpCanvas, xright, ytop, xwidth, ywidth, "north east", colors[i], legEntries[i]);
	}
	fpCanvas->						SaveAs((OutputFile + "(").c_str());
	fpCanvas->						Clear();

	// dE/dx vs momentum -------------------------------------------------------------------------------------------------------------------
	// Resolution
	pr0->ptge_mom_reso_XP->			SetTitle(";Momentum (MeV/c);Resolution (%)");
	pr0->ptge_mom_reso_XP->			GetXaxis()->SetLimits(-momrange, momrange);
	pr0->ptge_mom_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_mom_reso_XP->	DrawClone(i==0 ? "AP" : "P same");
		v_processes[i]->ptge_mom_reso_XP->	SetMarkerSize(7);
		fpLegend->AddEntry(v_processes[i]->ptge_mom_reso_XP,  legEntries[i].c_str(), "p");
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs(OutputFile.c_str());
	fpCanvas-> 						Clear();	

	// Resolution's mean
	pr0->ptge_mom_mean_XP->			SetTitle(";Momentum (MeV/c);Mean (ADC counts/cm)");
	pr0->ptge_mom_mean_XP->			GetXaxis()->SetLimits(-momrange, momrange);
	pr0->ptge_mom_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_mom_mean_XP->	DrawClone(i==0 ? "AP" : "P same");
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs(OutputFile.c_str());
	fpCanvas-> 						Clear();

	// dE/dx per time bins -------------------------------------------------------------------------------------------------------------------
	// Resolution
	pr0->ptge_dd_reso_XP->			SetTitle(";Drift time (timebins);Resolution (%)");
	pr0->ptge_dd_reso_XP->			GetXaxis()->SetLimits(0, 510);
	pr0->ptge_dd_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_dd_reso_XP->	DrawClone(i==0 ? "AP" : "P same");
		v_processes[i]->ptge_dd_reso_XP->	SetMarkerSize(7);
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs(OutputFile.c_str());
	fpCanvas-> 						Clear();

	// Resolution's mean
	pr0->ptge_dd_mean_XP->			SetTitle(";Drift time (timebins);Mean (ADC counts/cm)");
	pr0->ptge_dd_mean_XP->			GetXaxis()->SetLimits(0, 510);
	pr0->ptge_dd_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_dd_mean_XP->	DrawClone(i==0 ? "AP" : "P same");
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs(OutputFile.c_str());
	fpCanvas-> 						Clear();

	// dE/dx per phi bins -------------------------------------------------------------------------------------------------------------------
	// Resolution
	pr0->ptge_phi_reso_XP->			SetTitle(";#varphi (#circ);Resolution (%)");
	pr0->ptge_phi_reso_XP->			GetXaxis()->SetLimits(-90, 90);
	pr0->ptge_phi_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_phi_reso_XP->	DrawClone(i==0 ? "AP" : "P same");
		v_processes[i]->ptge_phi_reso_XP->	SetMarkerSize(7);
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs(OutputFile.c_str());
	fpCanvas-> 						Clear();

	// Resolution's mean
	pr0->ptge_phi_mean_XP->			SetTitle(";#varphi (#circ);Mean (ADC counts/cm)");
	pr0->ptge_phi_mean_XP->			GetXaxis()->SetLimits(-90, 90);
	pr0->ptge_phi_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_phi_mean_XP->	DrawClone(i==0 ? "AP" : "P same");
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs(OutputFile.c_str());
	fpCanvas-> 						Clear();

	// dE/dx per theta bins -------------------------------------------------------------------------------------------------------------------
	// Resolution
	pr0->ptge_theta_reso_XP->			SetTitle(";#theta (#circ);Resolution (%)");
	pr0->ptge_theta_reso_XP->			GetXaxis()->SetLimits(-90, 90);
	pr0->ptge_theta_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_theta_reso_XP->	DrawClone(i==0 ? "AP" : "P same");
		v_processes[i]->ptge_theta_reso_XP->	SetMarkerSize(7);
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs(OutputFile.c_str());
	fpCanvas-> 						Clear();

	// Resolution's mean
	pr0->ptge_theta_mean_XP->			SetTitle(";#theta (#circ);Mean (ADC counts/cm)");
	pr0->ptge_theta_mean_XP->			GetXaxis()->SetLimits(-90, 90);
	pr0->ptge_theta_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_theta_mean_XP->	DrawClone(i==0 ? "AP" : "P same");
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs((OutputFile + ")").c_str());
	fpCanvas-> 						Clear();

	for(Process *pr : v_processes) delete pr;
}
