#include "Draw.h"
#include "Misc_Functions.h"
#include "TFrame.h"
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



void Draw::Run(const Process &pr){
	std::cout << "drawout:" << pr.fdrawfile << std::endl;

	// Plot setup ------------------------------------------------------------------------------------------------------------------------------
	int xmax = 						pr.xmax;
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
	if(pr.fph1f_WF->GetMean() > xmax/2) invX = 0.4;
	pr.fph1f_WF->					Draw("HIST");
	pr.fph1f_XP->					Draw("HIST sames");
	PrintResolution(pr.fph1f_XP, fpCanvas, 0.65-invX, 0.58, kMagenta+2, "XP");
	PrintResolution(pr.fph1f_WF, fpCanvas, 0.65-invX, 0.25, kCyan+2, "WF");
	fpCanvas->						SaveAs((pr.fdrawfile + "(").c_str());

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
		pr.vmod_fph1f_WF[i]->GetMean() > xmax/2 ? invX = 0.4 : invX = 0;
		PrintResolution(pr.vmod_fph1f_WF[i], fpCanvas, 0.65-invX, 0.58, kCyan+2, "WF");
		PrintResolution(pr.vmod_fph1f_XP[i], fpCanvas, 0.65-invX, 0.25, kMagenta+2, "XP");
	}
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	for(int i = 16; i < 32; i++){
		fpCanvas->					cd(i+1-16);
		pr.vmod_fph1f_WF[i]->		Draw("HIST");
		pr.vmod_fph1f_XP[i]->		Draw("HIST same");

		if(pr.vmod_fph1f_WF[i]->GetEntries() < 100) continue;
		xMax = pr.vmod_fph1f_WF[i]->GetXaxis()->GetXmax();
		yMax = pr.vmod_fph1f_WF[i]->GetMaximum();
		pr.vmod_fph1f_WF[i]->GetMean() > xmax/2 ? invX = 0.4 : invX = 0;
		PrintResolution(pr.vmod_fph1f_WF[i], fpCanvas, 0.65-invX, 0.58, kCyan+2, "WF");
		PrintResolution(pr.vmod_fph1f_XP[i], fpCanvas, 0.65-invX, 0.25, kMagenta+2, "XP");
	}
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// dE/dx XP vs WF ---------------------------------------------------------------------------------------------------------------------------
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	gStyle->						SetOptStat(111111);
	gStyle->						SetStatX(0.33);
	gStyle->						SetStatY(0.95);
	pr.fph2f_WFXP->					Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

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
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_dd_mean_WF, 0, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_dd_mean_XP, 0, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_dd_mean_XP->			SetTitle(";Drift time (timebins);Mean (ADC counts/cm)");
	pr.ptge_dd_mean_XP->			GetXaxis()->SetLimits(0, 510);
	pr.ptge_dd_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	pr.ptge_dd_mean_XP->			Draw("APL");
	pr.ptge_dd_mean_WF->			Draw("PL same");
	fpLegend->						Draw();
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

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
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// 2D distribution module per module
	for(int i = 0; i < 32; i++){
		if(i%16 == 0){
			fpCanvas->				Clear();
			fpCanvas->				Divide(4,4);
		}
		fpCanvas->					cd(i%16+1);
		fpCanvas->cd(i%16+1)->		SetRightMargin(0.1);
		pr.vmod_fph2f_XPtmean[i]->	Draw("colz");
		if(i%16 == 15)fpCanvas->	SaveAs(pr.fdrawfile.c_str());
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
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_phi_mean_WF, 0, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_phi_mean_XP, 0, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_phi_mean_XP->			SetTitle(";#varphi (#circ);Mean (ADC counts/cm)");
	pr.ptge_phi_mean_XP->			GetXaxis()->SetLimits(-90, 90);
	pr.ptge_phi_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	pr.ptge_phi_mean_XP->			Draw("APL");
	pr.ptge_phi_mean_WF->			Draw("PL same");
	fpLegend->						Draw();
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// 2D distribution
	fpCanvas->						Clear();
	gPad->							SetTopMargin(0.02);
	gPad->							SetRightMargin(0.11);
	gStyle->						SetOptStat(0);
	gStyle->						SetOptFit(0);
	pr.fph2f_XPphi->				SetTitle(";#varphi (#circ);dE/dx with XP (ADC counts/cm)");
	pr.fph2f_XPphi->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

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
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_theta_mean_WF, 0, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_theta_mean_XP, 0, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_theta_mean_XP->			SetTitle(";#theta (#circ);Mean (ADC counts/cm)");
	pr.ptge_theta_mean_XP->			GetXaxis()->SetLimits(-90, 90);
	pr.ptge_theta_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	pr.ptge_theta_mean_XP->			Draw("APL");
	pr.ptge_theta_mean_WF->			Draw("PL same");
	fpLegend->						Draw();
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// 2D distribution
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.11);
	pr.fph2f_XPtheta->				SetTitle(";#theta (#circ);dE/dx with XP (ADC counts/cm)");
	pr.fph2f_XPtheta->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// dEdx per momentum bin --------------------------------------------------------------------------------------------------------------------
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
	fpCanvas->						SaveAs((pr.fdrawfile + "(").c_str());

	// Mean
	fpCanvas->						Clear();
	Graphic_setup(pr.ptge_mom_mean_WF, 0, 33, kCyan+2, 2, kCyan+2);
	Graphic_setup(pr.ptge_mom_mean_XP, 0, 47, kMagenta+2, 2, kMagenta+2);
	pr.ptge_mom_mean_XP->			SetTitle(";Momentum (MeV/c);Mean (ADC counts/cm)");
	pr.ptge_mom_mean_XP->			GetXaxis()->SetLimits(-momrange, momrange);
	pr.ptge_mom_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	pr.ptge_mom_mean_XP->			Draw("APL");
	pr.ptge_mom_mean_WF->			Draw("PL same");
	fpLegend->						Draw();
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// 2D distribution
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.11);
	pr.fph2f_XPmom->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph2f_WFmom->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

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
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());
	delete plegmom;

	// Geometry plots ----------------------------------------------------------------------------------------------------------------------------

	// dE/dx vs track length
	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	gPad->							SetTopMargin(0.02);
	gStyle->						SetStatX(0.87);
	pr.fph2f_XPlen->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	gPad->							SetLogz();
	pr.fph2f_XZ->					Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());
	gPad->							SetLogz(0);

	fpCanvas->						Clear();
	pr.fph2f_phitheta->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.02);
	pr.fph1f_phi->					Draw("HIST");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph1f_theta->				Draw("HIST");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	gStyle->						SetOptStat(0);
	pr.fph1f_trklen->				Draw("HIST");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph1f_chi2->					Draw("HIST");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	gPad->							SetRightMargin(0.13);
	pr.fph2f_chi2mom->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// Debugging ------------------------------------------------------------------------------------------------------------------------------
	fpCanvas->						Clear();
	pr.fph2f_pullemu->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph2f_chi2ndfphi->			Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph2f_momtheta->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph2f_momR->					Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph2f_chi2ndfR->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph2f_lentheta->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	fpCanvas->						Clear();
	pr.fph2f_lenphi->				Draw("colz");
	fpCanvas->						SaveAs(pr.fdrawfile.c_str());

	// dE/dx per momentum bin --------------------------------------------------------------------------------------------------------------------
	gPad->							SetTopMargin(0.08);
	gPad->							SetRightMargin(0.02);
	for(int i=0;i<nmombins;i++){
		fpCanvas->					Clear();
		pr.vmom_fph1f_WF[i]->				SetAxisRange(0, 1.1	*std::max({pr.vmom_fph1f_WF[i]->GetMaximum(), pr.vmom_fph1f_XP[i]->GetMaximum()}),	"Y");
		float invX = 				0;
		if(pr.vmom_fph1f_WF[i]->GetMean() > xmax/2) invX = 0.4;
		Graphic_setup(pr.vmom_fph1f_WF[i], 0.5, 1, kCyan+1, 2, kCyan-2, kCyan, 0.2);
		Graphic_setup(pr.vmom_fph1f_XP[i], 0.5, 1, kMagenta+1, 2, kMagenta-2, kMagenta, 0.2);
		pr.vmom_fph1f_WF[i]->				Draw("HIST");
		pr.vmom_fph1f_XP[i]->				Draw("HIST same");
		PrintResolution(pr.vmom_fph1f_XP[i], fpCanvas, 0.65-invX, 0.58, kMagenta+2, "XP");
		PrintResolution(pr.vmom_fph1f_WF[i], fpCanvas, 0.65-invX, 0.25, kCyan+2, "WF");
		fpCanvas->					SaveAs(pr.fdrawfile.c_str());
	}

	// dE/dx per theta bin --------------------------------------------------------------------------------------------------------------------
	gPad->							SetTopMargin(0.08);
	gPad->							SetRightMargin(0.02);
	for(int i=0;i<nthetabins;i++){
		fpCanvas->					Clear();
		pr.vtheta_fph1f_WF[i]->				SetAxisRange(0, 1.1	*std::max({pr.vtheta_fph1f_WF[i]->GetMaximum(), pr.vtheta_fph1f_XP[i]->GetMaximum()}),	"Y");
		float invX = 				0;
		if(pr.vtheta_fph1f_WF[i]->GetMean() > xmax/2) invX = 0.4;
		Graphic_setup(pr.vtheta_fph1f_WF[i], 0.5, 1, kCyan+1, 2, kCyan-2, kCyan, 0.2);
		Graphic_setup(pr.vtheta_fph1f_XP[i], 0.5, 1, kMagenta+1, 2, kMagenta-2, kMagenta, 0.2);
		pr.vtheta_fph1f_WF[i]->				Draw("");
		pr.vtheta_fph1f_XP[i]->				Draw("same");
		PrintResolution(pr.vtheta_fph1f_XP[i], fpCanvas, 0.65-invX, 0.58, kMagenta+2, "XP");
		PrintResolution(pr.vtheta_fph1f_WF[i], fpCanvas, 0.65-invX, 0.25, kCyan+2, "WF");
		fpCanvas->					SaveAs(pr.fdrawfile.c_str());
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
	fpCanvas->					SaveAs(pr.fdrawfile.c_str());


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
	fpCanvas->						SaveAs((pr.fdrawfile + ")").c_str());
	fpCanvas->						Clear();

	// Reset style
	gPad->							UseCurrentStyle();
}


void Draw::CompareRuns(const std::vector<Process*> &v_processes, const std::string &compfolderpath, const std::string &runsettagname){

	Process pr =					*v_processes[0];
	int ncomparisons =				(int)v_processes.size();
	std::string OutputFile =		compfolderpath + "/Comparison_" + runsettagname + pr.fcomment;
	for(int i=1;i<ncomparisons;i++) OutputFile += "_VS" + v_processes[i]->fcomment;
	OutputFile +=					".pdf";


	// Plot setup ------------------------------------------------------------------------------------------------------------------------------
	int xmax = 						pr.xmax;
	int momrange = 					pr.momrange;
	int nmombins = 					pr.nmombins;
	fpCanvas->						cd();
	gStyle->						SetOptStat(0);
	gStyle->						SetOptFit(0);
	gPad->							SetTopMargin(0.05);
	gPad->							SetRightMargin(0.03);
	float invX = 					0;

	// Legend
	fpLegend =						new TLegend(0.55+0.03*ncomparisons,0.83-0.03*ncomparisons,0.88,0.95);
	fpLegend->						SetTextSize(0.06-0.005*ncomparisons);
	fpLegend->						SetFillStyle(0);
	fpLegend->						SetTextColor(kBlue-1);

	float ampmaxXP = 0, ampmaxWF = 0;
	for(int i=0;i<ncomparisons;i++){
		Process pr_tmp = *v_processes[i];
		Graphic_setup(pr_tmp.fph1f_WF, 0.5, 1, colors[i], 2, colors[i], colors[i]-2, 0.2);
		Graphic_setup(pr_tmp.fph1f_XP, 0.5, 1, colors[i], 2, colors[i], colors[i]-2, 0.2);
		if(pr_tmp.fph1f_WF ->GetMaximum() > ampmaxWF) ampmaxWF = pr_tmp.fph1f_WF ->GetMaximum();
		if(pr_tmp.fph1f_XP ->GetMaximum() > ampmaxXP) ampmaxXP = pr_tmp.fph1f_XP ->GetMaximum();
		for(int j=0;j<32;j++){
			Graphic_setup(pr_tmp.vmod_fph1f_WF[j], 0.5, 1, colors[i], 1, colors[i], colors[i]-2, 0.2);
			Graphic_setup(pr_tmp.vmod_fph1f_XP[j], 0.5, 1, colors[i], 1, colors[i], colors[i]-2, 0.2);
		}
		Graphic_setup(pr_tmp.ptge_mom_reso_WF, 2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp.ptge_mom_reso_XP, 2, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp.ptge_mom_mean_WF, 0, markers[i], colors[i], 2, colors[i]);
		Graphic_setup(pr_tmp.ptge_mom_mean_XP, 0, markers[i], colors[i], 2, colors[i]);
	}

	// Global dE/dx plot -------------------------------------------------------------------------------------------------------------------
	pr.fph1f_XP->	SetAxisRange(0, 1.1	*ampmaxXP,	"Y");
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->fph1f_XP->	Draw(i==0 ? "HIST" : "HIST same");
		float ywidth = 0.4 - 0.05*ncomparisons;
		float ylow = 0.62 - (ywidth+0.01)*i;
		float xlow = 0.56+0.04*ncomparisons-invX;
		PrintResolution(v_processes[i]->fph1f_XP, fpCanvas, xlow, ylow, ywidth, colors[i], v_processes[i]->fcomment.c_str());
		std::cout << xlow << " " << ylow << " " << ywidth << std::endl;
	}
	fpCanvas->						SaveAs((OutputFile + "(").c_str());
	fpCanvas->						Clear();

	// dE/dx vs momentum -------------------------------------------------------------------------------------------------------------------
	// Resolution
	pr.ptge_mom_reso_XP->			SetTitle(";Momentum (MeV/c);Resolution (%)");
	pr.ptge_mom_reso_XP->			GetXaxis()->SetLimits(-momrange, momrange);
	pr.ptge_mom_reso_XP->			GetYaxis()->SetRangeUser(resomin, resomax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_mom_reso_XP->	DrawClone(i==0 ? "AP" : "P same");
		v_processes[i]->ptge_mom_reso_XP->	SetMarkerSize(7);
		fpLegend->AddEntry(v_processes[i]->ptge_mom_reso_XP,  v_processes[i]->fcomment.c_str(), "p");
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs(OutputFile.c_str());
	fpCanvas-> 						Clear();	

	// Resolution's mean
	pr.ptge_mom_mean_XP->			SetTitle(";Momentum (MeV/c);Mean (ADC counts/cm)");
	pr.ptge_mom_mean_XP->			GetXaxis()->SetLimits(-momrange, momrange);
	pr.ptge_mom_mean_XP->			GetYaxis()->SetRangeUser(meanmin, meanmax);
	for(int i=0;i<ncomparisons;i++){
		v_processes[i]->ptge_mom_mean_XP->	DrawClone(i==0 ? "APL" : "PL same");
	}
	fpLegend->						Draw();
	fpCanvas->						SaveAs((OutputFile + ")").c_str());
	fpCanvas-> 						Clear();
}
