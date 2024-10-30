#include "TrackModel.h"
#include "SetStyle.h"

#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>


void DrawTrackModel(){

	// Set Style
	int lw = 						4;
	TCanvas* pCanvas =				new TCanvas("pCanvas", "pCanvas", 1600, 1200);
	TStyle* ptstyle =				SetMyStyle();
	gROOT->							SetStyle(ptstyle->GetName());
	gStyle->						SetPadRightMargin(0.03);
	gStyle->						SetPadTopMargin(0.08);
	gPad->							UseCurrentStyle();

	int nbins = 					5e3;
	int peaking = 					412;
    float d =                       0;
    float phi =                     42;
    float RC =                      100;
    float drift =                   300;
	TrackModel *p_trackmodel = 		new TrackModel(d, phi, RC, drift, 0, 0);

	gStyle->SetOptStat(0);

	p_trackmodel->ptf1_Signal->SetNpx(nbins);
	p_trackmodel->ptf1_Signal->SetLineWidth(lw);
	p_trackmodel->ptf1_Signal->SetLineColor(kRed);
    p_trackmodel->ptf1_Signal->SetTitle(Form("d = %.0f mm | phi = %.0f#circ | RC = %.0f ns/mm^{2} | drift = %.0f mm", d, phi, RC, drift));
	p_trackmodel->ptf1_Signal->DrawClone();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.06);
    latex.SetTextAlign(13); // top right
    latex.DrawLatex(0.5, 0.8, Form("A_{max} = %.0f ADC counts", p_trackmodel->GetAmax()));
    latex.DrawLatex(0.5, 0.7, Form("T_{max} = %.0f ns", p_trackmodel->Gett_Amax()));
	pCanvas->SaveAs("TrackSignal.pdf(");

    pCanvas->Clear();
    p_trackmodel->ptf1_Charge->SetNpx(nbins);
    p_trackmodel->ptf1_Charge->SetLineWidth(lw);
    p_trackmodel->ptf1_Charge->SetLineColor(kBlue);
    p_trackmodel->ptf1_Charge->Draw();
    pCanvas->SaveAs("TrackSignal.pdf");

    pCanvas->Clear();
    p_trackmodel->ptf1_dETFdt->SetNpx(nbins);
    p_trackmodel->ptf1_dETFdt->SetLineWidth(lw);
    p_trackmodel->ptf1_dETFdt->SetLineColor(kGreen+2);
    p_trackmodel->ptf1_dETFdt->Draw();
    pCanvas->SaveAs("TrackSignal.pdf)");
}




// Constructor
TrackModel::TrackModel(const float &impact, const float &phi, const float &RC, const float &drift, 
                        const float &xpad, const float &ypad){
    // Settings
    fimpact =                   impact;
    fphirad =                   phi/180*M_PI;
    fdrift =                    drift;
    fRC =                       RC;
    fxpad =                     xpad;
    fypad =                     ypad;
    fm =                        std::tan(fphirad);
    fq =                        (cos(fphirad)*fxpad-sin(fphirad)*fxpad+fimpact)/cos(fphirad);

    // Defining the border of the pad considered
    fxleft =					fxpad - fSXWIDTH/2;
    fxright =					fxpad + fSXWIDTH/2;
    fylow =						fypad - fSYWIDTH/2;
    fyhigh =					fypad + fSYWIDTH/2;

    // Eletronics
    ws =                        2/fpeakingTime;
    expfactor =                 ws/(2*Q);
    arg =                       ws/2*std::sqrt(4-1/std::pow(Q,2));

    ptf1_dETFdt =		        new TF1(";time (ns); ADC counts/ns",   this, &TrackModel::dETFdt,  0, ftmax, 0);
    ptf1_Charge =		        new TF1(";time (ns); Charge (fC)",     this, &TrackModel::Charge,  0, ftmax, 0);

    ptf1_Convolution =	        new TF1Convolution(ptf1_Charge, ptf1_dETFdt, 0, ftmax, true);
    ptf1_Convolution->	        SetRange(-ftmax, ftmax);
    ptf1_Convolution->	        SetNofPointsFFT(fnconvolpoints);
    ptf1_Signal =		        new TF1("", *ptf1_Convolution, 0, ftmax, ptf1_Convolution->GetNpar());

    SetAmax(ptf1_Signal->GetMaximum());
    Sett_Amax(ptf1_Signal->GetMaximumX());
}


// Destructor
TrackModel::~TrackModel(){
    delete ptf1_dETFdt;
    delete ptf1_Charge;
    delete ptf1_Convolution;
    delete ptf1_Signal;
}


// Eletronic transfer function
Double_t TrackModel::dETFdt(Double_t *x, Double_t *par){
    Double_t t = x[0];

    if(t <= 0) return 0;
    else return 4096/120*(-ws*std::exp(-ws*t) + std::exp(-expfactor*t) * (expfactor*(std::cos(arg*t) - sinfactor*std::sin(arg*t)) + arg*(std::sin(arg*t) + sinfactor*std::cos(arg*t))))/fSETFMAX;
}




// Charge function
Double_t TrackModel::Charge(Double_t *x, Double_t *par){
    Double_t t = x[0];
    
    if(t <= 0) return 0;
    // Constants
    double sigma = std::sqrt(2 * t / fRC + fDt * fDt * fdrift);  // includes transverse diffusion
    double coeff1 = std::sqrt(2 * (1 + fm * fm) / M_PI) * sigma;
    double denom = 2 * (1 + fm * fm) * sigma * sigma;
    double sqrt2_m2_sigma = std::sqrt(2 * (1 + fm * fm)) * sigma;

    // Terms
    double term11 = std::exp(-std::pow(-fylow + fxright * fm + fq, 2) / denom);
    double term12 = std::exp(-std::pow(-fylow + fxleft * fm + fq, 2) / denom);
    double term13 = std::exp(-std::pow(-fyhigh + fxleft * fm + fq, 2) / denom);
    double term14 = std::exp(-std::pow(-fyhigh + fxright * fm + fq, 2) / denom);

    double term21 = (fylow - fxleft * fm - fq) * std::erf((-fylow + fxleft * fm + fq) / sqrt2_m2_sigma);
    double term22 = (fyhigh - fxleft * fm - fq) * std::erf((-fyhigh + fxleft * fm + fq) / sqrt2_m2_sigma);
    double term23 = (fylow - fxright * fm - fq) * std::erf((-fylow + fxright * fm + fq) / sqrt2_m2_sigma);
    double term24 = (fyhigh - fxright * fm - fq) * std::erf((-fyhigh + fxright * fm + fq) / sqrt2_m2_sigma);

    // Return result
    return flambdaG * std::sqrt(1 + fm * fm) / (2 * fm) * (coeff1 * (term11 - term12 + term13 - term14) + term21 - term22 - term23 + term24);
}
