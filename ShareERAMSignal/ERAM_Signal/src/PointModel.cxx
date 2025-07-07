#include "PointModel.h"
#include "SetStyle.h"

#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>

void DrawPointModel()
{
   // Set Style
   int lw = 4;
   TCanvas *pCanvas = new TCanvas("pCanvas", "pCanvas", 1600, 1200);
   TStyle *ptstyle = SetMyStyle();
   gROOT->SetStyle(ptstyle->GetName());
   gStyle->SetPadRightMargin(0.05);
   gStyle->SetPadTopMargin(0.05);
   gPad->UseCurrentStyle();

   int nbins = 5e3;
   int peaking = 412;
   PointModel *p_pointmodel = new PointModel(peaking, 11.28 / 2, 10.19 / 2, 11.28 / 2, 10.19 / 2, 120);

   std::string outputName = "Output_PDF/PointSignal_" + std::to_string(peaking) + ".pdf";

   p_pointmodel->ptf1_Charge->SetNpx(nbins);
   p_pointmodel->ptf1_Charge->SetLineWidth(lw);
   p_pointmodel->ptf1_Charge->SetLineColor(kRed);
   p_pointmodel->ptf1_CurrentScaled->SetNpx(nbins);
   p_pointmodel->ptf1_CurrentScaled->SetLineWidth(lw);
   p_pointmodel->ptf1_CurrentScaled->SetLineColor(kBlue);
   p_pointmodel->ptf1_ETF->SetNpx(nbins);
   p_pointmodel->ptf1_ETF->SetLineWidth(lw);
   p_pointmodel->ptf1_ETF->SetLineColor(kBlue);
   p_pointmodel->ptf1_dETFdt->SetNpx(nbins);
   p_pointmodel->ptf1_dETFdt->SetLineWidth(lw);
   p_pointmodel->ptf1_dETFdt->SetLineColor(kRed);
   gStyle->SetOptStat(0);

   p_pointmodel->ptf1_Charge->GetYaxis()->SetRangeUser(0, 1.1 * p_pointmodel->ptf1_Charge->GetMaximum());
   p_pointmodel->ptf1_Charge->Draw();
   pCanvas->SaveAs((outputName + "(").c_str());

   pCanvas->Clear();
   p_pointmodel->ptf1_CurrentScaled->Draw();
   pCanvas->SaveAs(outputName.c_str());

   pCanvas->Clear();
   p_pointmodel->ptf1_ETF->Draw();
   pCanvas->SaveAs(outputName.c_str());

   pCanvas->Clear();
   p_pointmodel->ptf1_dETFdt->Draw();
   pCanvas->SaveAs(outputName.c_str());

   pCanvas->Clear();
   p_pointmodel->ptf1_Signal->SetNpx(nbins);
   p_pointmodel->ptf1_Signal->SetLineWidth(lw);
   p_pointmodel->ptf1_Signal->SetLineColor(kOrange+7);
   p_pointmodel->ptf1_Signal2->SetNpx(nbins);
   p_pointmodel->ptf1_Signal2->SetLineWidth(lw);
   p_pointmodel->ptf1_Signal2->SetLineStyle(2);
   p_pointmodel->ptf1_Signal2->SetLineColor(kOrange+4);
   p_pointmodel->ptf1_SignalnoDirac->SetNpx(nbins);
   p_pointmodel->ptf1_SignalnoDirac->SetLineWidth(lw);
   p_pointmodel->ptf1_SignalnoDirac->SetLineColor(kGray + 2);
   p_pointmodel->ptf1_Signal->SetMaximum(
      1.1 * std::max(p_pointmodel->ptf1_Signal->GetMaximum(), p_pointmodel->ptf1_Signal2->GetMaximum()));
   p_pointmodel->ptf1_Signal->SetMinimum(
      1.1 * std::min(p_pointmodel->ptf1_Signal->GetMinimum(), p_pointmodel->ptf1_Signal2->GetMinimum()));
   p_pointmodel->ptf1_Signal->DrawClone();
   p_pointmodel->ptf1_Signal2->Draw("same");
   p_pointmodel->ptf1_SignalnoDirac->Draw("same");
   TLegend *plegend = new TLegend(0.5, 0.6, 0.9, 0.9);
   p_pointmodel->ptf1_Signal->SetLineStyle(2);
   plegend->AddEntry(p_pointmodel->ptf1_Signal, "I^{H} #otimes T^{H}", "l");
   plegend->AddEntry(p_pointmodel->ptf1_SignalnoDirac, "I #otimes T^{H}", "l");
   plegend->AddEntry(p_pointmodel->ptf1_Signal2, "Q^{H} #otimes #partial_{t}T^{H}", "l");
   plegend->Draw();
   pCanvas->SaveAs(outputName.c_str());
   
   pCanvas->Clear();
   p_pointmodel->ptf1_Signal->SetNpx(nbins);
   p_pointmodel->ptf1_Signal->SetLineWidth(lw);
   p_pointmodel->ptf1_Signal->SetLineStyle(1);
   p_pointmodel->ptf1_Signal->SetLineColor(kOrange + 7);
   p_pointmodel->ptf1_SignalnoDirac->SetNpx(nbins);
   p_pointmodel->ptf1_SignalnoDirac->SetLineWidth(lw);
   p_pointmodel->ptf1_SignalnoDirac->SetLineStyle(2);
   p_pointmodel->ptf1_SignalnoDirac->SetLineColor(kGray + 2);
   p_pointmodel->ptf1_ETFscaled->SetNpx(nbins);
   p_pointmodel->ptf1_ETFscaled->SetLineWidth(lw);
   p_pointmodel->ptf1_ETFscaled->SetLineStyle(2);
   p_pointmodel->ptf1_ETFscaled->SetLineColor(kBlue);
   p_pointmodel->ptf1_Signal->SetMaximum(
      1.1 * std::max(p_pointmodel->ptf1_Signal->GetMaximum(), p_pointmodel->ptf1_Signal2->GetMaximum()));
   p_pointmodel->ptf1_Signal->SetMinimum(
      1.1 * std::min(p_pointmodel->ptf1_Signal->GetMinimum(), p_pointmodel->ptf1_Signal2->GetMinimum()));
   p_pointmodel->ptf1_Signal->Draw();
   p_pointmodel->ptf1_SignalnoDirac->Draw("same");
   p_pointmodel->ptf1_ETFscaled->Draw("same");
   delete plegend;
   plegend = new TLegend(0.45, 0.6, 0.9, 0.9);
   plegend->SetTextSize(0.075);
   plegend->AddEntry(p_pointmodel->ptf1_Signal, "Complete signal", "l");
   plegend->AddEntry(p_pointmodel->ptf1_ETFscaled, "Pulse signal", "l");
   plegend->AddEntry(p_pointmodel->ptf1_SignalnoDirac, "Escape signal", "l");
   plegend->Draw();
   pCanvas->SaveAs((outputName + ")").c_str());
}

// Constructor
PointModel::PointModel(float peakingTime, float fxpad, float fypad, float xcharge, float ycharge, float fRC,
                       float fQdrop)
{
   SetPeakingTime(peakingTime);
   SetPadParam(fxpad, fypad, fRC);
   SetChargeParam(xcharge, ycharge, fQdrop);
   std::cout << ws << " " << expfactor << " " << arg << std::endl;

   // Defining the border of the pad consideBlue
   fxleft = fxpad - fSXWIDTH / 2;
   fxright = fxpad + fSXWIDTH / 2;
   fylow = fypad - fSYWIDTH / 2;
   fyhigh = fypad + fSYWIDTH / 2;

   ptf1_ETF = new TF1(";time (ns);ADC counts", this, &PointModel::ETF, 0, ftmax, 0);
   ptf1_ETFscaled = new TF1(";time (ns);ADC counts", this, &PointModel::ETFscaled, 0, ftmax, 0);
   ptf1_dETFdt = new TF1(";time (ns); ADC counts/ns", this, &PointModel::dETFdt, 0, ftmax, 0);
   ptf1_Charge = new TF1(";time (ns); Charge (fC)", this, &PointModel::Charge, 0, ftmax, 0);
   ptf1_Current = new TF1(";time (ns);Current (#muA)", this, &PointModel::Current, 0, ftmax, 0);
   ptf1_CurrentScaled = new TF1(";time (ns);Current (nA)", this, &PointModel::CurrentScaled, 0, ftmax, 0);
   ptf1_CurrentnoDirac = new TF1(";time (ns);Current (#muA)", this, &PointModel::CurrentnoDirac, 0, ftmax, 0);

   ptf1_Convolution = new TF1Convolution(ptf1_Current, ptf1_ETF, 0, ftmax, true);
   ptf1_Convolution->SetRange(-ftmax, ftmax);
   ptf1_Convolution->SetNofPointsFFT(fnconvolpoints);
   ptf1_Signal = new TF1(";time (ns);ADC counts", *ptf1_Convolution, 0, ftmax, ptf1_Convolution->GetNpar());

   ptf1_ConvolutionnoDirac = new TF1Convolution(ptf1_CurrentnoDirac, ptf1_ETF, 0, ftmax, true);
   ptf1_ConvolutionnoDirac->SetRange(-ftmax, ftmax);
   ptf1_ConvolutionnoDirac->SetNofPointsFFT(fnconvolpoints);
   ptf1_SignalnoDirac = new TF1("", *ptf1_ConvolutionnoDirac, 0, ftmax, ptf1_Convolution->GetNpar());

   ptf1_Convolution2 = new TF1Convolution(ptf1_Charge, ptf1_dETFdt, 0, ftmax, true);
   ptf1_Convolution2->SetRange(-ftmax, ftmax);
   ptf1_Signal2 = new TF1("", *ptf1_Convolution2, 0, ftmax, ptf1_Convolution2->GetNpar());

   ptf1_Convolution3 = new TF1Convolution(ptf1_Charge, ptf1_ETF, 0, ftmax, true);
   ptf1_Convolution3->SetRange(-ftmax, ftmax);
   ptf1_Signal3 = new TF1("", *ptf1_Convolution3, 0, ftmax, ptf1_Convolution3->GetNpar());
}

// Destructor
PointModel::~PointModel()
{
   delete ptf1_ETF;
   delete ptf1_ETFscaled;
   delete ptf1_dETFdt;
   delete ptf1_Charge;
   delete ptf1_Current;
   delete ptf1_CurrentScaled;
   delete ptf1_CurrentnoDirac;
   delete ptf1_Convolution;
   delete ptf1_ConvolutionnoDirac;
   delete ptf1_Convolution2;
   delete ptf1_Convolution3;
   delete ptf1_Signal;
}

// Settings
void PointModel::SetPeakingTime(float peakingTime)
{
   fpeakingTime = peakingTime;
   ws = 2 / fpeakingTime;
   expfactor = ws / (2 * Q);
   arg = ws / 2 * std::sqrt(4 - 1 / std::pow(Q, 2));
}

// Eletronic transfer function
Double_t PointModel::ETF(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t <= 0)
      return 0;
   else
      return 4096 / 120 *
             (std::exp(-ws * t) + std::exp(-expfactor * t) * (sinfactor * std::sin(arg * t) - std::cos(arg * t))) /
             fSETFMAX;
}
Double_t PointModel::ETFscaled(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t <= 0)
      return 0;
   else
      return Q0 * 4096 / 120 *
             (std::exp(-ws * t) + std::exp(-expfactor * t) * (sinfactor * std::sin(arg * t) - std::cos(arg * t))) /
             fSETFMAX;
}

Double_t PointModel::dETFdt(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t <= 0)
      return 0;
   else
      return 4096 / 120 *
             (-ws * std::exp(-ws * t) +
              std::exp(-expfactor * t) * (expfactor * (std::cos(arg * t) - sinfactor * std::sin(arg * t)) +
                                          arg * (std::sin(arg * t) + sinfactor * std::cos(arg * t)))) /
             fSETFMAX;
}

// Charge function
Double_t PointModel::Charge(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t <= 0)
      return 0;

   float sigma = std::sqrt(2 * t / fRC);
   auto erfRC = [sigma](float edge, float charge) { return std::erf((edge - charge) / (std::sqrt(2) * sigma)); };

   return fQdrop / 4 * (erfRC(fyhigh, fycharge) - erfRC(fylow, fycharge)) *
          (erfRC(fxright, fxcharge) - erfRC(fxleft, fxcharge));
}

// Current function
Double_t PointModel::Current(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t < 0)
      return 0;
   if (t <= 1 / fnconvolpoints)
      return fQdrop / (2 * ftmax) * fnconvolpoints; // normalized to mimic a dirac delta function

   float sigma = std::sqrt(2 * t / fRC);
   float norm = fQdrop / (4 * std::sqrt(2 * TMath::Pi()) * sigma * t);
   auto erfRC = [sigma](float edge, float charge) { return std::erf((edge - charge) / (std::sqrt(2) * sigma)); };
   auto expRC = [sigma, t](float edge, float charge) {
      return std::exp(-std::pow(edge - charge, 2) / (2 * std::pow(sigma, 2))) * (edge - charge);
   };

   return norm *
          ((erfRC(fyhigh, fycharge) - erfRC(fylow, fycharge)) * (expRC(fxleft, fxcharge) - expRC(fxright, fxcharge)) +
           (erfRC(fxright, fxcharge) - erfRC(fxleft, fxcharge)) * (expRC(fylow, fycharge) - expRC(fyhigh, fycharge)));
}

Double_t PointModel::CurrentScaled(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t < 0)
      return 0;
   if (t <= 1 / fnconvolpoints)
      return fQdrop / (2 * ftmax) * fnconvolpoints; // normalized to mimic a dirac delta function

   float sigma = std::sqrt(2 * t / fRC);
   float norm = fQdrop / (4 * std::sqrt(2 * TMath::Pi()) * sigma * t);
   auto erfRC = [sigma](float edge, float charge) { return std::erf((edge - charge) / (std::sqrt(2) * sigma)); };
   auto expRC = [sigma, t](float edge, float charge) {
      return std::exp(-std::pow(edge - charge, 2) / (2 * std::pow(sigma, 2))) * (edge - charge);
   };

   return 1e3 * norm *
          ((erfRC(fyhigh, fycharge) - erfRC(fylow, fycharge)) * (expRC(fxleft, fxcharge) - expRC(fxright, fxcharge)) +
           (erfRC(fxright, fxcharge) - erfRC(fxleft, fxcharge)) * (expRC(fylow, fycharge) - expRC(fyhigh, fycharge)));
}

Double_t PointModel::CurrentnoDirac(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t <= 0)
      return 0;

   float sigma = std::sqrt(2 * t / fRC);
   float norm = fQdrop / (4 * std::sqrt(2 * TMath::Pi()) * sigma * t);
   auto erfRC = [sigma](float edge, float charge) { return std::erf((edge - charge) / (std::sqrt(2) * sigma)); };
   auto expRC = [sigma, t](float edge, float charge) {
      return std::exp(-std::pow(edge - charge, 2) / (2 * std::pow(sigma, 2))) * (edge - charge);
   };

   return norm *
          ((erfRC(fyhigh, fycharge) - erfRC(fylow, fycharge)) * (expRC(fxleft, fxcharge) - expRC(fxright, fxcharge)) +
           (erfRC(fxright, fxcharge) - erfRC(fxleft, fxcharge)) * (expRC(fylow, fycharge) - expRC(fyhigh, fycharge)));
}