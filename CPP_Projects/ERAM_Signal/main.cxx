#include <iomanip>

#include "Util.h"
#include "SetStyle.h"

#include "PointModel.h"
#include "TrackModel.h"
#include "GeomFunctions.h"

int main(int argc, char const *argv[])
{
   float length = ComputeLength(0, 0);
   TrackModel trackmodel;
   int nconvpointsnominal = 1e4;
   trackmodel.SetNConvPoints(nconvpointsnominal);
   trackmodel.ComputeRealCharge(500, length, 0, 0, 120, 0, 310 / pow(10, 3.5));
   TF1 signal = *trackmodel.GetSignal();

   int nconvpointsnew = 500;
   trackmodel.SetNConvPoints(nconvpointsnew);
   trackmodel.ComputeRealCharge(500, length, 0, 0, 120, 0, 310 / pow(10, 3.5));
   TF1 signalCoarse = *trackmodel.GetSignal();

   TCanvas c1("c1", "c1", 1800, 1350);
   TStyle *ptstyle = SetMyStyle();
   gROOT->SetStyle(ptstyle->GetName());
   gStyle->SetPadLeftMargin(0.07);
   gStyle->SetPadRightMargin(0.04);
   gStyle->SetPadBottomMargin(0.06);
   gStyle->SetPadTopMargin(0.08);
   gPad->UseCurrentStyle();
   gStyle->SetOptStat(0);

   // Create two pads, one small at the bottom and one large at the top
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad1->SetBottomMargin(0); // Upper and lower pad are joined
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad1->Draw();
   pad2->Draw();

   // Draw signal and signalCoarse in the top pad
   pad1->cd();
   signal.SetLineWidth(2);
   signal.SetLineColor(kRed);
   signal.Draw();
   signalCoarse.SetLineWidth(2);
   signalCoarse.SetLineColor(kBlue);
   signalCoarse.SetLineStyle(2);
   signalCoarse.Draw("SAME");

   // Add legend to the top pad
   TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.85);
   std::string nominalLabel = "Nominal " + std::to_string(nconvpointsnominal);
   legend->AddEntry(&signal, nominalLabel.c_str(), "l");
   legend->AddEntry(&signalCoarse, std::to_string(nconvpointsnew).c_str(), "l");
   legend->Draw();

   // Compute the difference between signal and signalCoarse
   TH1F *hDiff = (TH1F *)signal.GetHistogram()->Clone("hDiff");
   hDiff->Add(signalCoarse.GetHistogram(), -1);

   // Draw the difference in the bottom pad
   pad2->cd();
   hDiff->SetLineColor(kViolet);
   hDiff->Draw();
   c1.SaveAs("Output_PDF/signal.pdf");
}