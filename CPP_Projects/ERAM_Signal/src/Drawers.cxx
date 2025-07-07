#include "TrackModel2.h"
#include "PointModel.h"
#include "SetStyle.h"

#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>

// Draw Punctual deposit model
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
   PointModel *p_pointmodel =
      new PointModel(peaking, 11.28 / 2, 10.19 / 2, 11.28 / 2, 10.19 / 2, 120);

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

   p_pointmodel->ptf1_Charge->GetYaxis()->SetRangeUser(
      0, 1.1 * p_pointmodel->ptf1_Charge->GetMaximum());
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
   p_pointmodel->ptf1_Signal->SetLineColor(kOrange + 7);
   p_pointmodel->ptf1_Signal2->SetNpx(nbins);
   p_pointmodel->ptf1_Signal2->SetLineWidth(lw);
   p_pointmodel->ptf1_Signal2->SetLineStyle(2);
   p_pointmodel->ptf1_Signal2->SetLineColor(kOrange + 4);
   p_pointmodel->ptf1_SignalnoDirac->SetNpx(nbins);
   p_pointmodel->ptf1_SignalnoDirac->SetLineWidth(lw);
   p_pointmodel->ptf1_SignalnoDirac->SetLineColor(kGray + 2);
   p_pointmodel->ptf1_Signal->SetMaximum(
      1.1 * std::max(p_pointmodel->ptf1_Signal->GetMaximum(),
                     p_pointmodel->ptf1_Signal2->GetMaximum()));
   p_pointmodel->ptf1_Signal->SetMinimum(
      1.1 * std::min(p_pointmodel->ptf1_Signal->GetMinimum(),
                     p_pointmodel->ptf1_Signal2->GetMinimum()));
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
      1.1 * std::max(p_pointmodel->ptf1_Signal->GetMaximum(),
                     p_pointmodel->ptf1_Signal2->GetMaximum()));
   p_pointmodel->ptf1_Signal->SetMinimum(
      1.1 * std::min(p_pointmodel->ptf1_Signal->GetMinimum(),
                     p_pointmodel->ptf1_Signal2->GetMinimum()));
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

// Draw Track deposit model
void DrawTrackModel2(double phi, double impact, double drift, double RC)
{
   // Set Style
   int lw = 4;
   TCanvas *pCanvas = new TCanvas("pCanvas", "pCanvas", 1600, 1200);
   TStyle *ptstyle = SetMyStyle();
   gROOT->SetStyle(ptstyle->GetName());
   gStyle->SetPadRightMargin(0.05);
   gStyle->SetPadTopMargin(0.05);
   gPad->UseCurrentStyle();

   std::string outputName =
      Form("Output_PDF/PointSignal_phi%.0f_impact%.0f_drift%.0f_RC%.0f.pdf", phi, impact,
           drift, RC);

   gStyle->SetOptStat(0);
   pCanvas->Divide(3, 3, 0, 0);

   int nbins = 5e3;
   int peaking = 412;
   double xwdith = 11.28;
   double ywdith = 10.19;
   double ymin = 1e6, ymax = 0;
   std::vector<TrackModel2 *> trackModels;
   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
         double tmp_xpad = xwdith / 2 + (i - 1) * xwdith;
         double tmp_ypad = ywdith / 2 + (j - 1) * ywdith;
         TrackModel2 *p_trackmodel =
            new TrackModel2(peaking, tmp_xpad, tmp_ypad, phi, impact, drift, RC);
         std::cout << "Amplitude: " << p_trackmodel->GetAmplitude() << std::endl;
         ymax = std::max(ymax, p_trackmodel->GetAmplitude());
         ymin = std::min(ymin, p_trackmodel->GetMinimum());
         delete p_trackmodel;
      }
   }
   double spacing = 0.1 * (ymax - ymin);
   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
         double tmp_xpad = xwdith / 2 + (i - 1) * xwdith;
         double tmp_ypad = ywdith / 2 + (j - 1) * ywdith;
         TrackModel2 *p_trackmodel =
            new TrackModel2(peaking, tmp_xpad, tmp_ypad, phi, impact, drift, RC);
         trackModels.push_back(p_trackmodel);
         pCanvas->cd(3 * (2 - j) + (i + 1));
         // gPad->SetMargin(0, 0, 0, 0);
         if (j == 0) {
            gPad->SetBottomMargin(0.15);
            p_trackmodel->ptf1_Signal->SetTitle(";time [ns]");
         }
         if (i == 0) {
            gPad->SetLeftMargin(0.15);
            p_trackmodel->ptf1_Signal->SetTitle(";;Signal [ADC]");
         }
         if (j == 0 && i == 0) {
            p_trackmodel->ptf1_Signal->SetTitle(";time [ns];Signal [ADC]");
         }
         p_trackmodel->ptf1_Signal->GetYaxis()->SetRangeUser(ymin - spacing,
                                                             ymax + spacing);
         p_trackmodel->ptf1_Signal->Draw();
      }
   }
   pCanvas->SaveAs(outputName.c_str());
   for (auto &trackModel : trackModels) {
      delete trackModel;
   }
   delete pCanvas;
   delete ptstyle;
}