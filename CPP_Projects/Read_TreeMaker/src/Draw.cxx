#include "Draw.h"
#include "Misc_Functions.h"
#include "TFrame.h"
#include "TPaveText.h"
#include "SetStyle.h"

#include <fstream>

// Constructor
Draw::Draw()
{
   fpStyle = SetMyStyle();
   gROOT->SetStyle(fpStyle->GetName());
   fpCanvas = new TCanvas("fpCanvas", "fpCanvas", 1800, 1350);
   gPad->UseCurrentStyle();
}

// Destructor
Draw::~Draw() {}

void Draw::Run(const std::string &filepath)
{
   TFile inputFile(filepath.c_str(), "READ");
   Process pr = *(Process *)inputFile.Get("Process");
   fRealpathPDF = foutputPDFFolder + "/" + pr.ffileName + ".pdf";
   std::cout << "drawout:" << fRealpathPDF << std::endl;

   // -----------------------------------------------------------------------------------------------------------------
   // Plot setup
   int dEdxmax = pr.dEdxmax;
   int absmomrange = pr.absmomrange;
   int momrange = pr.momrange;
   int nabsmombins = pr.nabsmombins;
   int nmombins = pr.nmombins;
   int nthetabins = pr.nthetabins;
   fpCanvas->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   fpLegend = new TLegend(0.15, 0.74, 0.35, 0.94);
   fpLegend->SetTextSize(0.08);
   fpLegend->SetFillStyle(0);
   fpLegend->SetTextColor(kBlack);
   TLegendEntry *pentryWF = fpLegend->AddEntry((TObject *)0, "WF", "p");
   pentryWF->SetMarkerColor(kCyan + 2);
   pentryWF->SetMarkerSize(9);
   pentryWF->SetMarkerStyle(33);
   TLegendEntry *pentryXP = fpLegend->AddEntry((TObject *)0, "XP", "p");
   pentryXP->SetMarkerColor(kMagenta + 2);
   pentryXP->SetMarkerSize(9);
   pentryXP->SetMarkerStyle(47);
   gPad->SetTopMargin(0.05);
   float invX = 0;

   // -----------------------------------------------------------------------------------------------------------------
   // Global dE/dx plot
   Graphic_setup(pr.fph1f_WF, 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan, 0.2);
   Graphic_setup(pr.fph1f_XP, 0.5, 1, kMagenta + 2, 2, kMagenta - 2, kMagenta, 0.2);
   pr.fph1f_WF->SetAxisRange(
      0, 1.1 * std::max({pr.fph1f_WF->GetMaximum(), pr.fph1f_XP->GetMaximum()}), "Y");
   if (pr.fph1f_WF->GetMean() > dEdxmax / 2)
      invX = 0.4;
   pr.fph1f_WF->Draw("HIST");
   pr.fph1f_XP->Draw("HIST sames");
   PrintResolution(pr.fph1f_XP, fpCanvas, 0.65 - invX, 0.58, kMagenta + 2, "XP");
   PrintResolution(pr.fph1f_WF, fpCanvas, 0.65 - invX, 0.25, kCyan + 2, "WF");
   fpCanvas->SaveAs((fRealpathPDF + "(").c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx ERAM by ERAM
   float maxdEdx = 0;
   for (TH1F *hist : pr.vmod_fph1f_WF)
      if (hist->GetMaximum() > maxdEdx)
         maxdEdx = hist->GetMaximum();
   for (TH1F *hist : pr.vmod_fph1f_XP)
      if (hist->GetMaximum() > maxdEdx)
         maxdEdx = hist->GetMaximum();

   float xMax = 0, yMax = 0, reso = 0, dreso = 0;
   fpCanvas->Clear();
   fpCanvas->Divide(4, 4);
   for (int i = 0; i < 32; i++) {
      pr.vmod_fph1f_WF[i]->SetAxisRange(0, 1.1 * maxdEdx, "Y");
      pr.vmod_fph1f_XP[i]->SetAxisRange(0, 1.1 * maxdEdx, "Y");
      Graphic_setup(pr.vmod_fph1f_WF[i], 0.5, 1, kCyan + 1, 1, kCyan - 2, kCyan, 0.2);
      Graphic_setup(pr.vmod_fph1f_XP[i], 0.5, 1, kMagenta + 2, 1, kMagenta - 2, kMagenta,
                    0.2);
   }
   for (int i = 0; i < 16; i++) {
      fpCanvas->cd(i + 1);
      pr.vmod_fph1f_WF[i]->Draw("HIST");
      pr.vmod_fph1f_XP[i]->Draw("HIST same");
      if (pr.vmod_fph1f_WF[i]->GetEntries() < 100)
         continue;
      xMax = pr.vmod_fph1f_WF[i]->GetXaxis()->GetXmax();
      yMax = pr.vmod_fph1f_WF[i]->GetMaximum();
      pr.vmod_fph1f_WF[i]->GetMean() > dEdxmax / 2 ? invX = 0.4 : invX = 0;
      PrintResolution(pr.vmod_fph1f_XP[i], fpCanvas, 0.65 - invX, 0.58, kMagenta + 2,
                      "XP");
      PrintResolution(pr.vmod_fph1f_WF[i], fpCanvas, 0.65 - invX, 0.25, kCyan + 2, "WF");
   }
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   for (int i = 16; i < 32; i++) {
      fpCanvas->cd(i + 1 - 16);
      pr.vmod_fph1f_WF[i]->Draw("HIST");
      pr.vmod_fph1f_XP[i]->Draw("HIST same");

      if (pr.vmod_fph1f_WF[i]->GetEntries() < 100)
         continue;
      xMax = pr.vmod_fph1f_WF[i]->GetXaxis()->GetXmax();
      yMax = pr.vmod_fph1f_WF[i]->GetMaximum();
      pr.vmod_fph1f_WF[i]->GetMean() > dEdxmax / 2 ? invX = 0.4 : invX = 0;
      PrintResolution(pr.vmod_fph1f_XP[i], fpCanvas, 0.65 - invX, 0.58, kMagenta + 2,
                      "XP");
      PrintResolution(pr.vmod_fph1f_WF[i], fpCanvas, 0.65 - invX, 0.25, kCyan + 2, "WF");
   }
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx XP vs WF
   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   gStyle->SetOptStat("merou");
   gStyle->SetStatX(0.33);
   gStyle->SetStatY(0.95);
   pr.fph2f_WFXP->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs X position
   // Resolution
   fpCanvas->Clear();
   fpCanvas->SetGrid(1, 1);
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_X_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_X_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_X_reso_XP->SetTitle(";Track X position (mm);dE/dx resolution (%)");
   pr.ptge_X_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_X_reso_XP->DrawClone("AP");
   pr.ptge_X_reso_WF->DrawClone("P same");
   pr.ptge_X_reso_XP->SetMarkerSize(7);
   pr.ptge_X_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   fpCanvas->SetGrid(1, 1);
   Graphic_setup(pr.ptge_X_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_X_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_X_mean_XP->SetTitle(";Track X position (mm);dE/dx (ADC counts/cm)");
   pr.ptge_X_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_X_mean_XP->Draw("AP");
   pr.ptge_X_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   gPad->SetTopMargin(0.02);
   gPad->SetRightMargin(0.1);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   // 2D distribution WF vs X
   pr.fph2f_WFX->Draw("colz");
   TGraphErrors *ptge_mom_WFX = Convert_TH2_TGE(pr.fph2f_WFX);
   TF1 linearFitWFX("linearFitWFX", "pol1", 0, 1000);
   ptge_mom_WFX->Fit(&linearFitWFX, "RQ");
   linearFitWFX.SetLineColor(kRed);
   linearFitWFX.Draw("same");
   TPaveText pavetext(0.5, 0.8, 0.88, 0.95, "NDC");
   pavetext.AddText(Form("y = %.2fx + %.2f", linearFitWFX.GetParameter(1),
                         linearFitWFX.GetParameter(0)));
   pavetext.SetFillColorAlpha(kWhite, 0.9);
   pavetext.SetBorderSize(1);
   pavetext.SetLineColor(kBlue - 1);
   pavetext.SetTextColor(kRed);
   pavetext.SetTextSize(0.07);
   pavetext.DrawClone();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution XP vs X
   pr.fph2f_XPX->Draw("colz");
   TGraphErrors *ptge_mom_XPX = Convert_TH2_TGE(pr.fph2f_XPX);
   TF1 linearFitXPX("linearFitXPX", "pol1", 0, 1000);
   ptge_mom_XPX->Fit(&linearFitXPX, "RQ");
   linearFitXPX.SetLineColor(kRed);
   linearFitXPX.Draw("same");
   pavetext.Clear();
   pavetext.AddText(Form("y = %.2fx + %.2f", linearFitXPX.GetParameter(1),
                         linearFitXPX.GetParameter(0)));
   pavetext.DrawClone();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs time bin
   // Resolution
   fpCanvas->Clear();
   fpCanvas->SetGrid(1, 1);
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_dt_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_dt_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_dt_reso_XP->SetTitle(";Drift time (timebins);dE/dx resolution (%)");
   pr.ptge_dt_reso_XP->GetXaxis()->SetLimits(0, 510);
   pr.ptge_dt_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_dt_reso_XP->DrawClone("AP");
   pr.ptge_dt_reso_WF->DrawClone("P same");
   pr.ptge_dt_reso_XP->SetMarkerSize(7);
   pr.ptge_dt_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   fpCanvas->SetGrid(1, 1);
   Graphic_setup(pr.ptge_dt_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_dt_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_dt_mean_XP->SetTitle(";Drift time (timebins);Mean (ADC counts/cm)");
   pr.ptge_dt_mean_XP->GetXaxis()->SetLimits(0, 510);
   pr.ptge_dt_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_dt_mean_XP->Draw("AP");
   pr.ptge_dt_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   gPad->SetTopMargin(0.02);
   gPad->SetRightMargin(0.1);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   // heatmap WF vs time bin
   pr.fph2f_WFdrift->Draw("colz");
   TGraphErrors *ptge_mom_WFdrift = Convert_TH2_TGE(pr.fph2f_WFdrift);
   TF1 fitWFtime("fitWFtime", "pol1", 80, 350);
   ptge_mom_WFdrift->Fit(&fitWFtime, "RQ");
   fitWFtime.SetLineColor(kRed);
   fitWFtime.Draw("same");
   pavetext.Clear();
   pavetext.AddText(
      Form("y = %.2fx + %.2f", fitWFtime.GetParameter(1), fitWFtime.GetParameter(0)));
   pavetext.DrawClone();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // heatmap XP vs time bin
   pr.fph2f_XPdrift->Draw("colz");
   TGraphErrors *ptge_mom_XPdrift = Convert_TH2_TGE(pr.fph2f_XPdrift);
   TF1 fitXPtime("fitXPtime", "pol1", 80, 350);
   ptge_mom_XPdrift->Fit(&fitXPtime, "RQ");
   fitXPtime.SetLineColor(kRed);
   fitXPtime.Draw("same");
   pavetext.Clear();
   pavetext.AddText(
      Form("y = %.2fx + %.2f", fitXPtime.GetParameter(1), fitXPtime.GetParameter(0)));
   pavetext.DrawClone();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // heatmap XP vs time bin module VS module
   pavetext.SetY1NDC(0.7);
   pavetext.SetY2NDC(0.9);
   std::vector<TF1> vlinearFit;
   for (int i = 0; i < 32; i++) {
      vlinearFit.push_back(TF1(Form("linearFit_%d", i), "pol1", 80, 350));
      vlinearFit[i].SetLineColor(kRed);
      vlinearFit[i].SetLineWidth(1);
      if (i % 16 == 0) {
         fpCanvas->Clear();
         fpCanvas->Divide(4, 4);
      }
      fpCanvas->cd(i % 16 + 1);
      fpCanvas->cd(i % 16 + 1)->SetRightMargin(0.1);
      pr.vmod_fph2f_XPdrift[i]->Draw("colz");
      pavetext.Clear();
      pr.vmod_fph2f_XPdrift[i]->Fit(&vlinearFit[i], "RQ");
      vlinearFit[i].Draw("same");
      pavetext.AddText(Form("y = %.2fx + %.2f", vlinearFit[i].GetParameter(1),
                            vlinearFit[i].GetParameter(0)));
      pavetext.DrawClone();
      if (i % 16 == 15)
         fpCanvas->SaveAs(fRealpathPDF.c_str());
   }

   // Time bin vs X track position
   fpCanvas->Clear();
   pr.fph2f_timeX->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs track length
   // Resolution
   fpCanvas->Clear();
   fpCanvas->SetGrid(1, 1);
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_trklen_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_trklen_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_trklen_reso_XP->SetTitle(";Track length (mm);dE/dx resolution (%)");
   pr.ptge_trklen_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_trklen_reso_XP->DrawClone("AP");
   pr.ptge_trklen_reso_WF->DrawClone("P same");
   pr.ptge_trklen_reso_XP->SetMarkerSize(7);
   pr.ptge_trklen_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   fpCanvas->SetGrid(1, 1);
   Graphic_setup(pr.ptge_trklen_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_trklen_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_trklen_mean_XP->SetTitle(";Track length (mm);Mean (ADC counts/cm)");
   pr.ptge_trklen_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_trklen_mean_XP->Draw("AP");
   pr.ptge_trklen_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   gPad->SetTopMargin(0.02);
   gStyle->SetStatX(0.87);
   pr.fph2f_XPlen->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs absolute track angle phi
   // Resolution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_absphi_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_absphi_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_absphi_reso_XP->SetTitle(";#varphi (#circ);dE/dx resolution (%)");
   pr.ptge_absphi_reso_XP->GetXaxis()->SetLimits(0, 90);
   pr.ptge_absphi_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_absphi_reso_XP->DrawClone("AP");
   pr.ptge_absphi_reso_WF->DrawClone("P same");
   pr.ptge_absphi_reso_XP->SetMarkerSize(7);
   pr.ptge_absphi_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_absphi_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_absphi_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_absphi_mean_XP->SetTitle(";#varphi (#circ);Mean (ADC counts/cm)");
   pr.ptge_absphi_mean_XP->GetXaxis()->SetLimits(0, 90);
   pr.ptge_absphi_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_absphi_mean_XP->Draw("AP");
   pr.ptge_absphi_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs track angle phi
   // Resolution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_phi_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_phi_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_phi_reso_XP->SetTitle(";#varphi (#circ);dE/dx resolution (%)");
   pr.ptge_phi_reso_XP->GetXaxis()->SetLimits(-90, 90);
   pr.ptge_phi_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_phi_reso_XP->DrawClone("AP");
   pr.ptge_phi_reso_WF->DrawClone("P same");
   pr.ptge_phi_reso_XP->SetMarkerSize(7);
   pr.ptge_phi_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_phi_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_phi_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_phi_mean_XP->SetTitle(";#varphi (#circ);Mean (ADC counts/cm)");
   pr.ptge_phi_mean_XP->GetXaxis()->SetLimits(-90, 90);
   pr.ptge_phi_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_phi_mean_XP->Draw("AP");
   pr.ptge_phi_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution
   fpCanvas->Clear();
   gPad->SetTopMargin(0.02);
   gPad->SetRightMargin(0.13);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   pr.fph2f_XPphi->SetTitle(";#varphi (#circ);dE/dx with XP (ADC counts/cm)");
   pr.fph2f_XPphi->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_lenphi->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_chi2ndfphi->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_momphi->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   gPad->SetRightMargin(0.02);
   pr.fph1f_phi->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // dE/dx vs absolute track angle theta
   // -----------------------------------------------------------------------------------------------------------------
   // Resolution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_abstheta_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_abstheta_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_abstheta_reso_XP->SetTitle(";#theta (#circ);dE/dx resolution (%)");
   pr.ptge_abstheta_reso_XP->GetXaxis()->SetLimits(0, 90);
   pr.ptge_abstheta_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_abstheta_reso_XP->DrawClone("AP");
   pr.ptge_abstheta_reso_WF->DrawClone("P same");
   pr.ptge_abstheta_reso_XP->SetMarkerSize(7);
   pr.ptge_abstheta_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_abstheta_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_abstheta_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_abstheta_mean_XP->SetTitle(";#theta (#circ);Mean (ADC counts/cm)");
   pr.ptge_abstheta_mean_XP->GetXaxis()->SetLimits(0, 90);
   pr.ptge_abstheta_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_abstheta_mean_XP->Draw("AP");
   pr.ptge_abstheta_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // dE/dx vs track angle theta
   // -----------------------------------------------------------------------------------------------------------------
   // Resolution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_theta_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_theta_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_theta_reso_XP->SetTitle(";#theta (#circ);dE/dx resolution (%)");
   pr.ptge_theta_reso_XP->GetXaxis()->SetLimits(-90, 90);
   pr.ptge_theta_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_theta_reso_XP->DrawClone("AP");
   pr.ptge_theta_reso_WF->DrawClone("P same");
   pr.ptge_theta_reso_XP->SetMarkerSize(7);
   pr.ptge_theta_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_theta_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_theta_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_theta_mean_XP->SetTitle(";#theta (#circ);Mean (ADC counts/cm)");
   pr.ptge_theta_mean_XP->GetXaxis()->SetLimits(-90, 90);
   pr.ptge_theta_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_theta_mean_XP->Draw("AP");
   pr.ptge_theta_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   pr.fph2f_XPtheta->SetTitle(";#theta (#circ);dE/dx with XP (ADC counts/cm)");
   pr.fph2f_XPtheta->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_lentheta->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_momtheta->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   pr.fph1f_theta->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // Both angles
   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   pr.fph2f_phitheta->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dEdx vs absolute momentum
   // Resolution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   gPad->SetTopMargin(0.03);
   Graphic_setup(pr.ptge_absmom_reso_WF, 3, 33, kCyan + 2, 1, kCyan + 2);
   Graphic_setup(pr.ptge_absmom_reso_XP, 3, 47, kMagenta + 2, 1, kMagenta + 2);
   pr.ptge_absmom_reso_XP->SetTitle(";Momentum (MeV/c);dE/dx resolution (%)");
   pr.ptge_absmom_reso_XP->GetXaxis()->SetLimits(0, absmomrange);
   pr.ptge_absmom_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_absmom_reso_XP->DrawClone("AP");
   pr.ptge_absmom_reso_WF->DrawClone("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs((fRealpathPDF + "(").c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_absmom_mean_WF, 0, 33, kCyan + 2, 1, kCyan + 2);
   Graphic_setup(pr.ptge_absmom_mean_XP, 0, 47, kMagenta + 2, 1, kMagenta + 2);
   pr.ptge_absmom_mean_XP->SetTitle(";Momentum (MeV/c);Mean (ADC counts/cm)");
   pr.ptge_absmom_mean_XP->GetXaxis()->SetLimits(0, absmomrange);
   pr.ptge_absmom_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_absmom_mean_XP->Draw("APL");
   pr.ptge_absmom_mean_WF->Draw("PL same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   pr.fph2f_WFabsmom->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_XPabsmom->Draw("colz");

   // -----------------------------------------------------------------------------------------------------------------
   // dEdx vs momentum
   // Resolution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   gPad->SetTopMargin(0.03);
   Graphic_setup(pr.ptge_mom_reso_WF, 3, 33, kCyan + 2, 1, kCyan + 2);
   Graphic_setup(pr.ptge_mom_reso_XP, 3, 47, kMagenta + 2, 1, kMagenta + 2);
   pr.ptge_mom_reso_XP->SetTitle(";Momentum (MeV/c);dE/dx resolution (%)");
   pr.ptge_mom_reso_XP->GetXaxis()->SetLimits(-momrange, momrange);
   pr.ptge_mom_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_mom_reso_XP->DrawClone("AP");
   pr.ptge_mom_reso_WF->DrawClone("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs((fRealpathPDF + "(").c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_mom_mean_WF, 0, 33, kCyan + 2, 1, kCyan + 2);
   Graphic_setup(pr.ptge_mom_mean_XP, 0, 47, kMagenta + 2, 1, kMagenta + 2);
   pr.ptge_mom_mean_XP->SetTitle(";Momentum (MeV/c);Mean (ADC counts/cm)");
   pr.ptge_mom_mean_XP->GetXaxis()->SetLimits(-momrange, momrange);
   pr.ptge_mom_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_mom_mean_XP->Draw("APL");
   pr.ptge_mom_mean_WF->Draw("PL same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   pr.fph2f_WFmom->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_XPmom->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   gStyle->SetOptStat(0);
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.fph1i_mom, 0.5, 1, kMagenta + 2, 2, kMagenta - 2, kMagenta, 0.2);
   Graphic_setup(pr.fph1i_mom_tHAT, 0.5, 1, kRed + 2, 2, kRed - 2, kRed, 0.2);
   Graphic_setup(pr.fph1i_mom_bHAT, 0.5, 1, kBlue + 2, 2, kBlue - 2, kBlue, 0.2);
   // pr.fph1i_mom->				SetAxisRange(0, 1e4, "Y");
   pr.fph1i_mom->Draw();
   pr.fph1i_mom_bHAT->Draw("same");
   pr.fph1i_mom_tHAT->Draw("same");
   TLegend *plegmom = new TLegend(0.7, 0.7, 0.9, 0.9);
   plegmom->AddEntry(pr.fph1i_mom, "Both HATs", "l");
   plegmom->AddEntry(pr.fph1i_mom_tHAT, "tHAT", "l");
   plegmom->AddEntry(pr.fph1i_mom_bHAT, "bHAT", "l");
   plegmom->SetTextSize(0.06);
   plegmom->SetTextColor(kBlue - 1);
   plegmom->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());
   delete plegmom;

   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   pr.fph2f_chi2mom->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // Pulls
   // -----------------------------------------------------------------------------------------------------------------
   // Pulls vs theta
   // Standard deviation
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_theta_std_pullmu, 2, 21, kBlue + 1, 2, kBlue + 1);
   Graphic_setup(pr.ptge_theta_std_pullelec, 2, 20, kOrange + 7, 2, kOrange + 7);
   pr.ptge_theta_std_pullmu->SetTitle(";#theta (#circ);Pulls standard deviation");
   float minpull = 0;
   float maxpull = (pr.ptge_theta_std_pullmu->GetHistogram()->GetMaximum() >
                    pr.ptge_theta_std_pullelec->GetHistogram()->GetMaximum())
                      ? pr.ptge_theta_std_pullmu->GetHistogram()->GetMaximum()
                      : pr.ptge_theta_std_pullelec->GetHistogram()->GetMaximum();
   float diffpull = maxpull - minpull;
   pr.ptge_theta_std_pullmu->GetXaxis()->SetLimits(-90, 90);
   pr.ptge_theta_std_pullmu->GetYaxis()->SetRangeUser(0, maxpull + 0.1 * diffpull);
   pr.ptge_theta_std_pullmu->DrawClone("AP");
   pr.ptge_theta_std_pullelec->DrawClone("P same");
   pr.ptge_theta_std_pullmu->SetMarkerSize(7);
   pr.ptge_theta_std_pullelec->SetMarkerSize(7);
   TLegend legpull(0.45, 0.75, 0.65, 0.95);
   legpull.SetTextAlign(22);
   legpull.SetNColumns(2);
   legpull.AddEntry(pr.ptge_theta_std_pullmu, " #mu^{#pm}", "p");
   legpull.AddEntry(pr.ptge_theta_std_pullelec, " e^{#pm}", "p");
   legpull.SetTextSize(0.06);
   legpull.SetTextColor(kBlue - 1);
   legpull.SetFillStyle(0);
   legpull.Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_theta_mean_pullmu, 2, 21, kBlue + 1, 2, kBlue + 1);
   Graphic_setup(pr.ptge_theta_mean_pullelec, 2, 20, kOrange + 7, 2, kOrange + 7);
   pr.ptge_theta_mean_pullmu->SetTitle(";#theta (#circ);Pulls mean");
   minpull = (pr.ptge_theta_mean_pullmu->GetHistogram()->GetMinimum() <
              pr.ptge_theta_mean_pullelec->GetHistogram()->GetMinimum())
                ? pr.ptge_theta_mean_pullmu->GetHistogram()->GetMinimum()
                : pr.ptge_theta_mean_pullelec->GetHistogram()->GetMinimum();
   maxpull = (pr.ptge_theta_mean_pullmu->GetHistogram()->GetMaximum() >
              pr.ptge_theta_mean_pullelec->GetHistogram()->GetMaximum())
                ? pr.ptge_theta_mean_pullmu->GetHistogram()->GetMaximum()
                : pr.ptge_theta_mean_pullelec->GetHistogram()->GetMaximum();
   pr.ptge_theta_mean_pullmu->GetXaxis()->SetLimits(-90, 90);
   pr.ptge_theta_mean_pullmu->GetYaxis()->SetRangeUser(minpull - 0.1 * diffpull,
                                                       maxpull + 0.1 * diffpull);
   pr.ptge_theta_mean_pullmu->Draw("AP");
   pr.ptge_theta_mean_pullelec->Draw("P same");
   legpull.SetX1NDC(0.15);
   legpull.SetX2NDC(0.35);
   legpull.SetY1NDC(0.1);
   legpull.SetY2NDC(0.3);
   legpull.Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // Pulls vs momentum
   // Standard deviation
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   Graphic_setup(pr.ptge_mom_std_pullmu, 2, 21, kBlue + 1, 2, kBlue + 1);
   Graphic_setup(pr.ptge_mom_std_pullelec, 2, 20, kOrange + 7, 2, kOrange + 7);
   pr.ptge_mom_std_pullmu->SetTitle(";Momentum (MeV/c);Pulls standard deviation");
   minpull = 0;
   maxpull = (pr.ptge_mom_std_pullmu->GetHistogram()->GetMaximum() >
              pr.ptge_mom_std_pullelec->GetHistogram()->GetMaximum())
                ? pr.ptge_mom_std_pullmu->GetHistogram()->GetMaximum()
                : pr.ptge_mom_std_pullelec->GetHistogram()->GetMaximum();
   pr.ptge_mom_std_pullmu->GetXaxis()->SetLimits(-momrange, momrange);
   pr.ptge_mom_std_pullmu->GetYaxis()->SetRangeUser(0, maxpull + 0.1 * diffpull);
   pr.ptge_mom_std_pullmu->DrawClone("AP");
   pr.ptge_mom_std_pullelec->DrawClone("P same");
   pr.ptge_mom_std_pullmu->SetMarkerSize(7);
   pr.ptge_mom_std_pullelec->SetMarkerSize(7);
   legpull.SetX1NDC(0.75);
   legpull.SetX2NDC(0.95);
   legpull.SetY1NDC(0.75);
   legpull.SetY2NDC(0.9);
   legpull.Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_mom_mean_pullmu, 2, 21, kBlue + 1, 2, kBlue + 1);
   Graphic_setup(pr.ptge_mom_mean_pullelec, 2, 20, kOrange + 7, 2, kOrange + 7);
   pr.ptge_mom_mean_pullmu->SetTitle(";Momentum (MeV/c);Pulls mean");
   minpull = (pr.ptge_mom_mean_pullmu->GetHistogram()->GetMinimum() <
              pr.ptge_mom_mean_pullelec->GetHistogram()->GetMinimum())
                ? pr.ptge_mom_mean_pullmu->GetHistogram()->GetMinimum()
                : pr.ptge_mom_mean_pullelec->GetHistogram()->GetMinimum();
   maxpull = (pr.ptge_mom_mean_pullmu->GetHistogram()->GetMaximum() >
              pr.ptge_mom_mean_pullelec->GetHistogram()->GetMaximum())
                ? pr.ptge_mom_mean_pullmu->GetHistogram()->GetMaximum()
                : pr.ptge_mom_mean_pullelec->GetHistogram()->GetMaximum();
   pr.ptge_mom_mean_pullmu->GetXaxis()->SetLimits(-momrange, momrange);
   pr.ptge_mom_mean_pullmu->GetYaxis()->SetRangeUser(minpull - 0.1 * diffpull,
                                                     maxpull + 0.1 * diffpull);
   pr.ptge_mom_mean_pullmu->Draw("AP");
   pr.ptge_mom_mean_pullelec->Draw("P same");
   legpull.SetX1NDC(0.75);
   legpull.SetX2NDC(0.95);
   legpull.SetY1NDC(0.25);
   legpull.SetY2NDC(0.45);
   legpull.Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   pr.fph2f_pullelecmu->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Pulls distribution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);
   gPad->SetTopMargin(0.05);
   Graphic_setup(pr.fph1f_pullmu, 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan, 0.2);
   Graphic_setup(pr.fph1f_pullelec, 0.5, 1, kOrange + 7, 2, kOrange - 2, kOrange, 0.2);
   Graphic_setup(pr.fph1f_pullproton, 0.5, 1, kRed + 1, 2, kRed - 2, kRed, 0.2);
   TLegend legpull2(0.75, 0.7, 0.9, 0.93);
   legpull2.AddEntry(pr.fph1f_pullmu, " #mu hyp.", "f");
   legpull2.AddEntry(pr.fph1f_pullelec, " e hyp.", "f");
   legpull2.AddEntry(pr.fph1f_pullproton, " p hyp.", "f");
   legpull2.SetTextSize(0.06);
   legpull2.SetTextColor(kBlue - 1);
   pr.fph1f_pullmu->SetXTitle("Pull");
   pr.fph1f_pullmu->Draw("HIST");
   pr.fph1f_pullelec->Draw("HIST same");
   pr.fph1f_pullproton->Draw("HIST same");
   TF1 tf1_mu = *Fit1Gauss(pr.fph1f_pullmu);
   TF1 tf1_elec = *Fit1Gauss(pr.fph1f_pullelec);
   float meanmu = tf1_mu.GetParameter(1);
   float dmeanmu = tf1_mu.GetParError(1);
   float meanelec = tf1_elec.GetParameter(1);
   float dmeanelec = tf1_elec.GetParError(1);
   float stdmu = tf1_mu.GetParameter(2);
   float dstdmu = tf1_mu.GetParError(2);
   float stdelec = tf1_elec.GetParameter(2);
   float dstdelec = tf1_elec.GetParError(2);
   float separation = GetSeparation(&tf1_mu, &tf1_elec);
   float dseparation = GetSeparationError(&tf1_mu, &tf1_elec);
   TPaveText pavetextmu(0.6, 0.4, 0.93, 0.55, "NDC");
   pavetextmu.SetTextSize(0.05);
   pavetextmu.SetTextColor(kCyan + 4);
   pavetextmu.SetFillColorAlpha(kCyan, 0.2);
   pavetextmu.SetLineWidth(2);
   pavetextmu.SetBorderSize(2);
   pavetextmu.SetLineColor(kCyan - 2);
   pavetextmu.AddText(Form("#mu_{#mu} = %.3f #pm %.3f", meanmu, dmeanmu));
   pavetextmu.AddText(Form("#sigma_{#mu} = %.3f #pm %.3f", stdmu, dstdmu));

   TPaveText pavetextelec(0.6, 0.2, 0.93, 0.35, "NDC");
   pavetextelec.SetTextSize(0.05);
   pavetextelec.SetTextColor(kOrange + 3);
   pavetextelec.SetFillColorAlpha(kOrange, 0.2);
   pavetextelec.SetLineWidth(2);
   pavetextelec.SetBorderSize(2);
   pavetextelec.SetLineColor(kOrange - 2);
   pavetextelec.AddText(Form("#mu_{e} = %.3f #pm %.3f", meanelec, dmeanelec));
   pavetextelec.AddText(Form("#sigma_{e} = %.3f #pm %.3f", stdelec, dstdelec));

   TLatex latexpull;
   latexpull.SetNDC();
   latexpull.SetTextSize(0.06);
   latexpull.SetTextColor(kBlue - 1);
   latexpull.DrawLatex(0.6, 0.6,
                       Form("S(#mu/e) = %.2f #pm %.2f", separation, dseparation));
   legpull2.Draw();
   pavetextmu.Draw();
   pavetextelec.Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Geometry plots
   // -----------------------------------------------------------------------------------------------------------------

   fpCanvas->Clear();
   gPad->SetTopMargin(0.02);
   gPad->SetLogz();
   pr.fph2f_XZ->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());
   gPad->SetLogz(0);
   
   fpCanvas->Clear();
   gPad->SetTopMargin(0.02);
   gPad->SetRightMargin(0.02);
   pr.fph1f_trklen->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph1f_chi2->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Debugging
   // -----------------------------------------------------------------------------------------------------------------

   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   pr.fph2f_momR->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_chi2ndfR->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // dE/dx VS momentum bin
   // -----------------------------------------------------------------------------------------------------------------
   gPad->SetTopMargin(0.08);
   gPad->SetRightMargin(0.02);
   for (int i = 0; i < nmombins; i++) {
      fpCanvas->Clear();
      pr.vmom_fph1f_WF[i]->SetAxisRange(0,
                                        1.1 *
                                           std::max({pr.vmom_fph1f_WF[i]->GetMaximum(),
                                                     pr.vmom_fph1f_XP[i]->GetMaximum()}),
                                        "Y");
      float invX = 0;
      if (pr.vmom_fph1f_WF[i]->GetMean() > dEdxmax / 2)
         invX = 0.4;
      Graphic_setup(pr.vmom_fph1f_WF[i], 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan, 0.2);
      Graphic_setup(pr.vmom_fph1f_XP[i], 0.5, 1, kMagenta + 1, 2, kMagenta - 2, kMagenta,
                    0.2);
      pr.vmom_fph1f_WF[i]->Draw("HIST");
      pr.vmom_fph1f_XP[i]->Draw("HIST same");
      PrintResolution(pr.vmom_fph1f_XP[i], fpCanvas, 0.65 - invX, 0.58, kMagenta + 2,
                      "XP");
      PrintResolution(pr.vmom_fph1f_WF[i], fpCanvas, 0.65 - invX, 0.25, kCyan + 2, "WF");
      fpCanvas->SaveAs(fRealpathPDF.c_str());
   }

   // dE/dx VS theta bin
   // -----------------------------------------------------------------------------------------------------------------
   gPad->SetTopMargin(0.08);
   gPad->SetRightMargin(0.02);
   for (int i = 0; i < nthetabins; i++) {
      fpCanvas->Clear();
      pr.vtheta_fph1f_WF[i]->SetAxisRange(
         0,
         1.1 * std::max({pr.vtheta_fph1f_WF[i]->GetMaximum(),
                         pr.vtheta_fph1f_XP[i]->GetMaximum()}),
         "Y");
      float invX = 0;
      if (pr.vtheta_fph1f_WF[i]->GetMean() > dEdxmax / 2)
         invX = 0.4;
      Graphic_setup(pr.vtheta_fph1f_WF[i], 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan, 0.2);
      Graphic_setup(pr.vtheta_fph1f_XP[i], 0.5, 1, kMagenta + 1, 2, kMagenta - 2,
                    kMagenta, 0.2);
      pr.vtheta_fph1f_WF[i]->Draw("");
      pr.vtheta_fph1f_XP[i]->Draw("same");
      PrintResolution(pr.vtheta_fph1f_XP[i], fpCanvas, 0.65 - invX, 0.58, kMagenta + 2,
                      "XP");
      PrintResolution(pr.vtheta_fph1f_WF[i], fpCanvas, 0.65 - invX, 0.25, kCyan + 2,
                      "WF");
      fpCanvas->SaveAs(fRealpathPDF.c_str());
   }

   // Drift velocity plots
   // -----------------------------------------------------------------------------------------------------------------
   fpCanvas->Clear();
   gPad->SetRightMargin(0.02);
   int bins = 100;
   int tmin = -999;
   int max = -999;
   for (int i = 0; i < bins; i++) {
      if (pr.fph1i_tminEP0->GetBinContent(i) > max) {
         tmin = i;
         max = pr.fph1i_tminEP0->GetBinContent(i);
      }
   }

   int tmax = -999;
   for (int i = 300; i < 300 + bins; i++) {
      if (pr.fph1i_tmaxBotCath->GetBinContent(i) > max) {
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
   // std::cout << "Drift time in bHATPC: " << Tdrift << " +/- " << dTdrift << " ns" <<
   // std::endl; float driftVel =			98.5/Tdrift;
   // float dtriftVel =			driftVel * dTdrift/Tdrift;
   // std::cout << "Drift velocity in bHATPC: " << driftVel << " +/- " << dtriftVel << "
   // cm/µs" << std::endl;

   // Tmin =						ptf1_tminEP2->GetParameter(1)*40/1000;
   // dTmin =						ptf1_tminEP2->GetParError(1)*40/1000;
   // Tmax =						ptf1_tmaxTopCath->GetParameter(1)*40/1000;
   // dTmax =						ptf1_tmaxTopCath->GetParError(1)*40/1000;
   // Tdrift =					(Tmax - Tmin);
   // dTdrift =					sqrt(dTmin*dTmin + dTmax*dTmax);
   // std::cout << "Drift time in tHATPC: " << Tdrift << " +/- " << dTdrift << " ns" <<
   // std::endl; driftVel =					98.5/Tdrift;
   // dtriftVel =					driftVel * dTdrift/Tdrift;
   // std::cout << "Drift velocity in tHATPC: " << driftVel << " +/- " << dtriftVel << "
   // cm/µs" << std::endl;

   pr.fph1i_tminBotCath->SetAxisRange(
      0,
      1.1 * std::max({pr.fph1i_tminBotCath->GetMaximum(), pr.fph1i_tminEP0->GetMaximum(),
                      pr.fph1i_tminEP1->GetMaximum(), pr.fph1i_tminTopCath->GetMaximum(),
                      pr.fph1i_tminEP2->GetMaximum(), pr.fph1i_tminEP3->GetMaximum()}),
      "Y");
   pr.fph1i_tminBotCath->SetLineWidth(2);
   pr.fph1i_tminBotCath->SetLineColor(kGreen + 2);
   pr.fph1i_tminEP0->SetLineWidth(2);
   pr.fph1i_tminEP0->SetLineColor(kBlue);
   pr.fph1i_tminEP1->SetLineWidth(2);
   pr.fph1i_tminEP1->SetLineColor(kRed);
   pr.fph1i_tminTopCath->SetLineWidth(2);
   pr.fph1i_tminTopCath->SetLineColor(kGreen - 6);
   pr.fph1i_tminEP2->SetLineWidth(2);
   pr.fph1i_tminEP2->SetLineColor(kBlue - 6);
   pr.fph1i_tminEP3->SetLineWidth(2);
   pr.fph1i_tminEP3->SetLineColor(kRed - 6);
   pr.fph1i_tminBotCath->Draw();
   pr.fph1i_tminEP0->Draw("same");
   pr.fph1i_tminEP1->Draw("same");
   pr.fph1i_tminTopCath->Draw("same");
   pr.fph1i_tminEP2->Draw("same");
   pr.fph1i_tminEP3->Draw("same");
   // ptf1_tminEP0->				SetLineColor(kOrange+2);
   // ptf1_tminEP2->				SetLineColor(kViolet+2);
   // ptf1_tminEP0->				Draw("same");
   // ptf1_tminEP2->				Draw("same");
   TLegend *legtmin = new TLegend(0.6, 0.5, 0.9, 0.9);
   legtmin->AddEntry(pr.fph1i_tminBotCath, "T_{min} (bottom cathode) ", "l");
   legtmin->AddEntry(pr.fph1i_tminEP0, "T_{min} EP0 ", "l");
   legtmin->AddEntry(pr.fph1i_tminEP1, "T_{min} EP1 ", "l");
   legtmin->AddEntry(pr.fph1i_tminTopCath, "T_{min} (top cathode) ", "l");
   legtmin->AddEntry(pr.fph1i_tminEP2, "T_{min} EP2 ", "l");
   legtmin->AddEntry(pr.fph1i_tminEP3, "T_{min} EP3 ", "l");
   legtmin->DrawClone();
   fpCanvas->Update();
   delete legtmin;
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   pr.fph1i_tmaxBotCath->SetAxisRange(
      0,
      1.1 * std::max({pr.fph1i_tmaxBotCath->GetMaximum(), pr.fph1i_tmaxEP0->GetMaximum(),
                      pr.fph1i_tmaxEP1->GetMaximum(), pr.fph1i_tmaxTopCath->GetMaximum(),
                      pr.fph1i_tmaxEP2->GetMaximum(), pr.fph1i_tmaxEP3->GetMaximum()}),
      "Y");
   pr.fph1i_tmaxBotCath->SetLineWidth(2);
   pr.fph1i_tmaxBotCath->SetLineColor(kGreen + 2);
   pr.fph1i_tmaxEP0->SetLineWidth(2);
   pr.fph1i_tmaxEP0->SetLineColor(kBlue);
   pr.fph1i_tmaxEP1->SetLineWidth(2);
   pr.fph1i_tmaxEP1->SetLineColor(kRed);
   pr.fph1i_tmaxTopCath->SetLineWidth(2);
   pr.fph1i_tmaxTopCath->SetLineColor(kGreen - 6);
   pr.fph1i_tmaxEP2->SetLineWidth(2);
   pr.fph1i_tmaxEP2->SetLineColor(kBlue - 6);
   pr.fph1i_tmaxEP3->SetLineWidth(2);
   pr.fph1i_tmaxEP3->SetLineColor(kRed - 6);
   pr.fph1i_tmaxBotCath->Draw();
   pr.fph1i_tmaxEP0->Draw("same");
   pr.fph1i_tmaxEP1->Draw("same");
   pr.fph1i_tmaxTopCath->Draw("same");
   pr.fph1i_tmaxEP2->Draw("same");
   pr.fph1i_tmaxEP3->Draw("same");
   // ptf1_tmaxBotCath->			SetLineColor(kOrange+2);
   // ptf1_tmaxTopCath->			SetLineColor(kViolet+2);
   // ptf1_tmaxTopCath->			Draw("same");
   // ptf1_tmaxBotCath->			Draw("same");
   TLegend *legtmax = new TLegend(0.15, 0.5, 0.4, 0.9);
   legtmax->AddEntry(pr.fph1i_tmaxBotCath, "T_{max} (bottom cathode) ", "l");
   legtmax->AddEntry(pr.fph1i_tmaxEP0, "T_{max} EP0 ", "l");
   legtmax->AddEntry(pr.fph1i_tmaxEP1, "T_{max} EP1 ", "l");
   legtmax->AddEntry(pr.fph1i_tmaxTopCath, "T_{max} (top cathode) ", "l");
   legtmax->AddEntry(pr.fph1i_tmaxEP2, "T_{max} EP2 ", "l");
   legtmax->AddEntry(pr.fph1i_tmaxEP3, "T_{max} EP3 ", "l");
   legtmax->DrawClone();
   fpCanvas->Update();
   delete legtmax;
   fpCanvas->SaveAs((fRealpathPDF + ")").c_str());
   fpCanvas->Clear();

   // Reset style
   gPad->UseCurrentStyle();
}

void Draw::Compare(const std::vector<std::string> &v_filepaths, const std::string &type)
{

   // Load processes
   std::vector<Process *> v_processes;
   for (int i = 0; i < v_filepaths.size(); i++) {
      std::cout << "File " << i << ": " << v_filepaths[i] << std::endl;
      TFile inputFile(v_filepaths[i].c_str(), "READ");
      Process *process = dynamic_cast<Process *>(inputFile.Get("Process"));
      v_processes.push_back(process);
   }

   Process *pr0 = v_processes[0];
   int ncomparisons = (int)v_processes.size();
   std::string OutputFile;
   if (type == "comments") {
      OutputFile =
         foutputComparisonFolder + "/Comparison_" + pr0->frun + pr0->ftag + pr0->fcomment;
      for (int i = 1; i < ncomparisons; i++)
         OutputFile += "_VS_" + v_processes[i]->fcomment;
      OutputFile += pr0->fcutslist + ".pdf";
   } else if (type == "cuts") {
      OutputFile = foutputComparisonFolder + "/" + pr0->fcomment + "/Comparison_" +
                   pr0->frun + pr0->ftag + pr0->fcomment + pr0->fcutslist;
      for (int i = 1; i < ncomparisons; i++)
         OutputFile += "_VS_" + v_processes[i]->fcutslist;
      OutputFile += ".pdf";
   } else if (type == "runs") {
      if (pr0->ffileName.find("dog1") != std::string::npos or
          pr0->ffileName.find("cosmics") != std::string::npos) {
         OutputFile =
            foutputComparisonFolder + "/Cosmics_Data_MC" + pr0->fcutslist + ".pdf";
      } else if (pr0->ffileName.find("beam") != std::string::npos or
                 pr0->ffileName.find("sandmu") != std::string::npos) {
         OutputFile = foutputComparisonFolder + "/Beam_Data_MC" + pr0->fcutslist + ".pdf";
      }
   }

   // Plot setup
   // -----------------------------------------------------------------------------------------------------------------
   int dEdxmax = pr0->dEdxmax;
   int momrange = pr0->momrange;
   int nmombins = pr0->nmombins;
   fpCanvas->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gPad->SetTopMargin(0.05);
   gPad->SetRightMargin(0.035);

   int resominmom = 2;
   int resomaxmom = 12;
   int dedxminmom = 350;
   int dedxmaxmom = 650;
   int resominX = 1;
   int resomaxX = 12;
   int dedxminX = 350;
   int dedxmaxX = 500;
   int resomintrklen = 4;
   int resomaxtrklen = 16;
   int dedxmintrklen = 300;
   int dedxmaxtrklen = 550;
   int resominphi = 2;
   int resomaxphi = 12;
   int dedxminphi = 350;
   int dedxmaxphi = 500;
   int resomintheta = 3;
   int resomaxtheta = 13;
   int dedxmintheta = 350;
   int dedxmaxtheta = 500;
   if (pr0->ffileName.find("dog1") != std::string::npos or
       pr0->ffileName.find("cosmics") != std::string::npos) {
      resominmom = 6;
      resomaxmom = 14;
      dedxminmom = 350;
      dedxmaxmom = 650;
      resominX = 4;
      resomaxX = 16;
      dedxminX = 300;
      dedxmaxX = 550;
      resomintrklen = 4;
      resomaxtrklen = 16;
      dedxmintrklen = 300;
      dedxmaxtrklen = 550;
      resominphi = 4;
      resomaxphi = 14;
      dedxminphi = 350;
      dedxmaxphi = 500;
      resomintheta = 4;
      resomaxtheta = 16;
      dedxmintheta = 300;
      dedxmaxtheta = 550;
   }

   // Legend
   fpLegend = new TLegend(0.65, 0.85 - 0.06 * ncomparisons, 0.9, 0.92);
   fpLegend->SetTextSize(0.06 - 0.005 * ncomparisons);
   // fpLegend->						SetFillStyle(0);
   fpLegend->SetTextColor(kBlue - 1);
   fpLegend->SetFillColorAlpha(kWhite, 0.8);

   float ampmax = 0;
   for (int i = 0; i < ncomparisons; i++) {
      Process *pr_tmp = v_processes[i];
      Graphic_setup(pr_tmp->fph1f_WF, 0.5, 1, colors[i] + 1, 2, colors[i] + 1, colors[i],
                    0.2);
      Graphic_setup(pr_tmp->fph1f_XP, 0.5, 1, colors[i] + 1, 2, colors[i] + 1, colors[i],
                    0.2);
      if (pr_tmp->fph1f_WF->GetMaximum() > ampmax)
         ampmax = pr_tmp->fph1f_WF->GetMaximum();
      if (pr_tmp->fph1f_XP->GetMaximum() > ampmax)
         ampmax = pr_tmp->fph1f_XP->GetMaximum();
      Graphic_setup(pr_tmp->ptge_absmom_reso_WF, 3, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_absmom_reso_XP, 3, markers[i], colors[i], 2, colors[i]);
      Graphic_setup(pr_tmp->ptge_absmom_mean_WF, 2, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_absmom_mean_XP, 2, markers[i], colors[i], 2, colors[i]);
      Graphic_setup(pr_tmp->ptge_dd_reso_WF, 3, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_dd_reso_XP, 3, markers[i], colors[i], 2, colors[i]);
      Graphic_setup(pr_tmp->ptge_dd_mean_WF, 2, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_dd_mean_XP, 2, markers[i], colors[i], 2, colors[i]);
      Graphic_setup(pr_tmp->ptge_trklen_reso_WF, 3, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_trklen_reso_XP, 3, markers[i], colors[i], 2,
                    colors[i]);
      Graphic_setup(pr_tmp->ptge_trklen_mean_WF, 2, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_trklen_mean_XP, 2, markers[i], colors[i], 2,
                    colors[i]);
      Graphic_setup(pr_tmp->ptge_absphi_reso_WF, 3, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_absphi_reso_XP, 3, markers[i], colors[i], 2, colors[i]);
      Graphic_setup(pr_tmp->ptge_absphi_mean_WF, 2, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_absphi_mean_XP, 2, markers[i], colors[i], 2, colors[i]);
      Graphic_setup(pr_tmp->ptge_abstheta_reso_WF, 3, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_abstheta_reso_XP, 3, markers[i], colors[i], 2,
                    colors[i]);
      Graphic_setup(pr_tmp->ptge_abstheta_mean_WF, 2, markers[i], colors[2 + i], 2,
                    colors[2 + i]);
      Graphic_setup(pr_tmp->ptge_abstheta_mean_XP, 2, markers[i], colors[i], 2,
                    colors[i]);
   }

   std::vector<std::string> legEntries;
   if (pr0->ffileName.find("dog1") != std::string::npos or
       pr0->ffileName.find("cosmics") != std::string::npos) {
      legEntries.push_back("Cosmic rays");
   } else if (pr0->ffileName.find("beam") != std::string::npos or
              pr0->ffileName.find("sandmu") != std::string::npos) {
      legEntries.push_back("Beam data");
   }
   legEntries.push_back("Simulation");

   // -----------------------------------------------------------------------------------------------------------------
   // Bethe Bloch fit on data
   const char *formula =
   "[0]/pow(x/sqrt(x*x+[5]*[5]),[3]) * ( [1] - pow(x/sqrt(x*x+[5]*[5]),[3]) - log([2] "
   "+ 1 / pow(x/[5], [4])) )";
   // Muon
   TF1 *bethefixmuon = new TF1("bethefixmuon", formula, 0, 2000, "");
   bethefixmuon->FixParameter(1, 3.62857e+00);
   bethefixmuon->FixParameter(2, 3.18209e-02);
   bethefixmuon->FixParameter(3, 2.07081e+00);
   bethefixmuon->FixParameter(4, 7.14413e-01);
   bethefixmuon->FixParameter(5, 105.658);
   pr0->ptge_absmom_mean_XP->Fit(bethefixmuon, "R", "", 100, 2000);
   pr0->ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetLineColor(kBlack);
   pr0->ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetLineWidth(2);
   pr0->ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetRange(0, 2000);
   // Proton
   TF1 *bethefixproton = (TF1 *)bethefixmuon->Clone("bethefixproton");
   bethefixproton->SetParameter(5, 938.272);
   // Electron
   TF1 *bethefixelectron = (TF1 *)bethefixmuon->Clone("bethefixelectron");
   bethefixelectron->SetParameter(5, 0.511);

   // -----------------------------------------------------------------------------------------------------------------
   // Global dE/dx plot
   pr0->fph1f_XP->SetAxisRange(0, 1.1 * ampmax, "Y");
   // for (int i = 1; i < ncomparisons; i++)
   //    v_processes[i]->fph1f_XP->Scale(pr0->fph1f_XP->Integral() /
   //                                    v_processes[i]->fph1f_XP->Integral());
   float maxcount =
      std::max(pr0->fph1f_XP->GetMaximum(), v_processes[1]->fph1f_XP->GetMaximum());
   std::vector<std::string> legEntriesDistrib;
   legEntriesDistrib.push_back(legEntries[0]);
   legEntriesDistrib.push_back("MC");
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->fph1f_XP->SetAxisRange(0, 1.1 * maxcount, "Y");
      v_processes[i]->fph1f_XP->Draw(i == 0 ? "HIST" : "HIST same");
      float xwidth = 0.23;
      float xright = 0.94;
      float ywidth = (0.92 - 0.15) / ncomparisons - 0.01;
      float ytop = 0.92 - i * (ywidth + 0.01);
      std::cout << i << " " << xright - xwidth << " < x < " << xright << " | "
                << ytop - ywidth << " < y < " << ytop << std::endl;
      PrintResolution(v_processes[i]->fph1f_XP, fpCanvas, xright, ytop, xwidth, ywidth,
                      "north east", colors[i], legEntriesDistrib[i]);
   }
   fpCanvas->SaveAs((OutputFile + "(").c_str());
   fpCanvas->Clear();

   // dE/dx vs momentum
   // -----------------------------------------------------------------------------------------------------------------
   // Resolution
   pr0->ptge_absmom_reso_XP->SetTitle(";Momentum (MeV/c);dE/dx resolution (%)");
   pr0->ptge_absmom_reso_XP->GetXaxis()->SetLimits(0, momrange * 1.05);
   pr0->ptge_absmom_reso_XP->GetYaxis()->SetRangeUser(resominmom, resomaxmom);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absmom_reso_XP->SetMarkerSize(7);
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_absmom_reso_XP->SetTitle(";Momentum (MeV/c);dE/dx resolution (%)");
   pr0->ptge_absmom_reso_XP->GetXaxis()->SetLimits(0, momrange * 1.05);
   pr0->ptge_absmom_reso_XP->GetYaxis()->SetRangeUser(resominmom, resomaxmom);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_reso_XP->SetMarkerSize(3);
      v_processes[i]->ptge_absmom_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absmom_reso_WF->DrawClone("P same");
      v_processes[i]->ptge_absmom_reso_XP->SetMarkerSize(7);
      v_processes[i]->ptge_absmom_reso_WF->SetMarkerSize(7);
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // Mean
   fpLegend->Clear();
   pr0->ptge_absmom_mean_XP->SetTitle(";Momentum (MeV/c);dE/dx (ADC counts/cm)");
   pr0->ptge_absmom_mean_XP->GetXaxis()->SetLimits(0, momrange * 1.05);
   pr0->ptge_absmom_mean_XP->GetYaxis()->SetRangeUser(dedxminmom, dedxmaxmom);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_absmom_mean_XP->SetTitle(";Momentum (MeV/c);dE/dx (ADC counts/cm)");
   pr0->ptge_absmom_mean_XP->GetXaxis()->SetLimits(0, momrange * 1.05);
   pr0->ptge_absmom_mean_XP->GetYaxis()->SetRangeUser(dedxminmom, dedxmaxmom);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absmom_mean_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // Heatmap
   gPad->SetRightMargin(0.11);
   pr0->fph2f_XPabsmom->GetYaxis()->SetRangeUser(200, 1000);
   pr0->fph2f_XPabsmom->Draw("colz");
   bethefixmuon->SetLineColor(kCyan);
   bethefixmuon->Draw("same");
   bethefixproton->SetLineColor(kRed);
   bethefixproton->Draw("same");
   bethefixelectron->SetLineColor(kGreen);
   bethefixelectron->Draw("same");
   TLegend legbethe(0.6, 0.7, 0.83, 0.9);
   legbethe.SetFillColorAlpha(kWhite, 0.8);
   legbethe.SetTextSize(0.06);
   legbethe.AddEntry(bethefixmuon, "Muon fit", "l");
   legbethe.AddEntry(bethefixproton, "Proton fit", "l");
   legbethe.AddEntry(bethefixelectron, "Electron fit", "l");
   legbethe.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   delete bethefixmuon;
   delete bethefixproton;
   delete bethefixelectron;
   fpCanvas->Clear();
   gPad->SetRightMargin(0.035);

   // dE/dx vs drift distance
   // -----------------------------------------------------------------------------------------------------------------
   // Resolution
   fpLegend->Clear();
   pr0->ptge_dd_reso_XP->SetTitle(";Drift distance (mm);dE/dx resolution (%)");
   pr0->ptge_dd_reso_XP->GetXaxis()->SetLimits(-50, 1050);
   pr0->ptge_dd_reso_XP->GetYaxis()->SetRangeUser(resominX, resomaxX);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_dd_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_dd_reso_XP->SetTitle(";Drift distance (mm);dE/dx resolution (%)");
   pr0->ptge_dd_reso_XP->GetXaxis()->SetLimits(-50, 1050);
   pr0->ptge_dd_reso_XP->GetYaxis()->SetRangeUser(resominX, resomaxX);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_dd_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_dd_reso_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // Mean
   fpLegend->Clear();
   pr0->ptge_dd_mean_XP->SetTitle(";Drift distance (mm);dE/dx (ADC counts/cm)");
   pr0->ptge_dd_mean_XP->GetXaxis()->SetLimits(-50, 1050);
   pr0->ptge_dd_mean_XP->GetYaxis()->SetRangeUser(dedxminX, dedxmaxX);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_dd_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_dd_mean_XP->SetTitle(";Drift distance (mm);dE/dx (ADC counts/cm)");
   pr0->ptge_dd_mean_XP->GetXaxis()->SetLimits(-50, 1050);
   pr0->ptge_dd_mean_XP->GetYaxis()->SetRangeUser(dedxminX, dedxmaxX);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_dd_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_dd_mean_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx VS track length
   // Resolution
   fpLegend->Clear();
   pr0->ptge_trklen_reso_XP->SetTitle(";Track length (mm);dE/dx resolution (%)");
   pr0->ptge_trklen_reso_XP->GetXaxis()->SetLimits(-50, 1850);
   pr0->ptge_trklen_reso_XP->GetYaxis()->SetRangeUser(resomintrklen, resomaxtrklen);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_trklen_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_trklen_reso_XP->SetTitle(";Track length (mm);dE/dx resolution (%)");
   pr0->ptge_trklen_reso_XP->GetXaxis()->SetLimits(-50, 1850);
   pr0->ptge_trklen_reso_XP->GetYaxis()->SetRangeUser(resomintrklen, resomaxtrklen);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_trklen_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_trklen_reso_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // Mean
   fpLegend->Clear();
   pr0->ptge_trklen_mean_XP->SetTitle(";Track length (mm);dE/dx (ADC counts/cm)");
   pr0->ptge_trklen_mean_XP->GetXaxis()->SetLimits(-50, 1850);
   pr0->ptge_trklen_mean_XP->GetYaxis()->SetRangeUser(dedxmintrklen, dedxmaxtrklen);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_trklen_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_trklen_mean_XP->SetTitle(";Track length (mm);dE/dx (ADC counts/cm)");
   pr0->ptge_trklen_mean_XP->GetXaxis()->SetLimits(-50, 1850);
   pr0->ptge_trklen_mean_XP->GetYaxis()->SetRangeUser(dedxmintrklen, dedxmaxtrklen);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_trklen_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_trklen_mean_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx VS phi bins
   // Resolution
   fpLegend->Clear();
   pr0->ptge_absphi_reso_XP->SetTitle(";#varphi (#circ);dE/dx resolution (%)");
   pr0->ptge_absphi_reso_XP->GetXaxis()->SetLimits(-3, 90);
   pr0->ptge_absphi_reso_XP->GetYaxis()->SetRangeUser(resominphi, resomaxphi);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absphi_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_absphi_reso_XP->SetTitle(";#varphi (#circ);dE/dx resolution (%)");
   pr0->ptge_absphi_reso_XP->GetXaxis()->SetLimits(-3, 90);
   pr0->ptge_absphi_reso_XP->GetYaxis()->SetRangeUser(resominphi, resomaxphi);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absphi_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absphi_reso_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // Mean
   fpLegend->Clear();
   pr0->ptge_absphi_mean_XP->SetTitle(";#varphi (#circ);dE/dx (ADC counts/cm)");
   pr0->ptge_absphi_mean_XP->GetXaxis()->SetLimits(-3, 90);
   pr0->ptge_absphi_mean_XP->GetYaxis()->SetRangeUser(dedxminphi, dedxmaxphi);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absphi_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_absphi_mean_XP->SetTitle(";#varphi (#circ);dE/dx (ADC counts/cm)");
   pr0->ptge_absphi_mean_XP->GetXaxis()->SetLimits(-3, 90);
   pr0->ptge_absphi_mean_XP->GetYaxis()->SetRangeUser(dedxminphi, dedxmaxphi);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absphi_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absphi_mean_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx VS theta bins
   // Resolution
   fpLegend->Clear();
   pr0->ptge_abstheta_reso_XP->SetTitle(";#theta (#circ);dE/dx resolution (%)");
   pr0->ptge_abstheta_reso_XP->GetXaxis()->SetLimits(-3, 90);
   pr0->ptge_abstheta_reso_XP->GetYaxis()->SetRangeUser(resomintheta, resomaxtheta);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_abstheta_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_abstheta_reso_XP->SetTitle(";#theta (#circ);dE/dx resolution (%)");
   pr0->ptge_abstheta_reso_XP->GetXaxis()->SetLimits(-3, 90);
   pr0->ptge_abstheta_reso_XP->GetYaxis()->SetRangeUser(resomintheta, resomaxtheta);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_abstheta_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_abstheta_reso_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // Mean
   fpLegend->Clear();
   pr0->ptge_abstheta_mean_XP->SetTitle(";#theta (#circ);dE/dx (ADC counts/cm)");
   pr0->ptge_abstheta_mean_XP->GetXaxis()->SetLimits(-3, 90);
   pr0->ptge_abstheta_mean_XP->GetYaxis()->SetRangeUser(dedxmintheta, dedxmaxtheta);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_abstheta_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP, legEntries[i].c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();
   // Add WF
   fpLegend->Clear();
   pr0->ptge_abstheta_mean_XP->SetTitle(";#theta (#circ);dE/dx (ADC counts/cm)");
   pr0->ptge_abstheta_mean_XP->GetXaxis()->SetLimits(-3, 90);
   pr0->ptge_abstheta_mean_XP->GetYaxis()->SetRangeUser(dedxmintheta, dedxmaxtheta);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_abstheta_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_abstheta_mean_WF->DrawClone("P same");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_XP,
                         ("XP " + legEntries[i]).c_str(), "p");
      fpLegend->AddEntry(v_processes[i]->ptge_absmom_reso_WF,
                         ("WF " + legEntries[i]).c_str(), "p");
   }
   fpLegend->Draw();
   fpCanvas->SaveAs((OutputFile + ")").c_str());

   for (Process *pr : v_processes)
      delete pr;
}
