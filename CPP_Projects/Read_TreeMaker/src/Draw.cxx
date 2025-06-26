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
   int nabsmommeanbins = pr.nabsmommeanbins;
   int nmombins = pr.nmombins;
   fpCanvas->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   fpLegend = new TLegend(0.8, 0.8, 0.98, 0.98);
   fpLegend->SetTextSize(0.08);
   fpLegend->SetFillStyle(1001);
   fpLegend->SetTextColor(kBlue - 1);
   fpLegend->SetBorderSize(2);
   fpLegend->SetLineColor(kBlue - 1);
   fpLegend->SetTextAlign(22);
   fpLegend->SetMargin(0.4);
   TLegendEntry *pentryWF = fpLegend->AddEntry((TObject *)0, "WF", "p");
   pentryWF->SetMarkerColor(kCyan + 2);
   pentryWF->SetMarkerSize(9);
   pentryWF->SetMarkerStyle(33);
   TLegendEntry *pentryXP = fpLegend->AddEntry((TObject *)0, "XP", "p");
   pentryXP->SetMarkerColor(kMagenta + 2);
   pentryXP->SetMarkerSize(9);
   pentryXP->SetMarkerStyle(47);
   gStyle->SetTitleSize(0.07, "xyz");
   gStyle->SetLabelSize(0.06, "xyz");
   // gStyle->SetNdivisions(505, "xyz");
   // gStyle->SetTitleOffset(0.8, "y");
   SetMarginH1((TPad *)gPad);
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
   PrintResolution(pr.fph1f_XP, fpCanvas, 0.65 - invX, 0.6, kMagenta + 2, "XP");
   PrintResolution(pr.fph1f_WF, fpCanvas, 0.65 - invX, 0.27, kCyan + 2, "WF");
   fpCanvas->SaveAs((fRealpathPDF + "(").c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx ERAM by ERAM
   fpCanvas->SetCanvasSize(fpCanvas->GetWw(), 1200);
   float maxdEdx = 0;
   for (TH1F *hist : pr.vmod_fph1f_WF)
      if (hist->GetMaximum() > maxdEdx)
         maxdEdx = hist->GetMaximum();
   for (TH1F *hist : pr.vmod_fph1f_XP)
      if (hist->GetMaximum() > maxdEdx)
         maxdEdx = hist->GetMaximum();

   float xMax = 0, yMax = 0, reso = 0, dreso = 0;
   fpCanvas->Clear();
   fpCanvas->Divide(4, 4, 0.0, 0.0);
   for (int i = 0; i < 32; i++) {
      pr.vmod_fph1f_WF[i]->SetAxisRange(0, 1.1 * maxdEdx, "Y");
      pr.vmod_fph1f_XP[i]->SetAxisRange(0, 1.1 * maxdEdx, "Y");
      Graphic_setup(pr.vmod_fph1f_WF[i], 0.5, 1, kCyan + 1, 1, kCyan - 2, kCyan, 0.2);
      Graphic_setup(pr.vmod_fph1f_XP[i], 0.5, 1, kMagenta + 2, 1, kMagenta - 2, kMagenta,
                    0.2);
   }

   std::vector<int> vSouth = {16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7};
   std::vector<int> vNorth = {24, 25, 26, 27, 28, 29, 30, 31,
                              8,  9,  10, 11, 12, 13, 14, 15};

   for (std::vector<int> vEndplate : {vSouth, vNorth}) {
      for (int i = 0; i < 16; i++) {
         fpCanvas->cd(i + 1);
         int eramIndex = vEndplate[i];
         int row = i / 4;
         int col = i % 4;
         // margins: left, right, bottom, top
         double left = 0, bottom = 0, myFontSize = 0.1;
         double leftPave = 0.05, rightPave = 0.35, bottomPave = 0.83, topPave = 0.93;
         double bottomResoShift = 0, leftResoShift = 0;
         if (col == 0) {
            left = 0.15;
            leftPave += 0.13;
            rightPave += 0.1;
            leftResoShift = 0.1;
         }
         if (row == 3) {
            bottom = 0.15;
            myFontSize -= 0.012;
            bottomResoShift = 0.1;
            bottomPave += 0.03;
            topPave += 0.025;
         }
         gPad->SetMargin(left, 0, bottom, 0);

         TPaveText pavetextMod(leftPave, bottomPave, rightPave, topPave, "NDC");
         pavetextMod.SetTextAlign(22);
         pavetextMod.SetTextSize(myFontSize);
         pavetextMod.SetFillStyle(1001);
         pavetextMod.SetFillColor(kWhite);
         pavetextMod.SetBorderSize(1);
         pavetextMod.AddText(pr.vmod_fph1f_WF[eramIndex]->GetTitle());

         pr.vmod_fph1f_WF[eramIndex]->SetTitle("");
         pr.vmod_fph1f_WF[eramIndex]->Draw("HIST");
         pr.vmod_fph1f_XP[eramIndex]->Draw("HIST same");
         if (pr.vmod_fph1f_WF[eramIndex]->GetEntries() < 100) {
            pavetextMod.DrawClone();
            continue;
         }
         xMax = pr.vmod_fph1f_WF[eramIndex]->GetXaxis()->GetXmax();
         yMax = pr.vmod_fph1f_WF[eramIndex]->GetMaximum();
         pr.vmod_fph1f_WF[eramIndex]->GetMean() > dEdxmax / 2 ? invX = 0.4 : invX = 0;
         float x0 = 0.43, y0 = 0.55, xwidth = x0 + 0.1;
         PrintResolution(pr.vmod_fph1f_XP[eramIndex], fpCanvas, x0 + leftResoShift - invX,
                         y0 + 0.5 * bottomResoShift, xwidth - 0.8 * leftResoShift,
                         0.4 - 0.5 * bottomResoShift, "south west", kMagenta + 2, "", 1,
                         myFontSize);
         PrintResolution(pr.vmod_fph1f_WF[eramIndex], fpCanvas, x0 + leftResoShift - invX,
                         0.1 + bottomResoShift, xwidth - 0.8 * leftResoShift,
                         0.4 - 0.5 * bottomResoShift, "south west", kCyan + 2, "", 1,
                         myFontSize);
         pavetextMod.DrawClone();
      }
      fpCanvas->SaveAs(fRealpathPDF.c_str());
   }
   fpCanvas->SetCanvasSize(fpCanvas->GetWw(), 1350);

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx XP vs WF
   fpCanvas->Clear();
   SetMarginH2((TPad *)gPad);
   gStyle->SetOptStat("merou");
   gStyle->SetStatX(0.36);
   gStyle->SetStatY(0.91);
   SetMarginH2((TPad *)gPad);
   pr.fph2f_WFXP->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dEdx vs Momentum
   // Resolution
   fpCanvas->Clear();
   SetMarginH1((TPad *)gPad);
   Graphic_setup(pr.ptge_absmom_reso_WF, 5, 33, kCyan + 2, 0, kCyan + 2);
   Graphic_setup(pr.ptge_absmom_reso_XP, 5, 47, kMagenta + 2, 0, kMagenta + 2);
   pr.ptge_absmom_reso_XP->SetTitle(";Momentum [MeV/c];dE/dx resolution [%]");
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
   pr.ptge_absmom_mean_XP->SetTitle(";Momentum [MeV/c];Mean dE/dx [ADC counts/cm]");
   pr.ptge_absmom_mean_XP->GetXaxis()->SetLimits(0, absmomrange);
   pr.ptge_absmom_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_absmom_mean_XP->Draw("APL");
   pr.ptge_absmom_mean_WF->Draw("PL same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // 2D distribution
   // WF
   fpCanvas->Clear();
   SetMarginH2((TPad *)gPad);
   gStyle->SetOptStat(0);
   pr.fph2f_WFabsmommean->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

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
   // Fit on MC
   pr.ptge_absmom_mean_XP->Fit(bethefixmuon, "R", "", 100, 2000);
   pr.ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetLineColor(kBlack);
   pr.ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetLineWidth(2);
   pr.ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetRange(0, 2000);
   // Electron
   TF1 *bethefixelectron = (TF1 *)bethefixmuon->Clone("bethefixelectron");
   bethefixelectron->SetParameter(5, 0.511);
   // Proton
   TF1 *bethefixproton = (TF1 *)bethefixmuon->Clone("bethefixproton");
   bethefixproton->SetParameter(5, 938.272);

   // XP
   fpCanvas->Clear();
   pr.fph2f_XPabsmommean->Draw("colz");
   bethefixmuon->SetLineColor(kCyan);
   bethefixmuon->SetLineStyle(9);
   bethefixmuon->Draw("same");
   bethefixelectron->SetLineColor(kOrange + 1);
   bethefixelectron->SetLineStyle(9);
   bethefixelectron->Draw("same");
   bethefixproton->SetLineColor(kRed);
   bethefixproton->SetLineStyle(9);
   TLegend legbethe(0.53, 0.76, 0.83, 0.91);
   legbethe.SetTextSize(0.06);
   legbethe.SetBorderSize(2);
   legbethe.SetLineColor(kBlue - 1);
   legbethe.AddEntry(bethefixmuon, "Muon fit", "l");
   legbethe.AddEntry(bethefixelectron, "Electron fit", "l");
   if (pr.frun.find("beam") != std::string::npos) {
      legbethe.SetY1(0.7);
      bethefixproton->Draw("same");
      legbethe.AddEntry(bethefixproton, "Proton fit", "l");
   }
   legbethe.Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());
   delete bethefixmuon;
   delete bethefixelectron;

   // -----------------------------------------------------------------------------------------------------------------
   // dEdx vs momentum
   // Resolution
   fpCanvas->Clear();
   SetMarginH1((TPad *)gPad);
   Graphic_setup(pr.ptge_mom_reso_WF, 3, 33, kCyan + 2, 1, kCyan + 2);
   Graphic_setup(pr.ptge_mom_reso_XP, 3, 47, kMagenta + 2, 1, kMagenta + 2);
   pr.ptge_mom_reso_XP->SetTitle(";Momentum [MeV/c];dE/dx resolution [%]");
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
   pr.ptge_mom_mean_XP->SetTitle(";Momentum [MeV/c];Mean dE/dx [ADC counts/cm]");
   pr.ptge_mom_mean_XP->GetXaxis()->SetLimits(-momrange, momrange);
   pr.ptge_mom_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_mom_mean_XP->Draw("APL");
   pr.ptge_mom_mean_WF->Draw("PL same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution
   fpCanvas->Clear();
   SetMarginH2((TPad *)gPad);
   pr.fph2f_WFmom->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_XPmom->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   SetMarginH1((TPad *)gPad);
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
   SetMarginH2((TPad *)gPad);
   pr.fph2f_chi2mom->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs X position
   // Resolution
   fpCanvas->Clear();
   fpCanvas->SetGrid(1, 1);
   SetMarginH1((TPad *)gPad);
   Graphic_setup(pr.ptge_X_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_X_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_X_reso_XP->SetTitle(";Track X position [mm];dE/dx resolution [%]");
   pr.ptge_X_reso_XP->GetXaxis()->SetLimits(-1000, 1000);
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
   pr.ptge_X_mean_XP->SetTitle(";Track X position [mm];dE/dx [ADC counts/cm]");
   pr.ptge_X_mean_XP->GetXaxis()->SetLimits(-1000, 1000);
   pr.ptge_X_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_X_mean_XP->Draw("AP");
   pr.ptge_X_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution WF vs X
   fpCanvas->Clear();
   SetMarginH2((TPad *)gPad);
   pr.fph2f_WFX->SetAxisRange(-1000, 1000, "X");
   pr.fph2f_WFX->Draw("colz");
   TGraphErrors *ptge_mom_WFX = Convert_TH2_TGE(pr.fph2f_WFX);
   TF1 linearFitWFX("linearFitWFX", "pol1", 0, 1000);
   ptge_mom_WFX->Fit(&linearFitWFX, "RQ");
   linearFitWFX.SetLineColor(kRed);
   // linearFitWFX.Draw("same");
   TPaveText pavetext(0.4, 0.75, 0.8, 0.9, "NDC");
   pavetext.AddText(Form("y = %.2fx + %.2f", linearFitWFX.GetParameter(1),
                         linearFitWFX.GetParameter(0)));
   pavetext.SetFillColorAlpha(kWhite, 0.9);
   pavetext.SetBorderSize(1);
   pavetext.SetLineColor(kBlue - 1);
   pavetext.SetTextColor(kRed);
   pavetext.SetTextSize(0.07);
   // pavetext.DrawClone();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution XP vs X
   pr.fph2f_XPX->SetAxisRange(-1000, 1000, "X");
   pr.fph2f_XPX->Draw("colz");
   TGraphErrors *ptge_mom_XPX = Convert_TH2_TGE(pr.fph2f_XPX);
   TF1 linearFitXPX("linearFitXPX", "pol1", 0, 1000);
   ptge_mom_XPX->Fit(&linearFitXPX, "RQ");
   linearFitXPX.SetLineColor(kRed);
   // linearFitXPX.Draw("same");
   pavetext.Clear();
   pavetext.AddText(Form("y = %.2fx + %.2f", linearFitXPX.GetParameter(1),
                         linearFitXPX.GetParameter(0)));
   // pavetext.DrawClone();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs Track length
   // Resolution
   fpCanvas->Clear();
   fpCanvas->SetGrid(1, 1);
   SetMarginH1((TPad *)gPad);
   Graphic_setup(pr.ptge_trklen_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_trklen_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_trklen_reso_XP->SetTitle(";Track length [mm];dE/dx resolution [%]");
   pr.ptge_trklen_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   pr.ptge_trklen_reso_XP->DrawClone("AP");
   pr.ptge_trklen_reso_WF->DrawClone("P same");
   pr.ptge_trklen_reso_XP->SetMarkerSize(7);
   pr.ptge_trklen_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Mean
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_trklen_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_trklen_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_trklen_mean_XP->SetTitle(";Track length [mm];Mean dE/dx [ADC counts/cm]");
   pr.ptge_trklen_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_trklen_mean_XP->Draw("AP");
   pr.ptge_trklen_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distributions
   SetMarginH2((TPad *)gPad);
   gStyle->SetStatX(0.87);
   pr.fph2f_XPlen->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   pr.fph2f_momlen->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   pr.fph2f_momncl->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs absolute Track angle phi
   // Resolution
   fpCanvas->Clear();
   SetMarginH1((TPad *)gPad);
   Graphic_setup(pr.ptge_absphi_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_absphi_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_absphi_reso_XP->SetTitle(";Track angle #varphi [#circ];dE/dx resolution [%]");
   pr.ptge_absphi_reso_XP->GetXaxis()->SetLimits(0, 93);
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
   pr.ptge_absphi_mean_XP->SetTitle(
      ";Track angle #varphi [#circ];Mean dE/dx [ADC counts/cm]");
   pr.ptge_absphi_mean_XP->GetXaxis()->SetLimits(0, 93);
   pr.ptge_absphi_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_absphi_mean_XP->Draw("AP");
   pr.ptge_absphi_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs Track angle phi
   // Resolution
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_phi_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_phi_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_phi_reso_XP->SetTitle(";Track angle #varphi [#circ];dE/dx resolution [%]");
   pr.ptge_phi_reso_XP->GetXaxis()->SetLimits(-93, 93);
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
   pr.ptge_phi_mean_XP->SetTitle(
      ";Track angle #varphi [#circ];Mean dE/dx [ADC counts/cm]");
   pr.ptge_phi_mean_XP->GetXaxis()->SetLimits(-93, 93);
   pr.ptge_phi_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_phi_mean_XP->Draw("AP");
   pr.ptge_phi_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution
   fpCanvas->Clear();
   SetMarginH2((TPad *)gPad);
   pr.fph2f_XPphi->SetTitle(";Track angle #varphi [#circ];dE/dx with XP [ADC counts/cm]");
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
   SetMarginH1((TPad *)gPad);
   pr.fph1f_phi->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // dE/dx vs absolute Track angle theta
   // -----------------------------------------------------------------------------------------------------------------
   // Resolution
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_abstheta_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_abstheta_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_abstheta_reso_XP->SetTitle(";Track angle #theta [#circ];dE/dx resolution [%]");
   pr.ptge_abstheta_reso_XP->GetXaxis()->SetLimits(0, 93);
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
   pr.ptge_abstheta_mean_XP->SetTitle(
      ";Track angle #theta [#circ];Mean dE/dx [ADC counts/cm]");
   pr.ptge_abstheta_mean_XP->GetXaxis()->SetLimits(0, 93);
   pr.ptge_abstheta_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_abstheta_mean_XP->Draw("AP");
   pr.ptge_abstheta_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // dE/dx vs Track angle theta
   // -----------------------------------------------------------------------------------------------------------------
   // Resolution
   fpCanvas->Clear();
   Graphic_setup(pr.ptge_theta_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(pr.ptge_theta_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   pr.ptge_theta_reso_XP->SetTitle(";Track angle #theta [#circ];dE/dx resolution [%]");
   pr.ptge_theta_reso_XP->GetXaxis()->SetLimits(-93, 93);
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
   pr.ptge_theta_mean_XP->SetTitle(
      ";Track angle #theta [#circ];Mean dE/dx [ADC counts/cm]");
   pr.ptge_theta_mean_XP->GetXaxis()->SetLimits(-93, 93);
   pr.ptge_theta_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   pr.ptge_theta_mean_XP->Draw("AP");
   pr.ptge_theta_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // 2D distribution
   fpCanvas->Clear();
   SetMarginH2((TPad *)gPad);
   pr.fph2f_XPtheta->SetTitle(
      ";Track angle #theta [#circ];dE/dx with XP [ADC counts/cm]");
   pr.fph2f_XPtheta->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_lentheta->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_momtheta->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   SetMarginH1((TPad *)gPad);
   pr.fph1f_theta->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // Both angles
   fpCanvas->Clear();
   SetMarginH2((TPad *)gPad);
   pr.fph2f_phitheta->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // Pulls
   // -----------------------------------------------------------------------------------------------------------------
   // Pulls vs theta
   // Standard deviation
   fpCanvas->Clear();
   SetMarginH1((TPad *)gPad);
   Graphic_setup(pr.ptge_theta_std_pullmu, 2, 21, kBlue + 1, 2, kBlue + 1);
   Graphic_setup(pr.ptge_theta_std_pullelec, 2, 20, kOrange + 7, 2, kOrange + 7);
   pr.ptge_theta_std_pullmu->SetTitle(
      ";Track angle #theta [#circ];Pulls standard deviation");
   float minpull = 0;
   float maxpull = (pr.ptge_theta_std_pullmu->GetHistogram()->GetMaximum() >
                    pr.ptge_theta_std_pullelec->GetHistogram()->GetMaximum())
                      ? pr.ptge_theta_std_pullmu->GetHistogram()->GetMaximum()
                      : pr.ptge_theta_std_pullelec->GetHistogram()->GetMaximum();
   float diffpull = maxpull - minpull;
   pr.ptge_theta_std_pullmu->GetXaxis()->SetLimits(-93, 93);
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
   pr.ptge_theta_mean_pullmu->SetTitle(";Track angle #theta [#circ];Pulls mean");
   minpull = (pr.ptge_theta_mean_pullmu->GetHistogram()->GetMinimum() <
              pr.ptge_theta_mean_pullelec->GetHistogram()->GetMinimum())
                ? pr.ptge_theta_mean_pullmu->GetHistogram()->GetMinimum()
                : pr.ptge_theta_mean_pullelec->GetHistogram()->GetMinimum();
   maxpull = (pr.ptge_theta_mean_pullmu->GetHistogram()->GetMaximum() >
              pr.ptge_theta_mean_pullelec->GetHistogram()->GetMaximum())
                ? pr.ptge_theta_mean_pullmu->GetHistogram()->GetMaximum()
                : pr.ptge_theta_mean_pullelec->GetHistogram()->GetMaximum();
   pr.ptge_theta_mean_pullmu->GetXaxis()->SetLimits(-93, 93);
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
   Graphic_setup(pr.ptge_mom_std_pullmu, 2, 21, kBlue + 1, 2, kBlue + 1);
   Graphic_setup(pr.ptge_mom_std_pullelec, 2, 20, kOrange + 7, 2, kOrange + 7);
   pr.ptge_mom_std_pullmu->SetTitle(";Momentum [MeV/c];Pulls standard deviation");
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
   pr.ptge_mom_mean_pullmu->SetTitle(";Momentum [MeV/c];Pulls mean");
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
   SetMarginH2((TPad *)gPad);
   pr.fph2f_pullelecmu->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Pulls distribution
   fpCanvas->Clear();
   SetMarginH1((TPad *)gPad);
   Graphic_setup(pr.fph1f_pullmu, 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan, 0.2);
   Graphic_setup(pr.fph1f_pullelec, 0.5, 1, kOrange + 7, 2, kOrange - 2, kOrange, 0.2);
   Graphic_setup(pr.fph1f_pullproton, 0.5, 1, kRed + 1, 2, kRed - 2, kRed, 0.2);
   TLegend legpull2(0.2, 0.7, 0.35, 0.9);
   legpull2.AddEntry(pr.fph1f_pullmu, "#mu?", "f");
   legpull2.AddEntry(pr.fph1f_pullelec, "e?", "f");
   legpull2.AddEntry(pr.fph1f_pullproton, "p?", "f");
   legpull2.SetTextSize(0.06);
   legpull2.SetTextColor(kBlue - 1);
   legpull2.SetLineColor(kBlue - 1);
   legpull2.SetBorderSize(1);
   legpull2.SetTextAlign(22);
   legpull2.SetMargin(0.5);
   pr.fph1f_pullmu->SetXTitle("Pull");
   float maxpull2 = std::max(
      pr.fph1f_pullmu->GetMaximum(),
      std::max(pr.fph1f_pullelec->GetMaximum(), pr.fph1f_pullproton->GetMaximum()));
   pr.fph1f_pullmu->GetYaxis()->SetRangeUser(0, maxpull2 + 0.1 * maxpull2);
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

   TPaveText pavetextmu(0.6, 0.75, 0.93, 0.9, "NDC");
   pavetextmu.SetTextSize(0.05);
   pavetextmu.SetTextColor(kCyan + 4);
   pavetextmu.SetFillColorAlpha(kCyan, 0.6);
   pavetextmu.SetLineWidth(2);
   pavetextmu.SetBorderSize(2);
   pavetextmu.SetTextAlign(22);
   pavetextmu.SetLineColor(kCyan - 2);
   pavetextmu.AddText(Form("#mu_{#mu} = %.3f #pm %.3f", meanmu, dmeanmu));
   pavetextmu.AddText(Form("#sigma_{#mu} = %.3f #pm %.3f", stdmu, dstdmu));

   TPaveText pavetextelec(0.6, 0.55, 0.93, 0.7, "NDC");
   pavetextelec.SetTextSize(0.05);
   pavetextelec.SetTextColor(kOrange + 3);
   pavetextelec.SetFillColorAlpha(kOrange, 0.6);
   pavetextelec.SetLineWidth(2);
   pavetextelec.SetBorderSize(2);
   pavetextelec.SetTextAlign(22);
   pavetextelec.SetLineColor(kOrange - 1);
   pavetextelec.AddText(Form("#mu_{e} = %.3f #pm %.3f", meanelec, dmeanelec));
   pavetextelec.AddText(Form("#sigma_{e} = %.3f #pm %.3f", stdelec, dstdelec));

   TPaveText paveTextSep(0.6, 0.4, 0.93, 0.5, "NDC");
   paveTextSep.SetTextSize(0.058);
   paveTextSep.SetTextColor(kBlue - 1);
   paveTextSep.SetFillColor(kWhite);
   paveTextSep.SetLineWidth(2);
   paveTextSep.SetBorderSize(2);
   paveTextSep.SetTextAlign(22);
   paveTextSep.SetLineColor(kRed);
   paveTextSep.AddText(Form("S(#mu/e) = %.2f #pm %.2f", separation, dseparation));

   legpull2.Draw();
   pavetextmu.Draw();
   pavetextelec.Draw();
   paveTextSep.Draw();
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // Geometry plots
   // -----------------------------------------------------------------------------------------------------------------
   pr.fph1f_trklen->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   pr.fph1f_chi2->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   pr.fph1f_Xpos->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   pr.fph1f_dir0->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   pr.fph1f_dir1->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   pr.fph1f_dir2->Draw("HIST");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   SetMarginH2((TPad *)gPad);
   gPad->SetLogz();
   pr.fph2f_XZ->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   pr.fph2f_YZ->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());
   gPad->SetLogz(0);

   // -----------------------------------------------------------------------------------------------------------------
   // Debugging

   fpCanvas->Clear();
   pr.fph2f_momR->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   fpCanvas->Clear();
   pr.fph2f_chi2ndfR->Draw("colz");
   fpCanvas->SaveAs(fRealpathPDF.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   TPaveText pavetextBins(0.65, 0.86, 0.98, 0.94, "NDC");
   pavetextBins.SetTextAlign(22);
   pavetextBins.SetTextSize(0.05);
   pavetextBins.SetFillStyle(1001);
   pavetextBins.SetFillColor(kWhite);
   pavetextBins.SetBorderSize(2);
   pavetextBins.SetLineWidth(2);

   // dE/dx VS momentum bin
   SetMarginH1((TPad *)gPad);
   gPad->SetTopMargin(0.03);
   for (int i = 0; i < pr.nabsmomresobins; i++) {
      fpCanvas->Clear();
      float invX = 0;
      if (pr.vabsmomreso_fph1f_WF[i]->GetMean() > dEdxmax / 2)
         invX = 0.4;
      pavetextBins.Clear();
      pavetextBins.SetX1NDC(0.65 - invX);
      pavetextBins.SetX2NDC(0.98 - invX);
      pavetextBins.AddText(pr.vabsmomreso_fph1f_WF[i]->GetTitle());
      pr.vabsmomreso_fph1f_WF[i]->SetTitle("");
      pr.vabsmomreso_fph1f_WF[i]->SetMaximum(
         1.1 * std::max({pr.vabsmomreso_fph1f_WF[i]->GetMaximum(),
                         pr.vabsmomreso_fph1f_XP[i]->GetMaximum()}));
      Graphic_setup(pr.vabsmomreso_fph1f_WF[i], 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan,
                    0.2);
      Graphic_setup(pr.vabsmomreso_fph1f_XP[i], 0.5, 1, kMagenta + 1, 2, kMagenta - 2,
                    kMagenta, 0.2);
      pr.vabsmomreso_fph1f_WF[i]->Draw("HIST");
      pr.vabsmomreso_fph1f_XP[i]->Draw("HIST same");
      PrintResolution(pr.vabsmomreso_fph1f_XP[i], fpCanvas, 0.65 - invX, 0.53,
                      kMagenta + 2, "XP");
      PrintResolution(pr.vabsmomreso_fph1f_WF[i], fpCanvas, 0.65 - invX, 0.2, kCyan + 2,
                      "WF");
      pavetextBins.Draw();
      fpCanvas->SaveAs(fRealpathPDF.c_str());
   }

   // dE/dx VS signed momentum bin
   // -----------------------------------------------------------------------------------------------------------------
   for (int i = 0; i < pr.nmombins; i++) {
      fpCanvas->Clear();
      float invX = 0;
      if (pr.vmom_fph1f_WF[i]->GetMean() > dEdxmax / 2)
         invX = 0.4;
      pavetextBins.Clear();
      pavetextBins.SetX1NDC(0.65 - invX);
      pavetextBins.SetX2NDC(0.98 - invX);
      pavetextBins.Clear();
      pavetextBins.AddText(pr.vmom_fph1f_WF[i]->GetTitle());
      pr.vmom_fph1f_WF[i]->SetTitle("");
      pr.vmom_fph1f_WF[i]->SetMaximum(1.1 *
                                      std::max({pr.vmom_fph1f_WF[i]->GetMaximum(),
                                                pr.vmom_fph1f_XP[i]->GetMaximum()}));
      Graphic_setup(pr.vmom_fph1f_WF[i], 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan, 0.2);
      Graphic_setup(pr.vmom_fph1f_XP[i], 0.5, 1, kMagenta + 1, 2, kMagenta - 2, kMagenta,
                    0.2);
      pr.vmom_fph1f_WF[i]->Draw("HIST");
      pr.vmom_fph1f_XP[i]->Draw("HIST same");
      PrintResolution(pr.vmom_fph1f_XP[i], fpCanvas, 0.65 - invX, 0.53, kMagenta + 2,
                      "XP");
      PrintResolution(pr.vmom_fph1f_WF[i], fpCanvas, 0.65 - invX, 0.2, kCyan + 2, "WF");
      pavetextBins.Draw();
      fpCanvas->SaveAs(fRealpathPDF.c_str());
   }

   // dE/dx VS absolute phi (resolution)
   // -----------------------------------------------------------------------------------------------------------------
   pavetextBins.SetTextSize(0.07);
   for (int i = 0; i < pr.nabsphiresobins; i++) {
      fpCanvas->Clear();
      float invX = 0;
      if (pr.vabsphireso_fph1f_WF[i]->GetMean() > dEdxmax / 2)
         invX = 0.4;
      pavetextBins.Clear();
      pavetextBins.SetX1NDC(0.65 - invX);
      pavetextBins.SetX2NDC(0.98 - invX);
      pavetextBins.Clear();
      pavetextBins.AddText(pr.vabsphireso_fph1f_WF[i]->GetTitle());
      pr.vabsphireso_fph1f_WF[i]->SetTitle("");
      pr.vabsphireso_fph1f_WF[i]->SetMaximum(
         1.1 * std::max({pr.vabsphireso_fph1f_WF[i]->GetMaximum(),
                         pr.vabsphireso_fph1f_XP[i]->GetMaximum()}));
      Graphic_setup(pr.vabsphireso_fph1f_WF[i], 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan,
                    0.2);
      Graphic_setup(pr.vabsphireso_fph1f_XP[i], 0.5, 1, kMagenta + 1, 2, kMagenta - 2,
                    kMagenta, 0.2);
      pr.vabsphireso_fph1f_WF[i]->Draw("HIST");
      pr.vabsphireso_fph1f_XP[i]->Draw("HIST same");
      PrintResolution(pr.vabsphireso_fph1f_XP[i], fpCanvas, 0.65 - invX, 0.53,
                      kMagenta + 2, "XP");
      PrintResolution(pr.vabsphireso_fph1f_WF[i], fpCanvas, 0.65 - invX, 0.2, kCyan + 2,
                      "WF");
      pavetextBins.Draw();
      fpCanvas->SaveAs(fRealpathPDF.c_str());
   }

   // dE/dx VS absolute theta (resolution)
   // -----------------------------------------------------------------------------------------------------------------
   for (int i = 0; i < pr.nabsthetaresobins; i++) {
      fpCanvas->Clear();
      float invX = 0;
      if (pr.vabsthetareso_fph1f_WF[i]->GetMean() > dEdxmax / 2)
         invX = 0.4;
      pavetextBins.Clear();
      pavetextBins.SetX1NDC(0.65 - invX);
      pavetextBins.SetX2NDC(0.98 - invX);
      pavetextBins.Clear();
      pavetextBins.AddText(pr.vabsthetareso_fph1f_WF[i]->GetTitle());
      pr.vabsthetareso_fph1f_WF[i]->SetTitle("");
      pr.vabsthetareso_fph1f_WF[i]->SetMaximum(
         1.1 * std::max({pr.vabsthetareso_fph1f_WF[i]->GetMaximum(),
                         pr.vabsthetareso_fph1f_XP[i]->GetMaximum()}));
      Graphic_setup(pr.vabsthetareso_fph1f_WF[i], 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan,
                    0.2);
      Graphic_setup(pr.vabsthetareso_fph1f_XP[i], 0.5, 1, kMagenta + 1, 2, kMagenta - 2,
                    kMagenta, 0.2);
      pr.vabsthetareso_fph1f_WF[i]->Draw("HIST");
      pr.vabsthetareso_fph1f_XP[i]->Draw("HIST same");
      PrintResolution(pr.vabsthetareso_fph1f_XP[i], fpCanvas, 0.65 - invX, 0.53,
                      kMagenta + 2, "XP");
      PrintResolution(pr.vabsthetareso_fph1f_WF[i], fpCanvas, 0.65 - invX, 0.2, kCyan + 2,
                      "WF");
      pavetextBins.Draw();
      fpCanvas->SaveAs(fRealpathPDF.c_str());
   }

   // Drift velocity plots
   // -----------------------------------------------------------------------------------------------------------------
   fpCanvas->Clear();
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
   // int ncomparisons = (int)v_processes.size();
   int ncomparisons = 2;
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
      // if (pr0->ffileName.find("dog1") != std::string::npos or
      //     pr0->ffileName.find("cosmics") != std::string::npos) {
      //    OutputFile =
      //       foutputComparisonFolder + "/Cosmics_Data_MC" + pr0->fcutslist + ".pdf";
      // } else if (pr0->ffileName.find("beam") != std::string::npos or
      //            pr0->ffileName.find("sandmu") != std::string::npos) {
      //    OutputFile = foutputComparisonFolder + "/Beam_Data_MC" + pr0->fcutslist +
      //    ".pdf";
      // }
      OutputFile = "Output_PDF/BeamCosmics_Data_MC" + pr0->fcutslist + ".pdf";
   }

   // Plot setup
   // -----------------------------------------------------------------------------------------------------------------
   int dEdxmax = pr0->dEdxmax;
   int momrange = pr0->momrange;
   int nmombins = pr0->nmombins;
   fpCanvas->cd();
   fpCanvas->SetGrid(1, 1);
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
   float resomintrklen = 4.5;
   int resomaxtrklen = 16;
   int dedxmintrklen = 300;
   int dedxmaxtrklen = 550;
   int resominphi = 2;
   int resomaxphi = 13;
   int dedxminphi = 350;
   int dedxmaxphi = 500;
   int resomintheta = 2;
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

   // Legends
   TLegend legendXP(0.65, 0.85 - 0.06 * ncomparisons, 0.9, 0.92);
   legendXP.SetTextSize(0.06 - 0.005 * ncomparisons);
   legendXP.SetTextColor(kBlue - 1);
   legendXP.SetFillColorAlpha(kWhite, 0.9);
   TLegend legendXPBC(0.65, 0.85 - 0.06 * (ncomparisons + 0.5), 0.9, 0.92);
   legendXPBC.SetTextSize(0.06 - 0.005 * ncomparisons);
   legendXPBC.SetTextColor(kBlue - 1);
   legendXPBC.SetFillColorAlpha(kWhite, 0.9);
   TLegend legendAll(0.65, 0.85 - 0.06 * (ncomparisons + 1), 0.9, 0.92);
   legendAll.SetTextSize(0.06 - 0.005 * (ncomparisons + 1));
   legendAll.SetTextColor(kBlue - 1);
   legendAll.SetFillColorAlpha(kWhite, 0.9);
   TLegend legendAllBC(0.65, 0.85 - 0.06 * (ncomparisons + 1), 0.9, 0.92);
   legendAllBC.SetTextSize(0.06 - 0.005 * (ncomparisons + 2));
   legendAllBC.SetTextColor(kBlue - 1);
   legendAllBC.SetFillColorAlpha(kWhite, 0.9);

   int resomarkersize = 3;
   int meanmarkersize = 3;
   for (int i = 0; i < ncomparisons + 2; i++) {
      Process *pr_tmp = v_processes[i];
      Graphic_setup(pr_tmp->vabsmomreso_fph1f_XP[7], 0.5, 1, colorsXP[i] + 1, 2,
                    colorsXP[i] + 1, colorsXP[i], 0.4);
      Graphic_setup(pr_tmp->ptge_absmom_reso_XP, resomarkersize, markers[i], colorsXP[i],
                    2, colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_absmom_mean_XP, meanmarkersize, markers[i], colorsXP[i],
                    2, colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_dd_reso_XP, resomarkersize, markers[i], colorsXP[i], 2,
                    colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_dd_mean_XP, meanmarkersize, markers[i], colorsXP[i], 2,
                    colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_trklen_reso_XP, resomarkersize, markers[i], colorsXP[i],
                    2, colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_trklen_mean_XP, meanmarkersize, markers[i], colorsXP[i],
                    2, colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_absphi_reso_XP, resomarkersize, markers[i], colorsXP[i],
                    2, colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_absphi_mean_XP, meanmarkersize, markers[i], colorsXP[i],
                    2, colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_abstheta_reso_XP, resomarkersize, markers[i],
                    colorsXP[i], 2, colorsXP[i]);
      Graphic_setup(pr_tmp->ptge_abstheta_mean_XP, meanmarkersize, markers[i],
                    colorsXP[i], 2, colorsXP[i]);

      Graphic_setup(pr_tmp->ptge_absmom_reso_WF, resomarkersize, markers[i], colorsWF[i],
                    2, colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_absmom_mean_WF, meanmarkersize, markers[i], colorsWF[i],
                    2, colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_dd_reso_WF, resomarkersize, markers[i], colorsWF[i], 2,
                    colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_dd_mean_WF, meanmarkersize, markers[i], colorsWF[i], 2,
                    colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_trklen_reso_WF, resomarkersize, markers[i], colorsWF[i],
                    2, colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_trklen_mean_WF, meanmarkersize, markers[i], colorsWF[i],
                    2, colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_absphi_reso_WF, resomarkersize, markers[i], colorsWF[i],
                    2, colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_absphi_mean_WF, meanmarkersize, markers[i], colorsWF[i],
                    2, colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_abstheta_reso_WF, resomarkersize, markers[i],
                    colorsWF[i], 2, colorsWF[i]);
      Graphic_setup(pr_tmp->ptge_abstheta_mean_WF, meanmarkersize, markers[i],
                    colorsWF[i], 2, colorsWF[i]);
   }

   // Legends fiiling
   std::vector<std::string> legEntries;
   // if (pr0->ffileName.find("dog1") != std::string::npos or
   //     pr0->ffileName.find("cosmics") != std::string::npos) {
   //    legEntries.push_back("Cosmic rays");
   // } else if (pr0->ffileName.find("beam") != std::string::npos or
   //            pr0->ffileName.find("sandmu") != std::string::npos) {
   //    legEntries.push_back("Beam data");
   // }
   legEntries.push_back("Beam data");
   legEntries.push_back("Simulation");
   legEntries.push_back("Cosmic rays");
   legEntries.push_back("Simulation");
   std::vector<TGraph> vdummyXPBC;
   std::vector<TGraph> vdummyWFBC;
   for (int i = 0; i < ncomparisons + 2; i++) {
      TGraph dummyXPBC;
      dummyXPBC.SetMarkerStyle(markers[i]);
      dummyXPBC.SetMarkerSize(7);
      dummyXPBC.SetMarkerColor(colorsXP[i]);
      vdummyXPBC.push_back(dummyXPBC);
      TGraph dummyWFBC;
      dummyWFBC.SetMarkerStyle(markers[i]);
      dummyWFBC.SetMarkerSize(7);
      dummyWFBC.SetMarkerColor(colorsWF[i]);
      vdummyWFBC.push_back(dummyWFBC);
   }
   legendXP.AddEntry(&vdummyXPBC[0], legEntries[0].c_str(), "P");
   legendXP.AddEntry(&vdummyXPBC[1], legEntries[1].c_str(), "P");
   legendXPBC.AddEntry(&vdummyXPBC[0], legEntries[0].c_str(), "P");
   legendXPBC.AddEntry(&vdummyXPBC[2], legEntries[2].c_str(), "P");
   legendXPBC.AddEntry(&vdummyXPBC[1], legEntries[1].c_str(), "P");
   legendAll.AddEntry(&vdummyXPBC[0], ("XP " + legEntries[0]).c_str(), "P");
   legendAll.AddEntry(&vdummyWFBC[0], ("WF " + legEntries[0]).c_str(), "P");
   legendAll.AddEntry(&vdummyXPBC[1], ("XP " + legEntries[1]).c_str(), "P");
   legendAll.AddEntry(&vdummyWFBC[1], ("WF " + legEntries[1]).c_str(), "P");
   legendAllBC.AddEntry(&vdummyXPBC[0], ("XP " + legEntries[0]).c_str(), "P");
   legendAllBC.AddEntry(&vdummyWFBC[0], ("WF " + legEntries[0]).c_str(), "P");
   legendAllBC.AddEntry(&vdummyXPBC[2], ("XP " + legEntries[2]).c_str(), "P");
   legendAllBC.AddEntry(&vdummyWFBC[2], ("WF " + legEntries[2]).c_str(), "P");
   legendAllBC.AddEntry(&vdummyXPBC[1], ("XP " + legEntries[1]).c_str(), "P");
   legendAllBC.AddEntry(&vdummyWFBC[1], ("WF " + legEntries[1]).c_str(), "P");

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
   // Fit on MC
   v_processes[1]->ptge_absmom_mean_XP->Fit(bethefixmuon, "R", "", 100, 2000);
   v_processes[1]->ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetLineColor(kBlack);
   v_processes[1]->ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetLineWidth(2);
   v_processes[1]->ptge_absmom_mean_XP->GetFunction("bethefixmuon")->SetRange(0, 2000);
   // Proton
   TF1 *bethefixproton = (TF1 *)bethefixmuon->Clone("bethefixproton");
   bethefixproton->SetParameter(5, 938.272);
   // Electron
   TF1 *bethefixelectron = (TF1 *)bethefixmuon->Clone("bethefixelectron");
   bethefixelectron->SetParameter(5, 0.511);

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs Track length fit
   const char *trklenresoform = "[0] + [1]*pow(x, [2])";
   TF1 trklenresofit("trklenresofit", trklenresoform, 0, 2000);
   trklenresofit.SetParameter(2, -0.5);
   trklenresofit.SetLineColor(kBlack);
   trklenresofit.SetLineStyle(9);
   pr0->ptge_trklen_reso_XP->Fit(&trklenresofit, "R", "", 200, 1700);
   pr0->ptge_trklen_reso_XP->GetFunction("trklenresofit")->SetLineColor(kBlack);
   pr0->ptge_trklen_reso_XP->GetFunction("trklenresofit")->SetLineWidth(2);
   pr0->ptge_trklen_reso_XP->GetFunction("trklenresofit")->SetRange(0, 2000);

   // -----------------------------------------------------------------------------------------------------------------
   // Global dE/dx plot (real data vs simulation)
   for (int i = 1; i < ncomparisons; i++)
      v_processes[i]->vabsmomreso_fph1f_XP[7]->Scale(
         pr0->vabsmomreso_fph1f_XP[7]->Integral() /
         v_processes[i]->vabsmomreso_fph1f_XP[7]->Integral());
   float maxcount = std::max(pr0->vabsmomreso_fph1f_XP[7]->GetMaximum(),
                             v_processes[1]->vabsmomreso_fph1f_XP[7]->GetMaximum());
   std::vector<std::string> legEntriesDistrib;
   legEntriesDistrib.push_back(legEntries[0]);
   legEntriesDistrib.push_back("Simulation");
   float xwidth = 0.3;
   float xright = 0.91;
   float ymax = 0.9;
   float ygap = 0.05;
   float ywidth = (ymax - 0.2) / ncomparisons - ygap;
   for (int i = ncomparisons - 1; i >= 0; i--) {
      v_processes[i]->vabsmomreso_fph1f_XP[7]->SetTitle("");
      v_processes[i]->vabsmomreso_fph1f_XP[7]->SetAxisRange(0, 1.1 * maxcount, "Y");
      v_processes[i]->vabsmomreso_fph1f_XP[7]->Draw(i == 1 ? "HIST" : "HIST same");
      float ytop = ymax - i * (ywidth + ygap);
      PrintResolution(v_processes[i]->vabsmomreso_fph1f_XP[7], fpCanvas, xright, ytop,
                      xwidth, ywidth, "north east", colors[i] + 1, legEntriesDistrib[i]);
   }
   fpCanvas->SaveAs((OutputFile + "(").c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // Global dE/dx plot (real data only)
   pr0->vabsmomreso_fph1f_XP[7]->SetAxisRange(200, 800, "X");
   pr0->vabsmomreso_fph1f_XP[7]->SetTitle("");
   int maxBin = pr0->vabsmomreso_fph1f_XP[7]->GetMaximumBin();
   float maxValue = pr0->vabsmomreso_fph1f_XP[7]->GetBinContent(maxBin);
   pr0->vabsmomreso_fph1f_XP[7]->SetAxisRange(0, 1.1 * maxValue, "Y");
   pr0->vabsmomreso_fph1f_XP[7]->Draw("HIST");
   PrintResolution(pr0->vabsmomreso_fph1f_XP[7], fpCanvas, xright, ymax, xwidth, ywidth,
                   "north east", colors[0] + 1, "T2K beam data");
   fpCanvas->SaveAs(OutputFile.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // Global dE/dx plot (horizontal and vertical tracks)
   int nbinsphireso = v_processes[2]->nabsphiresobins - 1;
   TH1F h1f_beam(*v_processes[0]->vabsphireso_fph1f_XP[0]);
   TH1F h1f_cosmics(*v_processes[2]->vabsphireso_fph1f_XP[nbinsphireso]);
   h1f_beam.Scale(1.0 / h1f_beam.Integral());
   h1f_cosmics.Scale(1.0 / h1f_cosmics.Integral());
   float maxbeam = h1f_beam.GetBinContent(h1f_beam.GetMaximumBin());
   float maxcosmics = h1f_cosmics.GetBinContent(h1f_cosmics.GetMaximumBin());
   maxValue = std::max(maxbeam, maxcosmics);
   Graphic_setup(&h1f_beam, 0.5, 1, colorsXP[0] + 1, 3, colorsXP[0] + 1, colorsXP[0],
                 0.6);
   Graphic_setup(&h1f_cosmics, 0.5, 1, kMagenta + 3, 3, kMagenta + 3, kMagenta + 2, 0.6);
   std::vector<std::string> legEntriesDistribphi;
   legEntriesDistribphi.push_back("Horizontal beam");
   legEntriesDistribphi.push_back("Vertical cosmics");
   h1f_cosmics.SetTitle("");
   h1f_cosmics.GetYaxis()->SetTitle("Normalized counts");
   h1f_cosmics.SetAxisRange(0, 1.1 * maxValue, "Y");
   h1f_cosmics.SetAxisRange(250, 800, "X");
   h1f_cosmics.Draw("HIST");
   h1f_beam.SetTitle("");
   h1f_beam.Draw("HIST same");
   float ytop = ymax;
   PrintResolution(&h1f_beam, fpCanvas, xright, ytop, xwidth, ywidth, "north east",
                   colors[0] + 1, legEntriesDistribphi[0]);
   ytop = ymax - (ywidth + ygap);
   PrintResolution(&h1f_cosmics, fpCanvas, xright, ytop, xwidth, ywidth, "north east",
                   kMagenta + 2, legEntriesDistribphi[1]);
   fpCanvas->SaveAs(OutputFile.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx vs momentum
   // Resolution
   pr0->ptge_absmom_reso_XP->SetTitle(";Momentum [MeV/c];dE/dx resolution [%]");
   pr0->ptge_absmom_reso_XP->GetXaxis()->SetLimits(0, momrange);
   pr0->ptge_absmom_reso_XP->GetYaxis()->SetRangeUser(resominmom, resomaxmom);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXP.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_absmom_reso_XP->SetTitle(";Momentum [MeV/c];dE/dx resolution [%]");
   pr0->ptge_absmom_reso_XP->GetXaxis()->SetLimits(0, momrange);
   pr0->ptge_absmom_reso_XP->GetYaxis()->SetRangeUser(resominmom, resomaxmom);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absmom_reso_WF->DrawClone("P same");
   }
   legendAll.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // Mean
   pr0->ptge_absmom_mean_XP->SetTitle(";Momentum [MeV/c];dE/dx [ADC counts/cm]");
   pr0->ptge_absmom_mean_XP->GetXaxis()->SetLimits(0, momrange);
   pr0->ptge_absmom_mean_XP->GetYaxis()->SetRangeUser(dedxminmom, dedxmaxmom);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXP.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_absmom_mean_XP->SetTitle(";Momentum [MeV/c];dE/dx [ADC counts/cm]");
   pr0->ptge_absmom_mean_XP->GetXaxis()->SetLimits(0, momrange);
   pr0->ptge_absmom_mean_XP->GetYaxis()->SetRangeUser(dedxminmom, dedxmaxmom);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absmom_mean_WF->DrawClone("P same");
   }
   legendAll.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // Heatmap
   fpCanvas->SetGrid(0, 0);
   gPad->SetRightMargin(0.11);
   pr0->fph2f_XPabsmommean->GetYaxis()->SetRangeUser(200, 1000);
   pr0->fph2f_XPabsmommean->GetYaxis()->SetTitle("dE/dx [ADC counts/cm]");
   pr0->fph2f_XPabsmommean->Draw("colz");
   bethefixmuon->SetLineColor(kCyan);
   bethefixmuon->SetLineStyle(9);
   bethefixmuon->Draw("same");
   bethefixproton->SetLineColor(kRed);
   bethefixproton->SetLineStyle(9);
   bethefixproton->Draw("same");
   bethefixelectron->SetLineColor(kOrange + 1);
   bethefixelectron->SetLineStyle(9);
   bethefixelectron->Draw("same");
   TLegend legbethe(0.57, 0.65, 0.87, 0.9);
   legbethe.SetFillColorAlpha(kWhite, 0.97);
   legbethe.SetTextSize(0.06);
   legbethe.AddEntry(bethefixmuon, "MC muon", "l");
   legbethe.AddEntry(bethefixproton, "MC proton", "l");
   legbethe.AddEntry(bethefixelectron, "MC electron", "l");
   legbethe.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   delete bethefixmuon;
   delete bethefixproton;
   delete bethefixelectron;
   gPad->SetRightMargin(0.035);
   fpCanvas->SetGrid(1, 1);

   // dE/dx vs drift distance
   // -----------------------------------------------------------------------------------------------------------------
   // Resolution
   pr0->ptge_dd_reso_XP->SetTitle(";Drift distance [mm];dE/dx resolution [%]");
   pr0->ptge_dd_reso_XP->GetXaxis()->SetLimits(0, 1000);
   pr0->ptge_dd_reso_XP->GetYaxis()->SetRangeUser(resominX, resomaxX);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_dd_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXP.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_dd_reso_XP->SetTitle(";Drift distance [mm];dE/dx resolution [%]");
   pr0->ptge_dd_reso_XP->GetXaxis()->SetLimits(0, 1000);
   pr0->ptge_dd_reso_XP->GetYaxis()->SetRangeUser(resominX, resomaxX);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_dd_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_dd_reso_WF->DrawClone("P same");
   }
   legendAll.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // Mean
   pr0->ptge_dd_mean_XP->SetTitle(";Drift distance [mm];dE/dx [ADC counts/cm]");
   pr0->ptge_dd_mean_XP->GetXaxis()->SetLimits(0, 1000);
   pr0->ptge_dd_mean_XP->GetYaxis()->SetRangeUser(dedxminX, dedxmaxX);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_dd_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXP.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_dd_mean_XP->SetTitle(";Drift distance [mm];dE/dx [ADC counts/cm]");
   pr0->ptge_dd_mean_XP->GetXaxis()->SetLimits(0, 1000);
   pr0->ptge_dd_mean_XP->GetYaxis()->SetRangeUser(dedxminX, dedxmaxX);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_dd_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_dd_mean_WF->DrawClone("P same");
   }
   legendAll.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx VS Track length
   // Resolution
   pr0->ptge_trklen_reso_XP->SetTitle(";Track length [mm];dE/dx resolution [%]");
   pr0->ptge_trklen_reso_XP->GetXaxis()->SetLimits(0, 1850);
   pr0->ptge_trklen_reso_XP->GetYaxis()->SetRangeUser(resomintrklen, resomaxtrklen);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_trklen_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   // -----------------------------------------------------------------------------------
   TPaveText paveparam(0.65, 0.38, 0.92, 0.63, "NDC");
   paveparam.SetTextColor(kBlue - 5);
   paveparam.SetFillColorAlpha(kWhite, 0.9);
   paveparam.SetTextAlign(12);
   paveparam.SetTextSize(0.045);
   paveparam.SetLineColor(kBlue - 5);
   paveparam.SetLineWidth(1);
   paveparam.SetShadowColor(0);
   std::string parnames[3] = {"#alpha", "#beta", "#gamma"};
   for (int i = 0; i < trklenresofit.GetNpar(); ++i) {
      TString line =
         TString::Format("%s = %.3g #pm %.2g", parnames[i].c_str(),
                         trklenresofit.GetParameter(i), trklenresofit.GetParError(i));
      paveparam.AddText(line);
   }
   TString chi2line = TString::Format(
      "#chi^{2}/ndf = %.2f / %d", trklenresofit.GetChisquare(), trklenresofit.GetNDF());
   paveparam.AddText(chi2line);
   paveparam.Draw();
   legendXP.AddEntry(&trklenresofit, "y = #alpha + #beta x^{#gamma}", "l");
   float previousY1 = legendXP.GetY1NDC();
   float previousY2 = legendXP.GetY2NDC();
   legendXP.SetY1NDC(1.5 * previousY1 - previousY2 / 2);
   legendXP.Draw();
   // -----------------------------------------------------------------------------------
   fpCanvas->SaveAs(OutputFile.c_str());
   legendXP.GetListOfPrimitives()->RemoveLast();
   legendXP.SetY1NDC(previousY1);
   // Add WF
   pr0->ptge_trklen_reso_XP->SetTitle(";Track length [mm];dE/dx resolution [%]");
   pr0->ptge_trklen_reso_XP->GetXaxis()->SetLimits(0, 1850);
   pr0->ptge_trklen_reso_XP->GetYaxis()->SetRangeUser(resomintrklen, resomaxtrklen);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_trklen_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_trklen_reso_WF->DrawClone("P same");
   }
   legendAll.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // Mean
   pr0->ptge_trklen_mean_XP->SetTitle(";Track length [mm];dE/dx [ADC counts/cm]");
   pr0->ptge_trklen_mean_XP->GetXaxis()->SetLimits(0, 1850);
   pr0->ptge_trklen_mean_XP->GetYaxis()->SetRangeUser(dedxmintrklen, dedxmaxtrklen);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_trklen_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXP.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_trklen_mean_XP->SetTitle(";Track length [mm];dE/dx [ADC counts/cm]");
   pr0->ptge_trklen_mean_XP->GetXaxis()->SetLimits(0, 1850);
   pr0->ptge_trklen_mean_XP->GetYaxis()->SetRangeUser(dedxmintrklen, dedxmaxtrklen);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_trklen_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_trklen_mean_WF->DrawClone("P same");
   }
   legendAll.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx VS phi
   // Resolution
   pr0->ptge_absphi_reso_XP->SetTitle(
      ";Track angle #varphi [#circ];dE/dx resolution [%]");
   pr0->ptge_absphi_reso_XP->GetXaxis()->SetLimits(0, 90);
   pr0->ptge_absphi_reso_XP->GetYaxis()->SetRangeUser(resominphi, resomaxphi);
   for (int i = 0; i < ncomparisons + 2; i++) {
      v_processes[i]->ptge_absphi_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXPBC.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_absphi_reso_XP->SetTitle(
      ";Track angle #varphi [#circ];dE/dx resolution [%]");
   pr0->ptge_absphi_reso_XP->GetXaxis()->SetLimits(0, 90);
   pr0->ptge_absphi_reso_XP->GetYaxis()->SetRangeUser(resominphi, resomaxphi);
   for (int i = 0; i < ncomparisons + 2; i++) {
      v_processes[i]->ptge_absphi_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absphi_reso_WF->DrawClone("P same");
   }
   legendAllBC.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // Mean
   pr0->ptge_absphi_mean_XP->SetTitle(
      ";Track angle #varphi [#circ];dE/dx [ADC counts/cm]");
   pr0->ptge_absphi_mean_XP->GetXaxis()->SetLimits(0, 90);
   pr0->ptge_absphi_mean_XP->GetYaxis()->SetRangeUser(dedxminphi, dedxmaxphi);
   for (int i = 0; i < ncomparisons + 2; i++) {
      v_processes[i]->ptge_absphi_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXPBC.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_absphi_mean_XP->SetTitle(
      ";Track angle #varphi [#circ];dE/dx [ADC counts/cm]");
   pr0->ptge_absphi_mean_XP->GetXaxis()->SetLimits(0, 90);
   pr0->ptge_absphi_mean_XP->GetYaxis()->SetRangeUser(dedxminphi, dedxmaxphi);
   for (int i = 0; i < ncomparisons + 2; i++) {
      v_processes[i]->ptge_absphi_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_absphi_mean_WF->DrawClone("P same");
   }
   legendAllBC.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // -----------------------------------------------------------------------------------------------------------------
   // dE/dx VS theta bins
   // Resolution
   pr0->ptge_abstheta_reso_XP->SetTitle(
      ";Track angle #theta [#circ];dE/dx resolution [%]");
   pr0->ptge_abstheta_reso_XP->GetXaxis()->SetLimits(0, 48);
   pr0->ptge_abstheta_reso_XP->GetYaxis()->SetRangeUser(resomintheta, resomaxtheta);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_abstheta_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXP.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_abstheta_reso_XP->SetTitle(
      ";Track angle #theta [#circ];dE/dx resolution [%]");
   pr0->ptge_abstheta_reso_XP->GetXaxis()->SetLimits(0, 48);
   pr0->ptge_abstheta_reso_XP->GetYaxis()->SetRangeUser(resomintheta, resomaxtheta);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_abstheta_reso_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_abstheta_reso_WF->DrawClone("P same");
   }
   legendAll.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // Mean
   pr0->ptge_abstheta_mean_XP->SetTitle(
      ";Track angle #theta [#circ];dE/dx [ADC counts/cm]");
   pr0->ptge_abstheta_mean_XP->GetXaxis()->SetLimits(0, 48);
   pr0->ptge_abstheta_mean_XP->GetYaxis()->SetRangeUser(dedxmintheta, dedxmaxtheta);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_abstheta_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
   }
   legendXP.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   // Add WF
   pr0->ptge_abstheta_mean_XP->SetTitle(
      ";Track angle #theta [#circ];dE/dx [ADC counts/cm]");
   pr0->ptge_abstheta_mean_XP->GetXaxis()->SetLimits(0, 48);
   pr0->ptge_abstheta_mean_XP->GetYaxis()->SetRangeUser(dedxmintheta, dedxmaxtheta);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_abstheta_mean_XP->DrawClone(i == 0 ? "AP" : "P same");
      v_processes[i]->ptge_abstheta_mean_WF->DrawClone("P same");
   }
   legendAll.Draw();
   fpCanvas->SaveAs((OutputFile + ")").c_str());

   for (Process *pr : v_processes)
      delete pr;
}

void Draw::CompareComments(const std::vector<std::string> &v_filepaths)
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
   std::string OutputFile = foutputComparisonFolder + "/Comparison_" + pr0->frun +
                            pr0->ftag + "_" + pr0->fcomment;
   for (int i = 1; i < ncomparisons; i++)
      OutputFile += "_VS_" + v_processes[i]->fcomment;
   OutputFile += pr0->fcutslist + ".pdf";
   std::cout << "Output file: " << OutputFile << std::endl;

   // Plot setup
   // -----------------------------------------------------------------------------------------------------------------
   int dEdxmax = pr0->dEdxmax;
   int momrange = pr0->momrange;
   int nmombins = pr0->nmombins;

   int resominmom = 2;
   int resomaxmom = 12;
   int dedxminmom = 350;
   int dedxmaxmom = 650;
   int resominX = 1;
   int resomaxX = 12;
   int dedxminX = 350;
   int dedxmaxX = 500;
   float resomintrklen = 4.5;
   int resomaxtrklen = 16;
   int dedxmintrklen = 300;
   int dedxmaxtrklen = 550;
   int resominphi = 2;
   int resomaxphi = 13;
   int dedxminphi = 350;
   int dedxmaxphi = 500;
   int resomintheta = 2;
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

   std::vector<int> colorsGain = {kBlue + 1, kGreen + 2, kOrange + 7};
   std::vector<int> markersGain = {20, 22, 23};
   for (int i = 0; i < ncomparisons; i++) {
      Graphic_setup(v_processes[i]->fph1f_WF, 0.5, 1, colorsGain[i], 2, colorsGain[i],
                    colorsGain[i], 0.2);
      Graphic_setup(v_processes[i]->ptge_absmom_reso_WF, 3, markersGain[i], colorsGain[i],
                    2, colorsGain[i]);
      Graphic_setup(v_processes[i]->ptge_absmom_mean_WF, 2, markersGain[i], colorsGain[i],
                    2, colorsGain[i]);
   }

   // Legends filling
   std::vector<std::string> legEntries;
   TLegend legGain(0.55, 0.75, 0.99, 0.99);
   std::vector<TGraph*> v_dummies;
   for (int i = 0; i < ncomparisons; i++) {
      legEntries.push_back(v_processes[i]->flegend);

      TGraph* dummyNoCorr = new TGraph();
      dummyNoCorr->SetMarkerStyle(markersGain[i]);
      dummyNoCorr->SetMarkerSize(7);
      dummyNoCorr->SetMarkerColor(colorsGain[i]);
      v_dummies.push_back(dummyNoCorr);
      legGain.AddEntry(v_dummies[i],
                       v_processes[i]->flegend.c_str(), "p");
   }
   legGain.SetBorderSize(2);
   legGain.SetLineColor(kBlue - 1);
   legGain.SetLineWidth(2);
   legGain.SetTextSize(0.06);
   legGain.SetTextColor(kBlue - 1);

   // -----------------------------------------------------------------------------------
   fpCanvas->cd();
   fpCanvas->SetGrid(1, 1);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gPad->SetBottomMargin(0.15);
   gPad->SetTopMargin(0.06);
   v_processes[0]->fph1f_WF->SetTitleOffset(0.8, "Y");
   // Global dE/dx plot (real data vs simulation)
   float maxcount = 0;
   for (int i = 0; i < ncomparisons; i++)
      maxcount =
         std::max(pr0->fph1f_WF->GetMaximum(), v_processes[i]->fph1f_WF->GetMaximum());
   float xwidth = 0.27;
   float xright = 0.93;
   float ymax = 0.91;
   float ygap = 0.03;
   float ywidth = (ymax - 0.15) / ncomparisons - ygap;

   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->fph1f_WF->SetTitle("");
      v_processes[i]->fph1f_WF->SetAxisRange(0, 1.1 * maxcount, "Y");
      float ytop = ymax - i * (ywidth + ygap);
      v_processes[i]->fph1f_WF->Draw(i == 0 ? "HIST" : "HIST same");
      PrintResolution(v_processes[i]->fph1f_WF, fpCanvas, xright, ytop, xwidth, ywidth,
                      "north east", colorsGain[i], v_processes[i]->flegend, 3, 0.045);
   }
   fpCanvas->SaveAs((OutputFile + "(").c_str());

   // -----------------------------------------------------------------------------------
   // dE/dx resolution as a function of momentum
   fpCanvas->Clear();
   gPad->SetRightMargin(0.05);
   gPad->SetLeftMargin(0.14);
   gPad->SetTopMargin(0.03);
   double minGainReso = 1e4, maxGainReso = 0;
   for (int i = 0; i < ncomparisons; i++) {
      maxGainReso = std::max(
         maxGainReso, v_processes[i]->ptge_absmom_reso_WF->GetHistogram()->GetMaximum());
      minGainReso = std::min(
         minGainReso, v_processes[i]->ptge_absmom_reso_WF->GetHistogram()->GetMinimum());
   }
   v_processes[0]->ptge_absmom_reso_WF->SetTitle(
      ";Momentum [MeV/c];dE/dx resolution [%]");
   v_processes[0]->ptge_absmom_reso_WF->GetXaxis()->SetLimits(0, momrange);
   v_processes[0]->ptge_absmom_reso_WF->GetYaxis()->SetRangeUser(6,
                                                                 17);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_reso_WF->DrawClone(i == 0 ? "AP" : "P same");
   }
   legGain.Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // -----------------------------------------------------------------------------------
   // dE/dx mean value as a function of momentum
   double minGainMean = 1e4, maxGainMean = 0;
   for (int i = 0; i < ncomparisons; i++) {
      maxGainMean = std::max(
         maxGainMean, v_processes[i]->ptge_absmom_mean_WF->GetHistogram()->GetMaximum());
      minGainMean = std::min(
         minGainMean, v_processes[i]->ptge_absmom_mean_WF->GetHistogram()->GetMinimum());
   }
   v_processes[0]->ptge_absmom_mean_WF->SetTitle(
      ";Momentum [MeV/c];dE/dx [ADC counts/cm]");
   v_processes[0]->ptge_absmom_mean_WF->GetXaxis()->SetLimits(0, momrange);
   v_processes[0]->ptge_absmom_mean_WF->GetYaxis()->SetRangeUser(minGainMean,
                                                                 maxGainMean);
   for (int i = 0; i < ncomparisons; i++) {
      v_processes[i]->ptge_absmom_mean_WF->DrawClone(i == 0 ? "AP" : "P same");
   }
   legGain.Draw();
   fpCanvas->SaveAs((OutputFile + ")").c_str());

   for (Process *pr : v_processes)
      delete pr;
   
   for (TGraph *dummy : v_dummies)
      delete dummy;
}
