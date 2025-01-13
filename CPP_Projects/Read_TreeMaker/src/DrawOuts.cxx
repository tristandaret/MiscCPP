#include "DrawOuts.h"
#include "Misc_Functions.h"
#include "TFrame.h"
#include "SetStyle.h"

#include <fstream>

// Constructor
Draw::Draw()
{
   fpStyle = SetMyStyle();
   gROOT->SetStyle(fpStyle->GetName());
   fpCanvas = new TCanvas("fpCanvas", "fpCanvas", 1800, 1350);
   gPad->UseCurrentStyle();
   gPad->SetTopMargin(0.03);
}

// Destructor
Draw::~Draw()
{
   delete fpCanvas;
   delete fpStyle;
   delete fpLegend;

   delete fph1f_WF;
   delete fph1f_XP;
   // delete fptf1_WF;
   // delete fptf1_XP;
   delete fpTGE_mean_WF;
   delete fpTGE_mean_XP;
   delete fpTGE_std_WF;
   delete fpTGE_std_XP;
   delete fpTGE_reso_WF;
   delete fpTGE_reso_XP;
}

void Draw::Run(const std::string &filepath)
{

   fpInputFile = TFile::Open(filepath.c_str());
   TTree *pTree = (TTree *)fpInputFile->Get("outTree");

   std::vector<int> channel2iD = {24, 30, 28, 19, 21, 13, 9,  2,  26, 17, 23, 29, 1,  10, 11, 3,   /*bottom HATPC*/
                                  47, 16, 14, 15, 42, 45, 37, 36, 20, 38, 7,  44, 43, 39, 41, 46}; /*top	HATPC*/

   int dEdxmax = 1300;
   int nentries = pTree->GetEntries();
   int nmombins = 101;
   int momrange = 3000;
   int mombinwidth = 2 * momrange / (nmombins - 1);
   int momindex = 0;

   // Vectors for dE/dx per ERAM module
   std::vector<TH1F *> vmod_ph1f_XP;
   std::vector<TH1F *> vmod_ph1f_WF;
   std::vector<TH2F *> vmod_ph2f_XPtmean;
   std::vector<TH2F *> vmod_ph2f_XPdrift;
   for (int i = 0; i < 32; i++) {
      vmod_ph1f_XP.push_back(new TH1F(Form("ph1f_XP_%d", i),
                                      Form("Energy loss in ERAM %i;dE/dx (ADC counts/cm);Count", channel2iD[i]), 100, 0,
                                      dEdxmax));
      vmod_ph1f_WF.push_back(new TH1F(Form("ph1f_WF_%d", i),
                                      Form("Energy loss in ERAM %i;dE/dx (ADC counts/cm);Count", channel2iD[i]), 100, 0,
                                      dEdxmax));
      vmod_ph2f_XPtmean.push_back(
         new TH2F(Form("ph2f_XPtmean_%d", i),
                  Form("Energy loss (XP) vs mean time in ERAM %i;mean time (ns);dE/dx (ADC counts/cm)", channel2iD[i]),
                  510, 0, 510, 100, 0, dEdxmax));
      vmod_ph2f_XPdrift.push_back(new TH2F(
         Form("ph2f_XPdrift_%d", i),
         Form("Energy loss (XP) vs drift time in ERAM %i;drift time (timebins);dE/dx (ADC counts/cm)", channel2iD[i]),
         510, 0, 510, 100, 0, dEdxmax));
   }

   // Vectors for dE/dx per momentum bin
   std::vector<TH1F *> v_ph1f_WF;
   std::vector<TH1F *> v_ph1f_XP;
   TGraphErrors *ptge_mean_WF = new TGraphErrors();
   TGraphErrors *ptge_mean_XP = new TGraphErrors();
   TGraphErrors *ptge_std_WF = new TGraphErrors();
   TGraphErrors *ptge_std_XP = new TGraphErrors();
   TGraphErrors *ptge_reso_WF = new TGraphErrors();
   TGraphErrors *ptge_reso_XP = new TGraphErrors();
   for (int i = 0; i < nmombins; i++) {
      int mommin = i * mombinwidth - momrange;
      int mommax = (i + 1) * mombinwidth - momrange;
      v_ph1f_WF.push_back(new TH1F(Form("v_ph1f_WF_%d_%d", mommin, mommax),
                                   Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count", mommin, mommax), 100,
                                   0, dEdxmax));
      v_ph1f_XP.push_back(new TH1F(Form("v_ph1f_XP_%d_%d", mommin, mommax),
                                   Form("Energy loss | %d < p < %d; dE/dx (ADC counts/cm); Count", mommin, mommax), 100,
                                   0, dEdxmax));
   }

   // Base
   TH1F *ph1f_WF = new TH1F("ph1f_WF", ";dE/dx(ADC counts/cm);Count", 100, 0, dEdxmax);
   TH1F *ph1f_XP = new TH1F("ph1f_XP", ";dE/dx(ADC counts/cm);Count", 100, 0, dEdxmax);
   // 2D with dE/dx
   TH2F *ph2f_WFXP =
      new TH2F("ph1f_WFXP", ";dE/dx with WF (ADC counts/cm);dE/dx with XP (ADC counts/cm)", 100, 0, 1000, 100, 0, 1000);
   TH2F *ph2f_XPdrift =
      new TH2F("ph2f_XPdrift", ";drift time (timebins);dE/dx with WF (ADC counts/cm)", 510, 0, 510, 100, 0, 1000);
   TH2F *ph2f_lenXP =
      new TH2F("ph1f_lenXP", ";track length (cm);dE/dx with XP (ADC counts/cm)", 171, 0, 170, 100, 0, 1000);
   // Momentum
   TH2F *ph2f_WFmom = new TH2F("ph2f_WFmom", ";momentum (MeV);dE/dx with WF (ADC counts/cm)", nmombins, -momrange,
                               momrange, 100, 0, 1000);
   TH2F *ph2f_XPmom = new TH2F("ph2f_XPmom", ";momentum (MeV);dE/dx with XP (ADC counts/cm)", nmombins, -momrange,
                               momrange, 100, 0, 1000);
   TH1I *ph1i_mom = new TH1I("ph1i_mom", ";momentum (MeV);Count", nmombins, -momrange, momrange);
   TH1I *ph1i_mom_tHAT = new TH1I("ph1i_mom_thAT", ";momentum (MeV);Count", nmombins, -momrange, momrange);
   TH1I *ph1i_mom_bHAT = new TH1I("ph1i_mom_bhAT", ";momentum (MeV);Count", nmombins, -momrange, momrange);
   // Geometry
   TH2F *ph2f_XZ = new TH2F("ph1f_XZ", ";X;Z", 100, -1, 1, 100, -1, 1);
   TH1F *ph1f_dirY = new TH1F("ph1f_dirY", ";Y direction;Count", 100, -1, 1);
   TH1F *ph1f_trklen = new TH1F("ph1f_trklen", ";track length (cm);Count", 171, 0, 170);
   TH1F *ph1f_chi2 = new TH1F("ph1f_chi2", ";#chi^{2};Count", 1000, 0, 1e4);
   TH2F *ph1f_chi2mom =
      new TH2F("ph1f_chi2mom", ";momentum (MeV);#chi^{2}", nmombins, -momrange, momrange, 1000, 0, 1e4);
   // Time in bHAT
   TH1I *ph1i_tminBotCath = new TH1I("ph1i_tminBotCath", ";time bin;Count", 510, 0, 510);
   TH1I *ph1i_tminEP0 = new TH1I("ph1i_tminEP0", "Start time in EP0;time bin;Count", 510, 0, 510);
   TH1I *ph1i_tminEP1 = new TH1I("ph1i_tminEP1", "Start time in EP1;time bin;Count", 510, 0, 510);
   TH1I *ph1i_tmaxBotCath = new TH1I("ph1i_tmaxBotCath", ";time bin;Count", 510, 0, 510);
   TH1I *ph1i_tmaxEP0 = new TH1I("ph1i_tmaxEP0", "End time in EP0;time bin;Count", 510, 0, 510);
   TH1I *ph1i_tmaxEP1 = new TH1I("ph1i_tmaxEP1", "End time in EP1;time bin;Count", 510, 0, 510);
   // Time in tHAT
   TH1I *ph1i_tminTopCath = new TH1I("ph1i_tminTopCath", "Start time in tHATPC;time bin;Count", 510, 0, 510);
   TH1I *ph1i_tminEP2 = new TH1I("ph1i_tminEP2", "Start time in EP2;time bin;Count", 510, 0, 510);
   TH1I *ph1i_tminEP3 = new TH1I("ph1i_tminEP3", "Start time in EP3;time bin;Count", 510, 0, 510);
   TH1I *ph1i_tmaxTopCath = new TH1I("ph1i_tmaxTopCath", "End time in tHATPC;time bin;Count", 510, 0, 510);
   TH1I *ph1i_tmaxEP2 = new TH1I("ph1i_tmaxEP2", "End time in EP2;time bin;Count", 510, 0, 510);
   TH1I *ph1i_tmaxEP3 = new TH1I("ph1i_tmaxEP3", "End time in EP3;time bin;Count", 510, 0, 510);

   Double_t wf;
   Double_t xp;
   Double_t dx;
   Double_t pos[3];
   Double_t dir[3];
   Double_t chi2;
   Double_t mom;
   Double_t APM;
   Double_t start_time;
   Double_t mean_time;
   Double_t end_time;
   Int_t ncl;
   Int_t endplate;
   Int_t eram_channel;
   Int_t eram_ID;
   pTree->SetBranchAddress("dEdx_WF", &wf);
   pTree->SetBranchAddress("dEdx_XP", &xp);
   pTree->SetBranchAddress("track_length", &dx);
   pTree->SetBranchAddress("pos", &pos);
   pTree->SetBranchAddress("dir", &dir);
   pTree->SetBranchAddress("chi2", &chi2);
   pTree->SetBranchAddress("mom", &mom);
   pTree->SetBranchAddress("avg_pad_mult", &APM);
   pTree->SetBranchAddress("starttime", &start_time);
   pTree->SetBranchAddress("mean_time", &mean_time);
   pTree->SetBranchAddress("endtime", &end_time);
   pTree->SetBranchAddress("ncl", &ncl);
   pTree->SetBranchAddress("EndPlate_ID", &endplate);
   pTree->SetBranchAddress("eram_channel", &eram_channel);
   pTree->SetBranchAddress("eram_ID", &eram_ID);

   // Cuts and output file
   int nclmin = 0, nclmax = 200, dirminY = -1, dirmaxY = 1, dxmin = 0, ddEdxmax = 2e3, ncroscut = 0;
   int apmcutlow = 0, apmcuthigh = 20, momcutlow = 0, momcuthigh = 1e6, tcutmin = 0, tcutmax = 1e3;
   int chi2max = 1e6;
   // nclmin = 32;							fileName += ("_" + std::to_string(nclmin)		+ "ncl");
   // apmcutlow = 2, apmcuthigh = 4;		fileName += ("_" + std::to_string(apmcutlow)	+ "apm" +
   // std::to_string(apmcuthigh)); momcutlow = 250, momcuthigh=600;		fileName += ("_" + std::to_string(momcutlow)	+
   // "mom" + std::to_string(momcuthigh));
   // momcutlow = 1e4;						fileName += ("_" + std::to_string(momcutlow)	+ "mom");
   // tcutmin = 0, tcutmax = 30;			fileName += ("_" + std::to_string(tcutmin)		+ "tmin" +
   // std::to_string(tcutmax)); nclmin=50, nclmax = 150;			fileName += ("_" + std::to_string(nclmin)		+ "ncl" +
   // std::to_string(nclmax));
   dxmin = 50, ddEdxmax = 150;
   fileName += ("_" + std::to_string(dxmin) + "dx" + std::to_string(ddEdxmax));
   dirminY = 0.7;
   dirmaxY = 1;
   fileName += ("_" + std::to_string(dirminY) + "dirY" + std::to_string(dirmaxY));
   chi2max = 1000;
   fileName += ("_chi2_" + std::to_string(chi2max));
   MakeMyDir("Output_PDF/" + type + "/" + run + "/" + comment);
   std::string OutputFile = "Output_PDF/" + type + "/" + run + "/" + comment + "/" + fileName + ".pdf";
   std::string OutputFile_Beg = OutputFile + "(";
   std::string OutputFile_End = OutputFile + ")";

   for (int i = 0; i < nentries; i++) {
      pTree->GetEntry(i);
      if (fabs(mom) < 1 || std::isnan(mom))
         continue;

      if (nclmin > ncl or ncl > nclmax)
         continue;
      if (dxmin > dx / 10 or ddEdxmax < dx / 10)
         continue;
      if (fabs(dir[1]) < dirminY or fabs(dir[1]) > dirmaxY)
         continue;
      if (dx < dxmin)
         continue;
      if (APM < apmcutlow or APM > apmcuthigh)
         continue;
      if (fabs(mom) < momcutlow or fabs(mom) > momcuthigh)
         continue;
      if (mean_time < tcutmin or mean_time > tcutmax)
         continue;
      if (chi2 > chi2max)
         continue;

      ph1f_WF->Fill(wf / 1.019);
      ph1f_XP->Fill(xp);
      ph2f_WFXP->Fill(wf / 1.019, xp);
      ph2f_lenXP->Fill(dx / 10, xp);
      ph2f_XPdrift->Fill(mean_time, xp);
      vmod_ph2f_XPdrift[eram_channel]->Fill(mean_time, xp);
      vmod_ph1f_WF[eram_channel]->Fill(wf / 1.019);
      vmod_ph1f_XP[eram_channel]->Fill(xp);
      vmod_ph2f_XPtmean[eram_channel]->Fill(mean_time, xp);

      ph2f_XZ->Fill(dir[0], dir[2]);
      ph1f_dirY->Fill(dir[1]);
      ph1f_trklen->Fill(dx / 10);
      ph1f_chi2->Fill(chi2);
      ph1f_chi2mom->Fill(mom, chi2);

      momindex = (int)std::round(mom / mombinwidth) + nmombins / 2;
      if (fabs(mom) < momrange) {
         v_ph1f_WF[momindex]->Fill(wf / 1.019);
         v_ph1f_XP[momindex]->Fill(xp);
      }
      ph1i_mom->Fill(mom);
      if (pos[1] < 0)
         ph1i_mom_bHAT->Fill(mom);
      if (pos[1] > 0)
         ph1i_mom_tHAT->Fill(mom);
      ph2f_WFmom->Fill(mom, wf / 1.019);
      ph2f_XPmom->Fill(mom, xp);

      if (eram_channel < 16) {
         if (endplate == 8) {
            ph1i_tminEP0->Fill(start_time);
            ph1i_tmaxEP0->Fill(end_time);
         }
         if (endplate == 4) {
            ph1i_tminEP1->Fill(start_time);
            ph1i_tmaxEP1->Fill(end_time);
         }
         if (endplate == 12) {
            ph1i_tminBotCath->Fill(start_time);
            ph1i_tmaxBotCath->Fill(end_time);
         }
      } else {
         if (endplate == 8) {
            ph1i_tminEP2->Fill(start_time);
            ph1i_tmaxEP2->Fill(end_time);
         }
         if (endplate == 4) {
            ph1i_tminEP3->Fill(start_time);
            ph1i_tmaxEP3->Fill(end_time);
         }
         if (endplate == 12) {
            ph1i_tminTopCath->Fill(start_time);
            ph1i_tmaxTopCath->Fill(end_time);
         }
      }
   }

   // TGraph filling
   int ivalid = 0;
   for (int i = 0; i < nmombins; i++) {
      int nentries_here = v_ph1f_WF[i]->GetEntries();
      if (nentries_here < 50)
         continue;
      fptf1_WF = Fit1Gauss(v_ph1f_WF[i]);
      fptf1_XP = Fit1Gauss(v_ph1f_XP[i]);

      float mean_WF = fptf1_WF->GetParameter(1);
      float mean_XP = fptf1_XP->GetParameter(1);
      float dmean_WF = fptf1_WF->GetParError(1);
      float dmean_XP = fptf1_XP->GetParError(1);

      float std_WF = fptf1_WF->GetParameter(2);
      float std_XP = fptf1_XP->GetParameter(2);
      float dstd_WF = fptf1_WF->GetParError(2);
      float dstd_XP = fptf1_XP->GetParError(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WF);
      float dreso_XP = GetResoError(fptf1_XP);

      if (mean_WF == 0 || mean_XP == 0)
         continue;

      ptge_mean_WF->SetPoint(ivalid, i * mombinwidth - momrange, mean_WF);
      ptge_mean_WF->SetPointError(ivalid, mombinwidth / 2, dmean_WF);
      ptge_mean_XP->SetPoint(ivalid, i * mombinwidth - momrange, mean_XP);
      ptge_mean_XP->SetPointError(ivalid, mombinwidth / 2, dmean_XP);

      ptge_std_WF->SetPoint(ivalid, i * mombinwidth - momrange, std_WF);
      ptge_std_WF->SetPointError(ivalid, mombinwidth / 2, dstd_WF);
      ptge_std_XP->SetPoint(ivalid, i * mombinwidth - momrange, std_XP);
      ptge_std_XP->SetPointError(ivalid, mombinwidth / 2, dstd_XP);

      ptge_reso_WF->SetPoint(ivalid, i * mombinwidth - momrange, reso_WF);
      ptge_reso_WF->SetPointError(ivalid, mombinwidth / 2, dreso_WF);
      ptge_reso_XP->SetPoint(ivalid, i * mombinwidth - momrange, reso_XP);
      ptge_reso_XP->SetPointError(ivalid, mombinwidth / 2, dreso_XP);
      ivalid++;
   }

   // 	// Bethe-Bloch curve
   // -------------------------------------------------------------------------------------------------------------------
   // float MeV =						5.9e-3/(224*1493/183);
   // const char *formula =			"[1]/pow(x/sqrt(x*x+[0]*[0]),[4])	*( [2] - pow(x/sqrt(x*x+[0]*[0]),[4]) - log([3]+
   // pow(x/[0], [5])) )/[6]"; TF1 *BBcurve =					new TF1(Form("dEdx_cosmics"), formula, 1, 50e3, "");
   // BBcurve->						SetParameters(105.658, 0.186543e-3, 5.382656, 0.004234, 2.028548, -0.994807, MeV);

   // Plot setup
   // ------------------------------------------------------------------------------------------------------------------------------
   fpCanvas->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   fpLegend = new TLegend(0.6, 0.78, 0.9, 0.91);
   fpLegend->SetTextSize(0.06);
   fpLegend->SetFillStyle(0);
   fpLegend->SetTextColor(kBlue - 1);
   fpLegend->AddEntry(ptge_reso_WF, "Waveforms sum", "p");
   fpLegend->AddEntry(ptge_reso_XP, "Crossed pads", "p");
   gPad->SetTopMargin(0.05);
   float invX = 0;

   // Global dE/dx plot
   // -------------------------------------------------------------------------------------------------------------------
   Graphic_setup(ph1f_WF, 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan, 0.2);
   Graphic_setup(ph1f_XP, 0.5, 1, kMagenta + 2, 2, kMagenta - 2, kMagenta, 0.2);
   ph1f_WF->SetAxisRange(0, 1.1 * std::max({ph1f_WF->GetMaximum(), ph1f_XP->GetMaximum()}), "Y");
   if (ph1f_WF->GetMean() > dEdxmax / 2)
      invX = 0.4;
   ph1f_WF->Draw("HIST");
   ph1f_XP->Draw("HIST sames");
   PrintResolution(ph1f_XP, fpCanvas, 0.65 - invX, 0.58, kMagenta + 2, "XP");
   PrintResolution(ph1f_WF, fpCanvas, 0.65 - invX, 0.25, kCyan + 2, "WF");
   fpCanvas->SaveAs(OutputFile_Beg.c_str());

   // dE/dx ERAM by ERAM
   // ----------------------------------------------------------------------------------------------------------------------
   float maxdEdx = 0;
   for (TH1F *hist : vmod_ph1f_WF)
      if (hist->GetMaximum() > maxdEdx)
         maxdEdx = hist->GetMaximum();
   for (TH1F *hist : vmod_ph1f_XP)
      if (hist->GetMaximum() > maxdEdx)
         maxdEdx = hist->GetMaximum();

   float xMax = 0, yMax = 0, reso = 0, dreso = 0;
   fpCanvas->Clear();
   fpCanvas->Divide(4, 4);
   for (int i = 0; i < 32; i++) {
      vmod_ph1f_WF[i]->SetAxisRange(0, 1.1 * maxdEdx, "Y");
      vmod_ph1f_XP[i]->SetAxisRange(0, 1.1 * maxdEdx, "Y");
      Graphic_setup(vmod_ph1f_WF[i], 0.5, 1, kCyan + 1, 1, kCyan - 2, kCyan, 0.2);
      Graphic_setup(vmod_ph1f_XP[i], 0.5, 1, kMagenta + 2, 1, kMagenta - 2, kMagenta, 0.2);
   }
   for (int i = 0; i < 16; i++) {
      fpCanvas->cd(i + 1);
      vmod_ph1f_WF[i]->Draw();
      vmod_ph1f_XP[i]->Draw("same");
      if (vmod_ph1f_WF[i]->GetEntries() < 100)
         continue;
      xMax = vmod_ph1f_WF[i]->GetXaxis()->GetXmax();
      yMax = vmod_ph1f_WF[i]->GetMaximum();
      vmod_ph1f_WF[i]->GetMean() > dEdxmax / 2 ? invX = 0.4 : invX = 0;
      PrintResolution(vmod_ph1f_WF[i], fpCanvas, 0.65 - invX, 0.58, kCyan + 2, "WF");
      PrintResolution(vmod_ph1f_XP[i], fpCanvas, 0.65 - invX, 0.25, kMagenta + 2, "XP");
   }
   fpCanvas->SaveAs(OutputFile.c_str());

   for (int i = 16; i < 32; i++) {
      fpCanvas->cd(i + 1 - 16);
      vmod_ph1f_WF[i]->Draw();
      vmod_ph1f_XP[i]->Draw("same");

      if (vmod_ph1f_WF[i]->GetEntries() < 100)
         continue;
      xMax = vmod_ph1f_WF[i]->GetXaxis()->GetXmax();
      yMax = vmod_ph1f_WF[i]->GetMaximum();
      vmod_ph1f_WF[i]->GetMean() > dEdxmax / 2 ? invX = 0.4 : invX = 0;
      PrintResolution(vmod_ph1f_WF[i], fpCanvas, 0.65 - invX, 0.58, kCyan + 2, "WF");
      PrintResolution(vmod_ph1f_XP[i], fpCanvas, 0.65 - invX, 0.25, kMagenta + 2, "XP");
   }
   fpCanvas->SaveAs(OutputFile.c_str());

   // dE/dx XP vs WF
   // ---------------------------------------------------------------------------------------------------------------------------
   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   gStyle->SetOptStat(111111);
   gStyle->SetStatX(0.33);
   gStyle->SetStatY(0.95);
   ph2f_WFXP->Draw("colz");
   fpCanvas->SaveAs(OutputFile.c_str());

   // dE/dx vs track length
   // ---------------------------------------------------------------------------------------------------------------------
   fpCanvas->Clear();
   gStyle->SetStatX(0.87);
   ph2f_lenXP->Draw("colz");
   fpCanvas->SaveAs(OutputFile.c_str());

   // dE/dx per time bin
   // -----------------------------------------------------------------------------------------------------------------------
   fpCanvas->Clear();
   gPad->SetRightMargin(0.02);
   gPad->SetRightMargin(0.13);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   ph2f_XPdrift->Draw("colz");
   TGraphErrors *ptge_WFdrift = Convert_TH2_TGE(ph2f_XPdrift);
   TF1 *linearFit = new TF1("linearFit", "pol1", 50, 250);
   ptge_WFdrift->Fit(linearFit, "RQ");
   linearFit->SetLineColor(kRed);
   linearFit->Draw("same");
   fpCanvas->SaveAs(OutputFile.c_str());

   for (int i = 0; i < 32; i++) {
      if (i % 16 == 0) {
         fpCanvas->Clear();
         fpCanvas->Divide(4, 4);
      }
      fpCanvas->cd(i % 16 + 1);
      fpCanvas->cd(i % 16 + 1)->SetRightMargin(0.1);
      vmod_ph2f_XPtmean[i]->Draw("colz");
      if (i % 16 == 15)
         fpCanvas->SaveAs(OutputFile.c_str());
   }

   // dE/dx vs momentum
   // -------------------------------------------------------------------------------------------------------------------
   // Resolution
   fpCanvas->Clear();
   gPad->SetRightMargin(0.03);
   Graphic_setup(ptge_reso_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(ptge_reso_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   ptge_reso_XP->SetTitle(";Momentum (MeV/c);Resolution (%)");
   ptge_reso_XP->GetXaxis()->SetLimits(-momrange, momrange);
   ptge_reso_XP->GetYaxis()->SetRangeUser(5, 30);
   ptge_reso_XP->DrawClone("AP");
   ptge_reso_WF->DrawClone("P same");
   ptge_reso_XP->SetMarkerSize(7);
   ptge_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs((OutputFile + "(").c_str());

   // Resolution's mean
   fpCanvas->Clear();
   Graphic_setup(ptge_mean_WF, 0, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(ptge_mean_XP, 0, 47, kMagenta + 2, 2, kMagenta + 2);
   ptge_mean_XP->SetTitle(";Momentum (MeV/c);Mean (ADC counts/cm)");
   ptge_mean_XP->GetXaxis()->SetLimits(-momrange, momrange);
   ptge_mean_XP->GetYaxis()->SetRangeUser(200, 900);
   ptge_mean_XP->Draw("APL");
   ptge_mean_WF->Draw("PL same");
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // Resolution's standard deviation
   fpCanvas->Clear();
   Graphic_setup(ptge_std_WF, 2, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(ptge_std_XP, 2, 47, kMagenta + 2, 2, kMagenta + 2);
   ptge_std_XP->SetTitle(";Momentum (MeV/c);Standard deviation (ADC counts/cm)");
   ptge_std_XP->GetXaxis()->SetLimits(-momrange, momrange);
   ptge_std_XP->GetYaxis()->SetRangeUser(25, 150);
   ptge_std_XP->Draw("AP");
   ptge_std_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());

   // Momentum plots
   // ------------------------------------------------------------------------------------------------------------------------------
   gPad->SetRightMargin(0.13);
   fpCanvas->Clear();
   ph2f_WFmom->Draw("colz");
   fpCanvas->SaveAs(OutputFile.c_str());

   fpCanvas->Clear();
   ph2f_XPmom->Draw("colz");
   fpCanvas->SaveAs(OutputFile.c_str());

   fpCanvas->Clear();
   gStyle->SetOptStat(0);
   gPad->SetRightMargin(0.03);
   Graphic_setup(ph1i_mom, 0.5, 1, kMagenta + 2, 2, kMagenta - 2, kMagenta, 0.2);
   Graphic_setup(ph1i_mom_tHAT, 0.5, 1, kRed + 2, 2, kRed - 2, kRed, 0.2);
   Graphic_setup(ph1i_mom_bHAT, 0.5, 1, kBlue + 2, 2, kBlue - 2, kBlue, 0.2);
   ph1i_mom->Draw();
   ph1i_mom_bHAT->Draw("same");
   ph1i_mom_tHAT->Draw("same");
   TLegend *plegmom = new TLegend(0.7, 0.7, 0.9, 0.9);
   plegmom->AddEntry(ph1i_mom, "Both HATs", "l");
   plegmom->AddEntry(ph1i_mom_tHAT, "tHAT", "l");
   plegmom->AddEntry(ph1i_mom_bHAT, "bHAT", "l");
   plegmom->SetTextSize(0.06);
   plegmom->SetTextColor(kBlue - 1);
   plegmom->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   delete plegmom;

   // dEdx per momentum bin
   // --------------------------------------------------------------------------------------------------------------------

   gPad->SetTopMargin(0.08);
   gPad->SetRightMargin(0.02);
   for (int i = 0; i < nmombins; i++) {
      fpCanvas->Clear();
      v_ph1f_WF[i]->SetAxisRange(0, 1.1 * std::max({v_ph1f_WF[i]->GetMaximum(), v_ph1f_XP[i]->GetMaximum()}), "Y");
      float invX = 0;
      if (v_ph1f_WF[i]->GetMean() > dEdxmax / 2)
         invX = 0.4;
      Graphic_setup(v_ph1f_WF[i], 0.5, 1, kCyan + 1, 2, kCyan - 2, kCyan, 0.2);
      Graphic_setup(v_ph1f_XP[i], 0.5, 1, kMagenta + 1, 2, kMagenta - 2, kMagenta, 0.2);
      v_ph1f_WF[i]->Draw("HIST");
      v_ph1f_XP[i]->Draw("HIST same");
      PrintResolution(v_ph1f_XP[i], fpCanvas, 0.65 - invX, 0.58, kMagenta + 2, "XP");
      PrintResolution(v_ph1f_WF[i], fpCanvas, 0.65 - invX, 0.25, kCyan + 2, "WF");
      fpCanvas->SaveAs(OutputFile.c_str());
   }

   // Geometry plots
   // ----------------------------------------------------------------------------------------------------------------------------
   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   ph2f_XZ->Draw("colz");
   fpCanvas->SaveAs(OutputFile.c_str());

   fpCanvas->Clear();
   gPad->SetRightMargin(0.02);
   ph1f_dirY->Draw("HIST");
   fpCanvas->SaveAs(OutputFile.c_str());

   fpCanvas->Clear();
   ph1f_trklen->Draw("HIST");
   fpCanvas->SaveAs(OutputFile.c_str());

   fpCanvas->Clear();
   gPad->SetRightMargin(0.07);
   ph1f_chi2->Draw("HIST");
   fpCanvas->SaveAs(OutputFile.c_str());

   fpCanvas->Clear();
   gPad->SetRightMargin(0.13);
   ph1f_chi2mom->Draw("colz");
   fpCanvas->SaveAs(OutputFile.c_str());

   // Drift velocity plots
   // ----------------------------------------------------------------------------------------------------------------------
   fpCanvas->Clear();
   int bins = 100;
   int tmin = -999;
   int max = -999;
   for (int i = 0; i < bins; i++) {
      if (ph1i_tminEP0->GetBinContent(i) > max) {
         tmin = i;
         max = ph1i_tminEP0->GetBinContent(i);
      }
   }

   int tmax = -999;
   for (int i = 300; i < 300 + bins; i++) {
      if (ph1i_tmaxBotCath->GetBinContent(i) > max) {
         tmax = i;
         max = ph1i_tmaxBotCath->GetBinContent(i);
      }
   }
   std::cout << std::fixed << std::setprecision(3);
   std::cout << "tmin, tmax: " << tmin << " " << tmax << std::endl;

   // ph1i_tminEP0->				Fit("gaus", "RQ", "", tmin-2, tmin+2);
   // ph1i_tminEP2->				Fit("gaus", "RQ", "", tmin-2, tmin+2);
   // ph1i_tmaxBotCath->			Fit("gaus", "RQ", "", tmax-4, tmax+4);
   // ph1i_tmaxTopCath->			Fit("gaus", "RQ", "", tmax-4, tmax+4);
   // TF1 *ptf1_tminEP0 =			ph1i_tminEP0->GetFunction("gaus");
   // TF1 *ptf1_tminEP2 =			ph1i_tminEP2->GetFunction("gaus");
   // TF1 *ptf1_tmaxBotCath =		ph1i_tmaxBotCath->GetFunction("gaus");
   // TF1 *ptf1_tmaxTopCath =		ph1i_tmaxTopCath->GetFunction("gaus");

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

   ph1i_tminBotCath->SetAxisRange(
      0,
      1.1 * std::max({ph1i_tminBotCath->GetMaximum(), ph1i_tminEP0->GetMaximum(), ph1i_tminEP1->GetMaximum(),
                      ph1i_tminTopCath->GetMaximum(), ph1i_tminEP2->GetMaximum(), ph1i_tminEP3->GetMaximum()}),
      "Y");
   ph1i_tminBotCath->SetLineWidth(2);
   ph1i_tminBotCath->SetLineColor(kGreen + 2);
   ph1i_tminEP0->SetLineWidth(2);
   ph1i_tminEP0->SetLineColor(kBlue);
   ph1i_tminEP1->SetLineWidth(2);
   ph1i_tminEP1->SetLineColor(kRed);
   ph1i_tminTopCath->SetLineWidth(2);
   ph1i_tminTopCath->SetLineColor(kGreen - 6);
   ph1i_tminEP2->SetLineWidth(2);
   ph1i_tminEP2->SetLineColor(kBlue - 6);
   ph1i_tminEP3->SetLineWidth(2);
   ph1i_tminEP3->SetLineColor(kRed - 6);
   ph1i_tminBotCath->Draw();
   ph1i_tminEP0->Draw("same");
   ph1i_tminEP1->Draw("same");
   ph1i_tminTopCath->Draw("same");
   ph1i_tminEP2->Draw("same");
   ph1i_tminEP3->Draw("same");
   // ptf1_tminEP0->				SetLineColor(kOrange+2);
   // ptf1_tminEP2->				SetLineColor(kViolet+2);
   // ptf1_tminEP0->				Draw("same");
   // ptf1_tminEP2->				Draw("same");
   TLegend *legtmin = new TLegend(0.6, 0.5, 0.9, 0.9);
   legtmin->AddEntry(ph1i_tminBotCath, "T_{min} (bottom cathode) ", "l");
   legtmin->AddEntry(ph1i_tminEP0, "T_{min} EP0 ", "l");
   legtmin->AddEntry(ph1i_tminEP1, "T_{min} EP1 ", "l");
   legtmin->AddEntry(ph1i_tminTopCath, "T_{min} (top cathode) ", "l");
   legtmin->AddEntry(ph1i_tminEP2, "T_{min} EP2 ", "l");
   legtmin->AddEntry(ph1i_tminEP3, "T_{min} EP3 ", "l");
   legtmin->DrawClone();
   fpCanvas->Update();
   delete legtmin;
   fpCanvas->SaveAs(OutputFile.c_str());

   fpCanvas->Clear();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   ph1i_tmaxBotCath->SetAxisRange(
      0,
      1.1 * std::max({ph1i_tmaxBotCath->GetMaximum(), ph1i_tmaxEP0->GetMaximum(), ph1i_tmaxEP1->GetMaximum(),
                      ph1i_tmaxTopCath->GetMaximum(), ph1i_tmaxEP2->GetMaximum(), ph1i_tmaxEP3->GetMaximum()}),
      "Y");
   ph1i_tmaxBotCath->SetLineWidth(2);
   ph1i_tmaxBotCath->SetLineColor(kGreen + 2);
   ph1i_tmaxEP0->SetLineWidth(2);
   ph1i_tmaxEP0->SetLineColor(kBlue);
   ph1i_tmaxEP1->SetLineWidth(2);
   ph1i_tmaxEP1->SetLineColor(kRed);
   ph1i_tmaxTopCath->SetLineWidth(2);
   ph1i_tmaxTopCath->SetLineColor(kGreen - 6);
   ph1i_tmaxEP2->SetLineWidth(2);
   ph1i_tmaxEP2->SetLineColor(kBlue - 6);
   ph1i_tmaxEP3->SetLineWidth(2);
   ph1i_tmaxEP3->SetLineColor(kRed - 6);
   ph1i_tmaxBotCath->Draw();
   ph1i_tmaxEP0->Draw("same");
   ph1i_tmaxEP1->Draw("same");
   ph1i_tmaxTopCath->Draw("same");
   ph1i_tmaxEP2->Draw("same");
   ph1i_tmaxEP3->Draw("same");
   // ptf1_tmaxBotCath->			SetLineColor(kOrange+2);
   // ptf1_tmaxTopCath->			SetLineColor(kViolet+2);
   // ptf1_tmaxTopCath->			Draw("same");
   // ptf1_tmaxBotCath->			Draw("same");
   TLegend *legtmax = new TLegend(0.15, 0.5, 0.4, 0.9);
   legtmax->AddEntry(ph1i_tmaxBotCath, "T_{max} (bottom cathode) ", "l");
   legtmax->AddEntry(ph1i_tmaxEP0, "T_{max} EP0 ", "l");
   legtmax->AddEntry(ph1i_tmaxEP1, "T_{max} EP1 ", "l");
   legtmax->AddEntry(ph1i_tmaxTopCath, "T_{max} (top cathode) ", "l");
   legtmax->AddEntry(ph1i_tmaxEP2, "T_{max} EP2 ", "l");
   legtmax->AddEntry(ph1i_tmaxEP3, "T_{max} EP3 ", "l");
   legtmax->DrawClone();
   fpCanvas->Update();
   delete legtmax;
   fpCanvas->SaveAs(OutputFile_End.c_str());
}

void Draw::Scan(const std::vector<std::string> &v_filepath)
{

   if (run == "phi")
      v_fvalues = {0, 10, 20, 30, 40, 45, 50, 60, 70, 80, 90};
   if (run == "drift")
      v_fvalues = {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 97};

   fpTGE_mean_WF = new TGraphErrors();
   fpTGE_mean_XP = new TGraphErrors();
   fpTGE_std_WF = new TGraphErrors();
   fpTGE_std_XP = new TGraphErrors();
   fpTGE_reso_WF = new TGraphErrors();
   fpTGE_reso_XP = new TGraphErrors();

   for (int i = 0; i < (int)v_filepath.size(); i++) {
      fpInputFile = TFile::Open(v_filepath[i].c_str());
      fpInputTree = (TTree *)fpInputFile->Get("outTree");
      int nentries = fpInputTree->GetEntries();

      fph1f_WF = new TH1F("fph1f_WF", "Energy loss with WF; dE/dx (ADC counts/cm); Count", 100, 0, 1300);
      fph1f_XP = new TH1F("fph1f_XP", "Energy loss with XP; dE/dx (ADC counts/cm); Count", 100, 0, 1300);

      Double_t wf;
      Double_t xp;
      fpInputTree->SetBranchAddress("dEdx_WF", &wf);
      fpInputTree->SetBranchAddress("dEdx_XP", &xp);

      for (int j = 0; j < nentries; j++) {
         fpInputTree->GetEntry(j);
         fph1f_WF->Fill(wf / 1.019);
         fph1f_XP->Fill(xp);
      }

      fptf1_WF = Fit1Gauss(fph1f_WF);
      fptf1_XP = Fit1Gauss(fph1f_XP);

      float mean_WF = fptf1_WF->GetParameter(1);
      float mean_XP = fptf1_XP->GetParameter(1);
      float dmean_WF = fptf1_WF->GetParError(1);
      float dmean_XP = fptf1_XP->GetParError(1);

      float std_WF = fptf1_WF->GetParameter(2);
      float std_XP = fptf1_XP->GetParameter(2);
      float dstd_WF = fptf1_WF->GetParError(2);
      float dstd_XP = fptf1_XP->GetParError(2);

      float reso_WF = std_WF / mean_WF * 100;
      float reso_XP = std_XP / mean_XP * 100;
      float dreso_WF = GetResoError(fptf1_WF);
      float dreso_XP = GetResoError(fptf1_XP);

      fpTGE_mean_WF->SetPoint(i, v_fvalues[i], mean_WF);
      fpTGE_mean_WF->SetPointError(i, 0, dmean_WF);
      fpTGE_mean_XP->SetPoint(i, v_fvalues[i], mean_XP);
      fpTGE_mean_XP->SetPointError(i, 0, dmean_XP);

      fpTGE_std_WF->SetPoint(i, v_fvalues[i], std_WF);
      fpTGE_std_WF->SetPointError(i, 0, dstd_WF);
      fpTGE_std_XP->SetPoint(i, v_fvalues[i], std_XP);
      fpTGE_std_XP->SetPointError(i, 0, dstd_XP);

      fpTGE_reso_WF->SetPoint(i, v_fvalues[i], reso_WF);
      fpTGE_reso_WF->SetPointError(i, 0, dreso_WF);
      fpTGE_reso_XP->SetPoint(i, v_fvalues[i], reso_XP);
      fpTGE_reso_XP->SetPointError(i, 0, dreso_XP);
   }

   // Drawing
   if (comment == "")
      comment = "reference";
   std::string OutputFile = "Output_PDF/" + type + "/" + run + "/" + comment + "/" + run + "Scan" + comment + ".pdf";
   fpCanvas->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   std::string xlabel;
   if (run == "phi")
      xlabel = "Track angle #varphi (#circ)";
   if (run == "drift")
      xlabel = "Drift distance (cm)";
   // Legend
   fpLegend = new TLegend(0.6, 0.78, 0.9, 0.91);
   fpLegend->SetTextSize(0.06);
   fpLegend->SetFillStyle(0);
   fpLegend->SetTextColor(kBlue - 1);
   fpLegend->AddEntry(fpTGE_reso_WF, "Waveforms sum", "p");
   fpLegend->AddEntry(fpTGE_reso_XP, "Crossed pads", "p");

   int xmin = -3;
   int dEdxmax;
   if (run == "phi")
      dEdxmax = 93;
   if (run == "drift")
      dEdxmax = 100;
   int resomin = 6;
   int resomax = 9;
   int meanmin = 300;
   int meanmax = 450;
   int stdmin = 15;
   int stdmax = 40;
   int meanpullmin = -1;
   int meanpullmax = 1;
   int stdpullmin = 0;
   int stdpullmax = 1;

   // Resolution
   Graphic_setup(fpTGE_reso_WF, 5, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(fpTGE_reso_XP, 4, 47, kMagenta + 2, 2, kMagenta + 2);
   fpTGE_reso_XP->SetTitle((";" + xlabel + ";Resolution (%)").c_str());
   fpTGE_reso_XP->GetXaxis()->SetLimits(xmin, dEdxmax);
   fpTGE_reso_XP->GetYaxis()->SetRangeUser(resomin, resomax);
   fpTGE_reso_XP->DrawClone("AP");
   fpTGE_reso_WF->DrawClone("P same");
   fpTGE_reso_XP->SetMarkerSize(7);
   fpTGE_reso_WF->SetMarkerSize(7);
   fpLegend->Draw();
   fpCanvas->SaveAs((OutputFile + "(").c_str());
   fpCanvas->Clear();

   // Resolution's mean
   Graphic_setup(fpTGE_mean_WF, 5, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(fpTGE_mean_XP, 4, 47, kMagenta + 2, 2, kMagenta + 2);
   fpTGE_mean_XP->SetTitle((";" + xlabel + ";Mean (ADC counts/cm)").c_str());
   fpTGE_mean_XP->GetXaxis()->SetLimits(xmin, dEdxmax);
   fpTGE_mean_XP->GetYaxis()->SetRangeUser(meanmin, meanmax);
   fpTGE_mean_XP->Draw("AP");
   fpTGE_mean_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // Resolution's standard deviation
   Graphic_setup(fpTGE_std_WF, 5, 33, kCyan + 2, 2, kCyan + 2);
   Graphic_setup(fpTGE_std_XP, 4, 47, kMagenta + 2, 2, kMagenta + 2);
   fpTGE_std_XP->SetTitle((";" + xlabel + ";Standard deviation (ADC counts/cm)").c_str());
   fpTGE_std_XP->GetXaxis()->SetLimits(xmin, dEdxmax);
   fpTGE_std_XP->GetYaxis()->SetRangeUser(stdmin, stdmax);
   fpTGE_std_XP->Draw("AP");
   fpTGE_std_WF->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs((OutputFile + ")").c_str());
}

void Draw::CompareScans(const std::vector<std::string> &v_tags, const std::string &method)
{

   if (run == "phi")
      v_fvalues = {0, 10, 20, 30, 40, 45, 50, 60, 70, 80, 90};
   else if (run == "drift")
      v_fvalues = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 98};
   else
      v_fvalues = {1};
   int ncomparisons = (int)v_comment.size();

   std::vector<TGraphErrors *> v_ptge_reso;
   std::vector<TGraphErrors *> v_ptge_mean;
   std::vector<TGraphErrors *> v_ptge_std;
   std::vector<std::vector<TH1F *>> array_ph1f_dedx;
   TF1 *ptf1_dedx;

   for (int j = 0; j < ncomparisons; j++) {
      v_ptge_reso.push_back(new TGraphErrors());
      v_ptge_mean.push_back(new TGraphErrors());
      v_ptge_std.push_back(new TGraphErrors());
      array_ph1f_dedx.push_back(std::vector<TH1F *>());
      for (int i = 0; i < (int)v_tags.size(); i++) {
         fpInputFile = TFile::Open(
            ("ROOT_files/" + type + "/" + run + "/" + v_comment[j] + "/" + type + v_tags[i] + v_comment[j] + ".root")
               .c_str());
         fpInputTree = (TTree *)fpInputFile->Get("outTree");
         int nentries = fpInputTree->GetEntries();

         array_ph1f_dedx[j].push_back(
            new TH1F(Form("array_ph1f_dedx_%s", v_tags[i].c_str()), ";dE/dx [ADC counts/cm];Count", 100, 0, 1300));

         Double_t dedx;
         Double_t dx;
         Double_t dir[3];
         Double_t mom;
         fpInputTree->SetBranchAddress(Form("dEdx_%s", method.c_str()), &dedx);
         fpInputTree->SetBranchAddress("track_length", &dx);
         fpInputTree->SetBranchAddress("dir", &dir);
         fpInputTree->SetBranchAddress("mom", &mom);

         for (int k = 0; k < nentries; k++) {
            fpInputTree->GetEntry(k);
            if (dx < 500)
               continue;
            if (fabs(dir[1] < 0.7))
               continue;
            if (mom < 250 or mom > 1000)
               continue;
            method == "XP" ? array_ph1f_dedx[j][i]->Fill(dedx) : array_ph1f_dedx[j][i]->Fill(dedx / 1.019);
         }

         ptf1_dedx = Fit1Gauss(array_ph1f_dedx[j][i]);

         float mean = ptf1_dedx->GetParameter(1);
         float dmean = ptf1_dedx->GetParError(1);

         float std = ptf1_dedx->GetParameter(2);
         float dstd = ptf1_dedx->GetParError(2);

         float reso = std / mean * 100;
         float dreso = GetResoError(ptf1_dedx);

         v_ptge_mean[j]->SetPoint(i, v_fvalues[i] - 1 + j, mean);
         v_ptge_mean[j]->SetPointError(i, 0, dmean);

         v_ptge_std[j]->SetPoint(i, v_fvalues[i] - 1 + j, std);
         v_ptge_std[j]->SetPointError(i, 0, dstd);

         v_ptge_reso[j]->SetPoint(i, v_fvalues[i] - 1 + j, reso);
         v_ptge_reso[j]->SetPointError(i, 0, dreso);
      }
      if (v_comment[j] == "")
         v_comment[j] = "_reference";
   }
   // Drawing
   std::string OutputFile = "Output_PDF/" + type + "/" + run + "/" + run + "Scan" + v_comment[0];
   for (int i = 1; i < ncomparisons; i++)
      OutputFile += "_VS" + v_comment[i];
   OutputFile += "_" + method + ".pdf";
   fpCanvas->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   std::string xlabel;
   if (run == "phi")
      xlabel = "Track angle #varphi (#circ)";
   if (run == "drift")
      xlabel = "Drift distance (cm)";
   // Legend
   fpLegend = new TLegend(0.55 + 0.03 * ncomparisons, 0.85 - 0.03 * ncomparisons, 0.9, 0.95);
   fpLegend->SetTextSize(0.06 - 0.005 * ncomparisons);
   fpLegend->SetFillStyle(0);
   fpLegend->SetTextColor(kBlue - 1);
   for (int i = 0; i < ncomparisons; i++)
      fpLegend->AddEntry(v_ptge_reso[i], (method + v_comment[i]).c_str(), "p");

   int xmin = -3;
   int dEdxmax;
   if (run == "phi")
      dEdxmax = 93;
   if (run == "drift")
      dEdxmax = 100;
   int resomin = 6;
   int resomax = 9;
   int meanmin = 300;
   int meanmax = 450;
   int stdmin = 15;
   int stdmax = 40;

   // Resolution
   for (int i = 0; i < ncomparisons; i++)
      Graphic_setup(v_ptge_reso[i], 4, markers[i], colors[i], 2, colors[i]);
   v_ptge_reso[0]->SetTitle((";" + xlabel + ";Resolution (%)").c_str());
   v_ptge_reso[0]->GetXaxis()->SetLimits(xmin, dEdxmax);
   v_ptge_reso[0]->GetYaxis()->SetRangeUser(resomin, resomax);
   v_ptge_reso[0]->DrawClone("AP");
   v_ptge_reso[0]->SetMarkerSize(7);
   for (int i = 1; i < ncomparisons; i++) {
      v_ptge_reso[i]->DrawClone("P same");
      v_ptge_reso[i]->SetMarkerSize(7);
   }
   fpLegend->Draw();
   fpCanvas->SaveAs((OutputFile + "(").c_str());
   fpCanvas->Clear();

   // Resolution's mean
   for (int i = 0; i < ncomparisons; i++)
      Graphic_setup(v_ptge_mean[i], 4, markers[i], colors[i], 2, colors[i]);
   v_ptge_mean[0]->SetTitle((";" + xlabel + ";Mean (ADC counts/cm)").c_str());
   v_ptge_mean[0]->GetXaxis()->SetLimits(xmin, dEdxmax);
   v_ptge_mean[0]->GetYaxis()->SetRangeUser(meanmin, meanmax);
   v_ptge_mean[0]->Draw("AP");
   for (int i = 1; i < ncomparisons; i++)
      v_ptge_mean[i]->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   // Resolution's standard deviation
   for (int i = 0; i < ncomparisons; i++)
      Graphic_setup(v_ptge_std[i], 4, markers[i], colors[i], 2, colors[i]);
   v_ptge_std[0]->SetTitle((";" + xlabel + ";Standard deviation (ADC counts/cm)").c_str());
   v_ptge_std[0]->GetXaxis()->SetLimits(xmin, dEdxmax);
   v_ptge_std[0]->GetYaxis()->SetRangeUser(stdmin, stdmax);
   v_ptge_std[0]->Draw("AP");
   for (int i = 1; i < ncomparisons; i++)
      v_ptge_std[i]->Draw("P same");
   fpLegend->Draw();
   fpCanvas->SaveAs(OutputFile.c_str());
   fpCanvas->Clear();

   gPad->SetTopMargin(0.05);
   for (int i = 0; i < (int)v_tags.size(); i++) {
      for (int j = 0; j < ncomparisons; j++) {
         array_ph1f_dedx[j][i]->GetFunction("gausn")->SetLineWidth(0);
         Graphic_setup(array_ph1f_dedx[j][i], 0.5, 1, colors2[j] + 1, 2, colors2[j] - 2, colors2[j], 0.2);
         j == 0 ? array_ph1f_dedx[j][i]->Draw() : array_ph1f_dedx[j][i]->Draw("same");
         std::string name;
         j == 0 ? name = "Gain correction: ON" : name = "Gain correction: OFF";
         PrintResolution(array_ph1f_dedx[j][i], fpCanvas, 0.65, 0.6 - j * 0.35, colors2[j], name.c_str());
      }
      fpCanvas->SaveAs((OutputFile + ")").c_str());
   }
}