#include "Includes.h"
#include "Util.h"
#include "SetStyle.h"

int main()
{

   std::string filenameold = "root_files/TreeMaker_MC_mu-_600MeV_x90_y105_z-248_phi-90_theta0_N100_ld.root";
   TFile fileold(filenameold.c_str(), "READ");
   std::string filenamenew = "root_files/TreeMaker_MC_mu-_600MeV_x90_y105_z-248_phi-90_theta0_N100_new.root";
   TFile filenew(filenamenew.c_str(), "READ");

   std::vector<TH1F *> v_tdiffbottom;
   std::vector<TH1F *> v_tdifftop;
   std::vector<TH1F *> v_qdiffbottom;
   std::vector<TH1F *> v_qdifftop;

   v_tdiffbottom.push_back((TH1F *)fileold.Get("fT2mT1bHAT"));
   v_tdiffbottom.push_back((TH1F *)filenew.Get("fT2mT1bHAT"));
   v_tdifftop.push_back((TH1F *)fileold.Get("fT2mT1tHAT"));
   v_tdifftop.push_back((TH1F *)filenew.Get("fT2mT1tHAT"));
   v_qdiffbottom.push_back((TH1F *)fileold.Get("fQ2overQ1bHAT"));
   v_qdiffbottom.push_back((TH1F *)filenew.Get("fQ2overQ1bHAT"));
   v_qdifftop.push_back((TH1F *)fileold.Get("fQ2overQ1tHAT"));
   v_qdifftop.push_back((TH1F *)filenew.Get("fQ2overQ1tHAT"));

   // Draw histograms
   TStyle *fpStyle = SetMyStyle();
   gROOT->SetStyle(fpStyle->GetName());
   TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1350);
   gPad->UseCurrentStyle();
   gStyle->SetPadTopMargin(0.1);
   c1->Divide(2, 2);

   c1->cd(1);
   v_tdifftop[0]->SetTitle("T_{2} - T_{1} (tHAT)");
   v_tdifftop[0]->SetName("Default");
   Graphic_setup(v_tdifftop[0], 0, 0, kBlue, 1, kBlue, kBlue, 0.2);
   v_tdifftop[0]->SetMaximum(std::max(v_tdifftop[0]->GetMaximum(), v_tdifftop[1]->GetMaximum()) * 1.1);
   v_tdifftop[0]->Draw("hist");
   gPad->Update();
   TPaveStats *pStatTdifftHATold = (TPaveStats *)v_tdifftop[0]->FindObject("stats");
   pStatTdifftHATold->SetTextColor(kBlue);
   pStatTdifftHATold->SetY1NDC(0.7);
   pStatTdifftHATold->SetY2NDC(0.9);
   v_tdifftop[1]->SetName("New RC values");
   Graphic_setup(v_tdifftop[1], 0, 0, kRed, 1, kRed, kRed, 0.2);
   v_tdifftop[1]->Draw("hist sames");
   gPad->Update();
   TPaveStats *pStatTdifftHATnew = (TPaveStats *)v_tdifftop[1]->FindObject("stats");
   pStatTdifftHATnew->SetTextColor(kRed);
   pStatTdifftHATnew->SetY1NDC(0.45);
   pStatTdifftHATnew->SetY2NDC(0.65);

   c1->cd(2);
   v_qdifftop[0]->SetTitle("Q_{2}/Q_{1} (tHAT)");
   v_qdifftop[0]->SetName("Default");
   Graphic_setup(v_qdifftop[0], 0, 0, kBlue, 1, kBlue, kBlue, 0.2);
   v_qdifftop[0]->SetMaximum(std::max(v_qdifftop[0]->GetMaximum(), v_qdifftop[1]->GetMaximum()) * 1.1);
   v_qdifftop[0]->Draw("hist");
   gPad->Update();
   TPaveStats *pStatQdifftHATold = (TPaveStats *)v_qdifftop[0]->FindObject("stats");
   pStatQdifftHATold->SetTextColor(kBlue);
   pStatQdifftHATold->SetY1NDC(0.7);
   pStatQdifftHATold->SetY2NDC(0.9);
   v_qdifftop[1]->SetName("New RC values");
   Graphic_setup(v_qdifftop[1], 0, 0, kRed, 1, kRed, kRed, 0.2);
   v_qdifftop[1]->Draw("hist sames");
   gPad->Update();
   TPaveStats *pStatQdifftHATnew = (TPaveStats *)v_qdifftop[1]->FindObject("stats");
   pStatQdifftHATnew->SetTextColor(kRed);
   pStatQdifftHATnew->SetY1NDC(0.45);
   pStatQdifftHATnew->SetY2NDC(0.65);

   c1->cd(3);
   v_tdiffbottom[0]->SetTitle("T_{2} - T_{1} (bHAT)");
   v_tdiffbottom[0]->SetName("Default");
   Graphic_setup(v_tdiffbottom[0], 0, 0, kBlue, 1, kBlue, kBlue, 0.2);
   v_tdiffbottom[0]->SetMaximum(std::max(v_tdiffbottom[0]->GetMaximum(), v_tdiffbottom[1]->GetMaximum()) * 1.1);
   v_tdiffbottom[0]->Draw("hist");
   gPad->Update();
   TPaveStats *pStatTdiffbHATold = (TPaveStats *)v_tdiffbottom[0]->FindObject("stats");
   pStatTdiffbHATold->SetTextColor(kBlue);
   pStatTdiffbHATold->SetY1NDC(0.7);
   pStatTdiffbHATold->SetY2NDC(0.9);
   v_tdiffbottom[1]->SetName("New RC values");
   Graphic_setup(v_tdiffbottom[1], 0, 0, kRed, 1, kRed, kRed, 0.2);
   v_tdiffbottom[1]->Draw("hist sames");
   gPad->Update();
   TPaveStats *pStatTdiffbHATnew = (TPaveStats *)v_tdiffbottom[1]->FindObject("stats");
   pStatTdiffbHATnew->SetTextColor(kRed);
   pStatTdiffbHATnew->SetY1NDC(0.45);
   pStatTdiffbHATnew->SetY2NDC(0.65);

   c1->cd(4);
   v_qdiffbottom[0]->SetTitle("Q_{2}/Q_{1} (bHAT)");
   v_qdiffbottom[0]->SetName("Default");
   Graphic_setup(v_qdiffbottom[0], 0, 0, kBlue, 1, kBlue, kBlue, 0.2);
   v_qdiffbottom[0]->SetMaximum(std::max(v_qdiffbottom[0]->GetMaximum(), v_qdiffbottom[1]->GetMaximum()) * 1.1);
   v_qdiffbottom[0]->Draw("hist");
   gPad->Update();
   TPaveStats *pStatQdiffbHATold = (TPaveStats *)v_qdiffbottom[0]->FindObject("stats");
   pStatQdiffbHATold->SetTextColor(kBlue);
   pStatQdiffbHATold->SetY1NDC(0.7);
   pStatQdiffbHATold->SetY2NDC(0.9);
   v_qdiffbottom[1]->SetName("New RC values");
   Graphic_setup(v_qdiffbottom[1], 0, 0, kRed, 1, kRed, kRed, 0.2);
   v_qdiffbottom[1]->Draw("hist sames");
   gPad->Update();
   TPaveStats *pStatQdiffbHATnew = (TPaveStats *)v_qdiffbottom[1]->FindObject("stats");
   pStatQdiffbHATnew->SetTextColor(kRed);
   pStatQdiffbHATnew->SetY1NDC(0.45);
   pStatQdiffbHATnew->SetY2NDC(0.65);

   c1->SaveAs("BasicRatios.pdf");

   // Delete
   delete c1;
   delete fpStyle;
   for (auto hist : v_tdiffbottom)
      delete hist;
   for (auto hist : v_tdifftop)
      delete hist;
   for (auto hist : v_qdiffbottom)
      delete hist;
   for (auto hist : v_qdifftop)
      delete hist;

   return 0;
}