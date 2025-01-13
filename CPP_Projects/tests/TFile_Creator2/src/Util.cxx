#include "Util.h"

// Convert TH2 to TGraphErrors
TGraphErrors *Convert_TH2_TGE(const TH2 *pTH2)
{
   TGraphErrors *gr = new TGraphErrors();
   for (int i = 1; i < pTH2->GetXaxis()->GetNbins(); ++i) {

      TH1D *temp_h = pTH2->ProjectionY(Form("projections_bin_%i", i), i, i);

      double x = pTH2->GetXaxis()->GetBinCenter(i);
      double y = temp_h->GetBinCenter(temp_h->GetMaximumBin());

      float start = -1.;
      float end = -1.;
      float max = temp_h->GetMaximum();

      for (Int_t bin = 0; bin < temp_h->GetXaxis()->GetNbins(); ++bin) {
         if (start == -1. && temp_h->GetBinContent(bin) >= max / 2.)
            start = temp_h->GetBinCenter(bin);

         if (end == -1. && start != -1. && temp_h->GetBinContent(bin) <= max / 2.)
            end = temp_h->GetBinCenter(bin);
      }

      float e = end - start;

      gr->SetPoint(gr->GetN(), x, y);
      gr->SetPointError(gr->GetN() - 1, 0, e / 2.);
   }
   return gr;
}

int MakeMyDir(const std::string &DirName)
{
   std::string MakeCommand = "mkdir -p " + DirName;
   return system(MakeCommand.c_str());
}