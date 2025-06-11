#include <iostream>

#include <TCanvas.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TPaveStats.h>
#include <TStyle.h>

int main()
{
   TRandom3 randGen;

   // Create a canvas
   TCanvas *canvas = new TCanvas("canvas", "Canvas for Histograms", 800, 600);

   // Create histograms
   TH1D hGaussian1("hGaussian1", "Gaussian Distribution 1;X;Frequency", 100, -5, 5);
   TH1D hNormal("hNormal", "Normal Distribution;X;Frequency", 100, -5, 5);
   TH1D hNormal2("hNormal2", "Normal Distribution 2;X;Frequency", 100, -5, 5);
   TH1D hGaussian2("hGaussian2", "Gaussian Distribution 2;X;Frequency", 100, -5, 5);

   // Fill histograms with random data
   for (int i = 0; i < 1000; ++i) {
      hGaussian1.Fill(randGen.Gaus(0, 1));   // Mean 0, StdDev 1
      hNormal.Fill(randGen.Uniform(-5, 5));  // Uniform distribution
      hNormal2.Fill(randGen.Uniform(-5, 5)); // Uniform distribution
      hGaussian2.Fill(randGen.Gaus(2, 0.5)); // Mean 2, StdDev 0.5
   }

   gStyle->SetOptStat(111111);

   std::string outputFileName = "histograms.pdf";
   hGaussian1.Draw();
   canvas->SaveAs((outputFileName + "(").c_str());

   hNormal.SetLineColor(kRed);
   hNormal.Draw("hist");
   gPad->Update();
   TPaveStats *stats1 = (TPaveStats *)hNormal.GetListOfFunctions()->FindObject("stats");
   stats1->SetY1NDC(0.3);
   stats1->SetY2NDC(0.5);

   hNormal2.SetLineColor(kBlue);
   hNormal2.Draw("SAMES");
   gPad->Update();
   TPaveStats *stats2 =
   (TPaveStats*)hNormal2.GetListOfFunctions()->FindObject("stats");
   stats2->SetY1NDC(0.5);
   stats2->SetY2NDC(0.7);
   canvas->SaveAs(outputFileName.c_str());

   hGaussian2.Draw();
   canvas->SaveAs((outputFileName + ")").c_str());

   return 0;
   // Clean up
   delete canvas;
   // delete stats;
}
