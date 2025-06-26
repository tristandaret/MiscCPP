#include <iostream>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPaveStats.h>
#include <TAxis.h>
#include <TLegendEntry.h>
#include <TStyle.h>

#include "ROOT_file.h"
#include "Tree_classes.h"

void Reader(const char *filename)
{
   TROOT root("MyRoot", "My Root Title");
   TFile *inputFile = TFile::Open(filename, "READ");

   if (!inputFile || inputFile->IsZombie()) {
      cerr << "Error: Cannot open input file." << endl;
      return;
   }

   TTree *treeT = dynamic_cast<TTree *>(inputFile->Get("treeT"));

   if (!treeT) {
      cerr << "Error: Cannot find treeT in the input file." << endl;
      inputFile->Close();
      return;
   }

   EventT *eventT = nullptr;
   treeT->SetBranchAddress("events", &eventT);

   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   std::string OutputFile = "../Evt_Display.pdf";
   std::string OutputFile_Beg = OutputFile + "(";
   std::string OutputFile_End = OutputFile + ")";
   TCanvas *pTCanvas = new TCanvas("Evt_Displays", "Evt_Displays", 1800, 1200);
   pTCanvas->cd();

   const int n_evt = treeT->GetEntries();
   int i_evt = 0;
   for (int entry = 0; entry < treeT->GetEntries(); entry++) {
      treeT->GetEntry(entry);
      std::cout << "Event's ID: #" << eventT->eventid << endl;
      int i_subtrk = 0;

      for (const TrackT &track : eventT->tracks) {
         std::cout << "  Track's ID: #" << track.trackid << " | " << track.n_subtracks
                   << " subtracks" << endl;
         std::cout << "  Track's parameters [mm]: X = " << track.params[0]
                   << " | Y = " << track.params[1] << " | Z = " << track.params[2]
                   << " | Vx = " << track.params[3] << " | Vy = " << track.params[4]
                   << " | Vz = " << track.params[5] << " | Curve = " << track.params[6]
                   << endl;

         for (const SubtrackT &subtrack : track.subtracks) {
            std::cout << "    Subtrack's ERAM ID: #" << subtrack.moduleid << endl;
            TGraph *tg_subtrk = new TGraph();
            std::string name_tg = Form("tg_evt%i_trk%i_sbtk%i", eventT->eventid,
                                       track.trackid, subtrack.moduleid);
            tg_subtrk->SetNameTitle(name_tg.c_str(), name_tg.c_str());

            int iC = 0;
            float avg_x = 0;
            for (const ClusterT &cluster : subtrack.clusters) {
               std::cout << "      Cluster's position: X = " << cluster.x
                         << " | Y = " << cluster.y << " | Z = " << cluster.z << endl;
               tg_subtrk->SetPoint(iC, cluster.z, cluster.y);
               avg_x += cluster.x;
               iC++;
            }
            avg_x /= (iC + 1);

            tg_subtrk->SetMarkerSize(2);
            tg_subtrk->SetMarkerColor(i_subtrk + 2);
            if (avg_x > 0)
               tg_subtrk->SetMarkerStyle(20);
            else
               tg_subtrk->SetMarkerStyle(24);
            tg_subtrk->GetXaxis()->SetLimits(-3000, -1800);
            tg_subtrk->SetMinimum(300);
            tg_subtrk->SetMaximum(1200);
            if (i_subtrk == 0) {
               tg_subtrk->SetNameTitle(
                  Form("tg_evt%i", eventT->eventid),
                  Form("Event %i Display;Z [mm];Y [mm]", eventT->eventid));
               tg_subtrk->DrawClone("ap");
            } else
               tg_subtrk->DrawClone("p same");
            delete tg_subtrk;
            tg_subtrk = nullptr;
            i_subtrk++;
         }
      }

      if (i_evt == 0)
         pTCanvas->SaveAs(OutputFile_Beg.c_str());
      else if (i_evt == n_evt - 1)
         pTCanvas->SaveAs(OutputFile_End.c_str());
      else
         pTCanvas->SaveAs(OutputFile.c_str());
      std::cout << " " << endl;
      i_evt++;
   }
   inputFile->Close();
}