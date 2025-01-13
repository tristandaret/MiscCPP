#include <iomanip>

#include <TROOT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "Util.h"
#include "Misc_Functions.h"
#include "SetStyle.h"

#include "TrackModel.h"
#include "LUTMaker.h"

void modifyBranches();

int main()
{
   TrackModel trackmodel;
   LUTMaker lutmaker;

   // => Exact: 356.422 | Amplitude: 358 | length: 10.1899 | d: -0.620518 | phi: -89.9624 | RC: 112 | X: 489
   // Approx: 382.929 | Scalefactor: 1.06963 | relative difference: 7.43679
   float amplitude = 358;
   float Dt = 310 / pow(10, 3.5);
   float drift = 489;
   float RC = 112;
   float impact = -0.620518;
   float phi = -89.9624;
   float length = 10.1899;
   // => Exact: 744.208 | Amplitude: 697 | length: 10.1899 | d: -0.618579 | phi: 89.9842 | RC: 112 | X: 485.9
   // Approx: 745.181 | Scalefactor: 1.06913 | relative difference: 0.130738
   // float amplitude = 697;
   // float Dt = 310/pow(10, 3.5);
   // float drift = 485.9;
   // float RC = 112;
   // float impact = -0.618579;
   // float phi = 89.9842;
   // float length = 10.1899;

   float diag = lutmaker.GetDiag();

   lutmaker.ComputeLengthMap();
   // lutmaker.DrawLengthMap();
   // lutmaker.MakeLUT();
   // lutmaker.DrawLUT(120, 0, 310/pow(10, 3.5));

   // lutmaker.LoadLUT("Output_LUT/dEdx_XP_LUT_v2.root");
   // float scalefactorLUT = lutmaker.GetFactorFromLUT(Dt, RC, drift, fabs(impact), fabs(phi));
   // float recochargeLUT = amplitude*scalefactorLUT;

   // trackmodel.ComputeRealCharge(amplitude, length, fabs(impact), fabs(phi), RC, drift, Dt);
   // float recochargeModel = trackmodel.GetRealCharge();
   // float scalefactorModel = 1/trackmodel.GetAmplitudeLoss();

   // std::cout << "Amplitude:                    " << amplitude << std::endl;
   // std::cout << "Dt:                           " << Dt << std::endl;
   // std::cout << "Drift:                        " << drift << std::endl;
   // std::cout << "RC:                           " << RC << std::endl;
   // std::cout << "Impact:                       " << impact << std::endl;
   // std::cout << "Phi:                          " << phi << std::endl;
   // std::cout << "Length:                       " << length << std::endl;
   // std::cout << "Scale factor LUT:             " << scalefactorLUT << std::endl;
   // std::cout << "Scale factor Model:           " << scalefactorModel << std::endl;
   // std::cout << "Relative difference:          " << (scalefactorModel-scalefactorLUT)/scalefactorLUT*100 << "%" <<
   // std::endl; std::cout << "Reconstructed charge (LUT):   " << recochargeLUT << std::endl; std::cout <<
   // "Reconstructed charge (Model): " << recochargeModel << std::endl; std::cout << "Relative difference:          " <<
   // (recochargeModel-recochargeLUT)/recochargeLUT*100 << "%" << std::endl;

   // lutmaker.DrawDiffInterpolExact(Dt, RC, drift, 0, 90, 0, diag/2);

   lutmaker.PrintLUTcorners();

   std::cout << lutmaker.GetDiag() << std::endl;
   for (int i = 0; i < lutmaker.GetVecPhi().size(); i++)
      std::cout << i << " " << lutmaker.GetVecPhi()[i] << std::endl;

   // modifyBranches();

   return 0;
}

void modifyBranches()
{
   TFile *file = TFile::Open("Output_LUT/LUT.root", "READ");
   TTree *tree = (TTree *)file->Get("LUT");

   // Deactivate branches that should not be copied
   tree->SetBranchStatus("length", 0);
   tree->SetBranchStatus("d", 0);
   tree->SetBranchStatus("phi", 0);
   tree->SetBranchStatus("Dt", 0);

   // Create a new tree by cloning the original tree without the deactivated branches
   TFile *newFile = TFile::Open("Output_LUT/dEdx_XP_LUT_v2.root", "RECREATE");
   TTree *newTree = tree->CloneTree(0);

   // Reactivate the branches to be renamed and set their addresses
   tree->SetBranchStatus("d", 1);
   tree->SetBranchStatus("phi", 1);
   tree->SetBranchStatus("Dt", 1);
   double d, phi;
   float Dt;
   tree->SetBranchAddress("d", &d);
   tree->SetBranchAddress("phi", &phi);
   tree->SetBranchAddress("Dt", &Dt);

   // Create new branches in the new tree
   double impact, angle;
   float transvDiff;
   newTree->Branch("impact", &impact);
   newTree->Branch("angle", &angle);
   newTree->Branch("transvDiff", &transvDiff);

   // Copy data from old branches to new branches
   Long64_t nentries = tree->GetEntries();
   for (Long64_t i = 0; i < nentries; i++) {
      tree->GetEntry(i);
      impact = d;
      angle = phi;
      transvDiff = Dt;
      newTree->Fill();
   }

   // Write the new tree to the new file
   newTree->Write("treeLUTdEdx");
   newFile->Close();
   file->Close();

   delete file;
   delete newFile;
}