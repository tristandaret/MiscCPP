#include "LUTs.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm> // for std::find
#include <iterator>  // for std::begin, std::end

ClassImp(EramInfo)

   /* ERAM MAPS READING
      ------------------------------------------------------------------------------------------------------------------
    */
   // Default constructor
   Reconstruction::ERAMMaps::ERAMMaps()
{
   fFile = "$HOME/Documents/Code/CPP/CPP_Projects/Read_ERAM_MAP/Data/ERAMinfo.root";

   std::cout << "ERAM maps: LOADING " << fFile << std::endl;
   Load();
   std::cout << "ERAM maps: LOADED" << std::endl;
}

Reconstruction::ERAMMaps::ERAMMaps(const std::string &file)
{
   fFile = file;

   std::cout << "ERAM maps: LOADING " << fFile << std::endl;
   Load();
   std::cout << "ERAM maps: LOADED" << std::endl;
}

// Destructor
Reconstruction::ERAMMaps::~ERAMMaps()
{
   pFile->Close();
   delete pFile;
   for (int i = 0; i < (int)fID.size(); i++)
      fID[i] = 0;
   fID.clear();
   for (int i = 0; i < 34; i++) {
      for (int j = 0; j < 36; j++) {
         for (int k = 0; k < 32; k++) {
            fGain[i][j][k] = 0;
            fRC[i][j][k] = 0;
            fResolution[i][j][k] = 0;
         }
      }
   }
}

// Load the ERAM maps
void Reconstruction::ERAMMaps::Load()
{
   pFile = TFile::Open(fFile.c_str(), "READ");
   pTree = (TTree *)pFile->Get("Infotree");
   pBranch = pTree->GetBranch("EramInfo");
   int nentries = pTree->GetEntries();
   pTree->SetMakeClass(1); // Because ROOT version in ND280 is old
   pTree->SetBranchAddress("Id", &fid, &pBranch);
   pTree->SetBranchAddress("Position", &fpos, &pBranch);
   pTree->SetBranchAddress("XX", &fx, &pBranch);
   pTree->SetBranchAddress("YY", &fy, &pBranch);
   pTree->SetBranchAddress("RC", &frc, &pBranch);
   pTree->SetBranchAddress("Gain", &fgain, &pBranch);
   pTree->SetBranchAddress("Resolution", &fres, &pBranch);
   pTree->GetVal(0);

   // Initializing all values to 0
   for (int i = 0; i < 34; i++)
      for (int j = 0; j < 36; j++)
         for (int k = 0; k < 32; k++) {
            setGain(i, j, k, 0);
            setRC(i, j, k, 0);
            setResolution(i, j, k, 0);
         }

   // Filling the maps
   for (int i = 0; i < nentries; i++) {
      pTree->GetEntry(i);
      if (fid == 12) { // Not mounted @ JPARC but used in CERN22 Mockup
         setGain(32, fx, fy, fgain);
         setRC(32, fx, fy, frc);
         setResolution(32, fx, fy, fres);
      }
      if (fid == 18) { // Not mounted @ JPARC but used for CERN22 prototype
         setGain(33, fx, fy, fgain);
         setRC(33, fx, fy, frc);
         setResolution(33, fx, fy, fres);
      }
      if (fpos > 31)
         continue;
      setGain(fpos, fx, fy, fgain);
      setRC(fpos, fx, fy, frc);
      setResolution(fpos, fx, fy, fres);
   }

   // Fill holes in the maps
   FillHoles();

   if (!verbose)
      return;
   float meanGain = 0;
   float meanRC = 0;
   int iD = 0;
   for (int i = 0; i < 34; i++) {
      iD = ID(i);
      for (int j = 0; j < 36; j++) {
         for (int k = 0; k < 32; k++) {
            meanGain += Gain(i, j, k);
            meanRC += RC(i, j, k);
         }
      }
      meanGain /= 1152;
      meanRC /= 1152;
   }
}

/* Private functions */
// ERAMs on endplates 1&3 are X-flipped compared to their maps in the file
void Reconstruction::ERAMMaps::setGain(const int &position, const int &iX, const int &iY, const double &gain)
{
   fGain[position][iX][iY] = gain;
}

void Reconstruction::ERAMMaps::setRC(const int &position, const int &iX, const int &iY, const double &RC)
{
   fRC[position][iX][iY] = RC;
}

void Reconstruction::ERAMMaps::setResolution(const int &position, const int &iX, const int &iY,
                                             const double &resolution)
{
   fResolution[position][iX][iY] = resolution;
}

void Reconstruction::ERAMMaps::FillHoles()
{
   for (int i = 0; i < 32; i++) {
      for (int iY = 0; iY < 32; iY++) {
         for (int iX = 0; iX < 36; iX++) {
            double gain = Gain(i, iX, iY);
            double rc = RC(i, iX, iY);

            v_sides.assign(4, 0);
            if (gain == 0) {
               if (iX > 0)
                  v_sides[0] = Gain(i, iX - 1, iY);
               if (iX < 35)
                  v_sides[1] = Gain(i, iX + 1, iY);
               if (iY > 0)
                  v_sides[2] = Gain(i, iX, iY - 1);
               if (iY < 31)
                  v_sides[3] = Gain(i, iX, iY + 1);
               float n_sides = 0;
               for (int i = 0; i < 4; i++)
                  if (v_sides[i] != 0) { // additionnal step to discard empty neighbours
                     gain += v_sides[i];
                     n_sides++;
                  }
               gain /= n_sides;
               setGain(i, iX, iY, gain);
               // std::cout <<"ERAM#" << std::setw(2) << ID(i) << std::setw(2) << " (" << std::setw(2) << i << "): " <<
               // "Gain hole in (iX,iY) = (" << iX << "," << iY << ") | value reset at " << Gain(i,iX,iY) << std::endl;
            }

            v_sides.assign(4, 0);
            if (rc == 0) {
               if (iX > 0)
                  v_sides[0] = RC(i, iX - 1, iY);
               if (iX < 35)
                  v_sides[1] = RC(i, iX + 1, iY);
               if (iY > 0)
                  v_sides[2] = RC(i, iX, iY - 1);
               if (iY < 31)
                  v_sides[3] = RC(i, iX, iY + 1);
               float n_sides = 0;
               for (int i = 0; i < 4; i++)
                  if (v_sides[i] != 0) { // additionnal step to discard empty neighbours
                     rc += v_sides[i];
                     n_sides++;
                  }
               rc /= n_sides;
               setRC(i, iX, iY, rc);
               // std::cout <<"ERAM#" << std::setw(2) << ID(i) << std::setw(2) << " (" << std::setw(2) << i << "): " <<
               // "RC   hole in (iX,iY) = (" << iX << "," << iY << ") | value reset at " << RC(i,iX,iY) << std::endl;
            }
            v_sides.clear();
         } // iY
      }    // iX
   }
}