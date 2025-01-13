#ifndef LUT_H
#define LUT_H

#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class EramInfo : public TObject {
public:
   int Id;
   int Position;
   int XX;
   int YY;
   double RC;
   double Gain;
   double Resolution;

   int Endplate;
   bool InbTPC;
   bool IntTPC;

   EramInfo() {}
   virtual ~EramInfo() {}

   ClassDef(EramInfo, 2)
};

namespace Reconstruction {

/* ERAM MAPS READING
 * ---------------------------------------------------------------------------------------------------------------------------
 */
class ERAMMaps {
public:
   /* Constructor */
   ERAMMaps();
   ERAMMaps(const std::string &file);
   virtual ~ERAMMaps();

   int ID(const int &position) { return channel2iD[position]; }
   double RC(const int &position, const int &iX, const int &iY) { return fRC[position][iX][iY]; }
   double Gain(const int &position, const int &iX, const int &iY) { return fGain[position][iX][iY]; }
   double Resolution(const int &position, const int &iX, const int &iY) { return fResolution[position][iX][iY]; }

private:
   bool verbose = true;
   std::string fFile;
   std::vector<int> fID;
   std::vector<float> v_sides;

   double fGain[34][36][32]; // [position][iX][iY]
   double fRC[34][36][32];
   double fResolution[34][36][32];

   void Load();
   void setGain(const int &position, const int &iX, const int &iY, const double &gain);
   void setRC(const int &position, const int &iX, const int &iY, const double &RC);
   void setResolution(const int &position, const int &iX, const int &iY, const double &resolution);
   void FillHoles();

   std::vector<int> channel2iD = {24, 30, 28, 19, 21, 13, 9,  2,  26, 17, 23, 29, 1,  10, 11, 3,  /*bottom HATPC*/
                                  47, 16, 14, 15, 42, 45, 37, 36, 20, 38, 7,  44, 43, 39, 41, 46, /*top    HATPC*/
                                  12, 18}; // CERN22 MockUp and prototype

   //------------------------------Data Members-----------------------//
   TFile *pFile;
   TTree *pTree;
   TBranch *pBranch;

   // Leaves
   int fid;
   int fpos;
   int fx;
   int fy;
   double frc;
   double fgain;
   double fres;
};
} // namespace Reconstruction

#endif