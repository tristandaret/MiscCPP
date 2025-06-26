#ifndef LUT_MAKER_H
#define LUT_MAKER_H

#include <cmath>
#include <vector>
#include <string>

#include "SetStyle.h"

#include "TrackModel.h"

class LUTMaker {
public:
   LUTMaker();
   virtual ~LUTMaker();

   // Geometry
   float GetX(const double &y, const double &phi_rad, const double &d)
   {
      return (y - (d - sin(phi_rad) * SXWIDTH / 2 + cos(phi_rad) * SYHEIGHT / 2) /
                     cos(phi_rad)) /
             tan(phi_rad);
   }
   float GetY(const double &x, const double &phi_rad, const double &d)
   {
      return tan(phi_rad) * x +
             (d - sin(phi_rad) * SXWIDTH / 2 + cos(phi_rad) * SYHEIGHT / 2) /
                cos(phi_rad);
   }
   float Slope(const double &phi_rad) { return tan(phi_rad); }
   float Intercept(const double &phi_rad, const double &d)
   {
      return (cos(phi_rad) * SYHEIGHT / 2 - sin(phi_rad) * SXWIDTH / 2 + d) /
             cos(phi_rad);
   }

   double ComputeLength(const double &phi_rad, const double &d);
   void ComputeLengthMap();
   void MakeLUT();

   void LoadLUT(std::string LUTpath);
   float GetFactorFromLUT(const double &transvDiff, const double &RC, const double &drift,
                          const double &d, const double &phi);

   // Getters
   std::vector<std::vector<float>> GetLengthMap() { return arr_length; }
   std::vector<double> GetVecPhi() { return v_phi; }
   std::vector<double> GetVecD() { return v_impact; }
   std::vector<double> GetVecZ() { return v_drift; }
   float GetDiag() { return PAD_DIAG; }

   // Setters
   void SetPeakingtime(const double &peakingTime)
   {
      fp_trackmodel->SetPeakingTime(peakingTime);
   }
   void SetTransDiff(const double &transvDiffwithB, const double &transvDiffwithoutB)
   {
      v_Dt.clear();
      v_Dt.push_back(transvDiffwithB / pow(10, 3.5)); // with B field
      v_Dt.push_back(transvDiffwithoutB / pow(10, 3.5)); // without B field
      stepSizeTrans = fabs(transvDiffwithB - transvDiffwithoutB);
   }

   // Drawing
   void DrawDiffInterpolExact(const double &transvDiff, const double &RC,
                              const double &drift, const double &phimin,
                              const double &phimax, const double &dmin,
                              const double &dmax);
   void DrawLengthMap();
   void DrawLUT(const double &RC, const double &drift, const double &transvDiff);

   // Debugging
   void PrintSettings();

private:
   // Graphics
   TCanvas *c1;
   TStyle *ptstyle;

   // Pad dimensions
   static constexpr float SXWIDTH = 11.28;  // mm
   static constexpr float SYHEIGHT = 10.19; // mm
   // const float PAD_DIAG = 				std::sqrt(std::pow(SXWIDTH, 2) + std::pow(SYHEIGHT,
   // 2));
   static constexpr float PAD_DIAG = 15.2011; // sqrt(pow(11.28,2) + pow(10.19, 2))

   // Track model
   TrackModel *fp_trackmodel;

   // Vectors of discrete values for each dimension
   std::vector<double> v_phi, v_impact, v_drift, v_RC, v_Dt;

   static constexpr int RCmin = 80;
   static constexpr int RCmax = 160;
   static constexpr int DtwithB=286;
   static constexpr int DtwithoutB=323;
   // Number of discrete steps in each dimension of the Look Up Table
   static const int SNSTEPS_TRANS = 2;
   static const int SNSTEPS_RC = 17;
   static const int SNSTEPS_DRIFT = 21;
   static const int SNSTEPS_D = 250;
   static const int SNSTEPS_PHI = 250;

   // Length map
   std::vector<std::vector<float>> arr_length =
      std::vector<std::vector<float>>(SNSTEPS_PHI, std::vector<float>(SNSTEPS_D, 0));

   // LUT loader
   std::string fLUTpath = "Output_LUT/LUT_test.root";
   // float LUTValues[SNSTEPS_PHI][SNSTEPS_D];
   static float LUTValues[SNSTEPS_TRANS][SNSTEPS_RC][SNSTEPS_DRIFT][SNSTEPS_D]
                         [SNSTEPS_PHI];
   float stepSizeTrans = 37; // 323->286, only 2 values
   static constexpr float sSTEP_RC = 5;
   static constexpr float sSTEP_PHI = 90. / (SNSTEPS_PHI - 1);
   static constexpr float sSTEP_IMPACT = (PAD_DIAG / 2) / (SNSTEPS_D - 1);
   static constexpr float sSTEP_DRIFT = 1000. / (SNSTEPS_DRIFT - 1);
};

#endif