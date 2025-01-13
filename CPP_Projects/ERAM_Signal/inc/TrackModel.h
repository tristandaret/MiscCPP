#ifndef TRACK_MODEL_H
#define TRACK_MODEL_H

#include <iostream>
#include <cmath>
#include <TF1.h>

class TrackModel {
public:
   TrackModel();
   virtual ~TrackModel();

   // Initialize parameters and variables
   void SetParameters(const double &length, const double &impact, const double &phi, const double &RC,
                      const double &drift, const double &Dt);
   void SetVariables(const double &signal);
   void SetNConvPoints(const double &nconvpoints) { fnconvolpoints = nconvpoints; }
   void SetTMax(const double &tmax) { ftmax = tmax; }

   // Compute the signal amplitude loss due to charge spreading
   void ComputeAmplitudeLoss();
   void ComputeAmplitudeLoss(const double &length, const double &impact, const double &phi, const double &RC,
                             const double &drift, const double &Dt);
   double GetAmplitudeLoss() const { return famplitudeloss; }

   // Recompute the deposited charge based on model
   // distances in mm, angles in degrees
   void ComputeRealCharge();
   void ComputeRealCharge(const double &ADCmax, const double &length, const double &impact, const double &phi,
                          const double &RC, const double &drift, const double &Dt);
   double GetRealCharge() const { return fchargereal; }

   double time = 0;

   // Debug
   double ComputeLength(const double &phi_rad, const double &d);
   float GetX(const double &y, const double &phi_rad, const double &d)
   {
      return (y - (d - sin(phi_rad) * SXWIDTH / 2 + cos(phi_rad) * SYHEIGHT / 2) / cos(phi_rad)) / tan(phi_rad);
   }
   float GetY(const double &x, const double &phi_rad, const double &d)
   {
      return tan(phi_rad) * x + (d - sin(phi_rad) * SXWIDTH / 2 + cos(phi_rad) * SYHEIGHT / 2) / cos(phi_rad);
   }

private:
   // Convolution
   double fnconvolpoints = 750;
   double ftmax = 750; // Maximum time

   // Eletronics transfer function
   Double_t ETF(Double_t *x, Double_t *par);
   Double_t dETFdt(Double_t *x, Double_t *par);
   TF1 *ptf1_dETFdt;
   double fpeakingTime = 412.;
   const double Q = 2. / 3.;
   const double ws = 2 / fpeakingTime;
   const double expfactor = ws / (2 * Q);
   const double arg = ws / 2 * std::sqrt(4 - 1 / std::pow(Q, 2));
   const double sinfactor = std::sqrt((2 * Q - 1) / (2 * Q + 1));
   const double fnormelec = 4096. / 120. / 0.16233963; // 4096 max ADC for 120 fC | 0.16... max of ETF
   double ETFmax = 0;

   // Charge function
   Double_t Charge(Double_t *x, Double_t *par);
   // Track inputs
   double fDt = 0; // conversion to sqrt(mm) Transverse diffusion coefficent
   double fRC = 0;
   double flength = 0;
   double fimpact = 0;
   double fphirad = 0;
   double fdrift = 0;
   double fm = 0; // slope
   double fq = 0; // intercept
   // Computation intermediate variables
   double fsignalreal = 0;
   double fchargemodel = 0;
   double fsignalmodel = 0;
   double famplitudeloss = 0;
   double fchargereal = 0;
   // Pad parameters
   const double fxpad = 0.; // pad x position
   const double fypad = 0.; // pad y position
   const double fxleft = -11.28 / 2;
   const double fxright = 11.28 / 2;
   const double fylow = -10.19 / 2;
   const double fyhigh = 10.19 / 2;

   // Debug
   static constexpr float SXWIDTH = 11.28;  // mm
   static constexpr float SYHEIGHT = 10.19; // mm
   std::vector<int> fRCvalues = {112, 158};
   std::vector<int> fRCtype = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  /*bHAT*/
                               1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1}; /*tHAT*/
   static constexpr float PAD_DIAG = 15.2011;                                   // sqrt(pow(11.28,2) + pow(10.19, 2))
   static const int SNSTEPS_TRANS = 2;
   static const int SNSTEPS_RC = 2;
   static const int SNSTEPS_DRIFT = 101;
   static const int SNSTEPS_D = 250;
   static const int SNSTEPS_PHI = 250;

   static constexpr float sSTEP_TRANS = 40;                    // 310->350, only 2 values
   static constexpr float sSTEP_RC = 46;                       // 112->158, only 2 values
   static constexpr float sSTEP_PHI = 90. / (SNSTEPS_PHI - 1); // small shift because of numerical instabilities
   static constexpr float sSTEP_D = (PAD_DIAG / 2) / (SNSTEPS_D - 1);
   static constexpr float sSTEP_DRIFT = 1000. / (SNSTEPS_DRIFT - 1);
};

#endif