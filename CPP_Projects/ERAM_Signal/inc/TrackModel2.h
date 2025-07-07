#ifndef TRACK2_MODEL_H
#define TRACK2_MODEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <TMath.h>
#include <TF1.h>
#include <TF1Convolution.h>

class TrackModel2 {
public:
   TrackModel2(double peakingTime, double xpad, double ypad, double phi, double impact,
               double drift, double RC);
   virtual ~TrackModel2();

   // Set parameters
   void SetPeakingTime(const double &peakingTime)
   {
      fpeakingTime = peakingTime;
      ws = 2 / fpeakingTime;
      expfactor = ws / (2 * Q);
      arg = ws / 2 * std::sqrt(4 - 1 / std::pow(Q, 2));
   }

   void SetPadParam(double xpad, double ypad, double RC)
   {
      fxpad = xpad;
      fypad = ypad;
      fxleft = fxpad - fSXWIDTH / 2;
      fxright = fxpad + fSXWIDTH / 2;
      fylow = fypad - fSYWIDTH / 2;
      fyhigh = fypad + fSYWIDTH / 2;
      fRC = RC;
   }
   void SetTrackParams(double phi, double impact, double drift);
   void SetFunctions();

   // Getters
   double GetSlope() const { return fm; }
   double GetIntercept() const { return fq; }
   double GetAmplitude() const
   {
      return ptf1_Signal->GetMaximum();
   }
   double GetMinimum() const
   {
      return ptf1_Signal->GetMinimum();
   }

   // Signal functions
   Double_t dETFdt(Double_t *x, Double_t *par);
   Double_t Charge(Double_t *x, Double_t *par);

   TF1 *ptf1_dETFdt;
   TF1 *ptf1_Charge;
   TF1 *ptf1_Signal;
   TF1Convolution *ptf1_Convolution;

private:
   double ftmax = 3200;              // Maximum time

   // Eletronics transfer function
   double fpeakingTime = 412.;
   const double Q = 2. / 3.;
   double ws = 2 / fpeakingTime;
   double expfactor = ws / (2 * Q);
   double arg = ws / 2 * std::sqrt(4 - 1 / std::pow(Q, 2));
   const double sinfactor = std::sqrt((2 * Q - 1) / (2 * Q + 1));
   const double fnormelec =
      4096. / 120. / 0.16233963; // 4096 max ADC for 120 fC | 0.16... max of ETF
   double ETFmax = 0;

   // Track parameters
   double fphi = 0, fphirad = 0, fimpact = 0, fdrift = 0;
   double fm = 0, fq = 0; // Slope and charge for the track
   double fDt = 286 / pow(10, 3.5);

   constexpr static double Ne = 10;              // electrons per cm
   constexpr static double qe = 1.602176634e-4; // Elementary charge in fC
   constexpr static double gain = 1493.0;        // Gain factor for the signal

   // Charge function
   double fxpad = 0, fypad = 0, fRC = 0;
   double fxleft = 0, fxright = 0, fylow = 0, fyhigh = 0;
   double linearCharge = Ne * qe * gain; // Linear charge in Coulombs per cm
   constexpr static double fSXWIDTH = 11.28;
   constexpr static double fSYWIDTH = 10.19;
   double fxcenter = fSXWIDTH / 2; 
   double fycenter = fSYWIDTH / 2; // Center of the pad
   
   // Convolution
   double fnconvolpoints = 1e4 + 1;
};

#endif