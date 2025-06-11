#ifndef POINT_MODEL_H
#define POINT_MODEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <TMath.h>
#include <TF1.h>
#include <TF1Convolution.h>

void DrawPointModel();

class PointModel {
public:
   PointModel(float peakingTime, float xpad, float ypad, float xcharge, float ycharge, float RC, float Qdrop = Q0);
   virtual ~PointModel();

   // Set parameters
   void SetPeakingTime(float peakingTime);
   void SetPadParam(float xpad, float ypad, float RC)
   {
      fxpad = xpad;
      fypad = ypad;
      fRC = RC;
   }
   void SetChargeParam(float xcharge, float ycharge, float Qdrop)
   {
      fxcharge = xcharge;
      fycharge = ycharge;
      fQdrop = Qdrop;
   }

   // Signal functions
   Double_t ETF(Double_t *x, Double_t *par);
   Double_t ETFscaled(Double_t *x, Double_t *par);
   Double_t dETFdt(Double_t *x, Double_t *par);
   Double_t Charge(Double_t *x, Double_t *par);
   Double_t Current(Double_t *x, Double_t *par);
   Double_t CurrentScaled(Double_t *x, Double_t *par);
   Double_t CurrentnoDirac(Double_t *x, Double_t *par);

   TF1 *ptf1_ETF;
   TF1 *ptf1_ETFscaled;
   TF1 *ptf1_dETFdt;
   TF1 *ptf1_Charge;
   TF1 *ptf1_Current;
   TF1 *ptf1_CurrentScaled;
   TF1 *ptf1_CurrentnoDirac;
   TF1 *ptf1_Signal;
   TF1 *ptf1_SignalnoDirac;
   TF1 *ptf1_Signal2;
   TF1 *ptf1_Signal3;

   TF1Convolution *ptf1_Convolution;
   TF1Convolution *ptf1_ConvolutionnoDirac;
   TF1Convolution *ptf1_Convolution2;
   TF1Convolution *ptf1_Convolution3;

private:
   float ftmax = 5e3; // Maximum time
   constexpr static float Q0 = 54; // Charge in fC

   // Eletronics transfer function
   float fpeakingTime = 412;
   float ws;
   float arg;
   float expfactor;

   constexpr static float Q = 2. / 3.;
   const float sinfactor = std::sqrt((2 * Q - 1) / (2 * Q + 1));

   constexpr static float fSETFMAX = 0.162345;

   // Charge function
   float fxpad = 0, fypad = 0, fxcharge = 0, fycharge = 0, fRC = 0, fQdrop = 0;
   float fxleft = 0, fxright = 0, fylow = 0, fyhigh = 0;
   constexpr static float fSXWIDTH = 11.28;
   constexpr static float fSYWIDTH = 10.19;
   float erfx, erfy;

   // Convolution
   float fnconvolpoints = 1e4 + 1;
};

#endif