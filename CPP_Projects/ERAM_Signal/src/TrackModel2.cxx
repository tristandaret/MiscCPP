#include "TrackModel2.h"
#include "SetStyle.h"

// Constructor
TrackModel2::TrackModel2(double peakingTime, double xpad, double ypad, double phi,
                         double impact, double drift, double RC)
{
   SetPeakingTime(peakingTime);
   SetPadParam(xpad, ypad, RC);
   SetTrackParams(phi, impact, drift);
   SetFunctions();
}

// Destructor
TrackModel2::~TrackModel2()
{
   delete ptf1_dETFdt;
   delete ptf1_Charge;
   delete ptf1_Convolution;
   delete ptf1_Signal;
}

// Settings

void TrackModel2::SetTrackParams(double phi, double impact, double drift)
{
   fdrift = drift;
   fimpact = impact;
   // Clamp phi to avoid numerical instabilities
   fphi = (phi < 1e-5) ? 1e-5 : ((phi > 90 - 1e-5) ? (90 - 1e-5) : phi);
   fphirad = phi / 180 * M_PI;
   // Safeguards for numerical instabilities
   if (phi < 1e-5)
      fphirad = 1e-5 / 180 * M_PI;
   if (phi > 90 - 1e-5)
      fphirad = (90 - 1e-5) / 180 * M_PI;
   // Set parameters
   fm = std::tan(fphirad);
   fq = (cos(fphirad) * fycenter - sin(fphirad) * fxcenter + fimpact) / cos(fphirad);
}

void TrackModel2::SetFunctions(){
   // Initialize the transfer function
   ptf1_dETFdt = new TF1("ptf1_dETFdt", this, &TrackModel2::dETFdt, 0, ftmax, 0);
   ptf1_dETFdt->SetNpx(1000);
   ptf1_dETFdt->SetLineColor(kRed);
   ptf1_dETFdt->SetLineWidth(2);

   // Initialize the charge function
   ptf1_Charge = new TF1("ptf1_Charge", this, &TrackModel2::Charge, 0, ftmax, 0);
   ptf1_Charge->SetNpx(1000);
   ptf1_Charge->SetLineColor(kBlue);
   ptf1_Charge->SetLineWidth(2);

   // Initialize the convolution function
   ptf1_Convolution = new TF1Convolution(ptf1_dETFdt, ptf1_Charge, -ftmax, ftmax, true);
   ptf1_Convolution->SetNofPointsFFT(fnconvolpoints);
   ptf1_Signal =
      new TF1("", ptf1_Convolution, 0, ftmax, ptf1_Convolution->GetNpar());
}



// Eletronic transfer function
Double_t TrackModel2::dETFdt(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t <= 0)
      return 0;
   else
      return fnormelec *
             (-ws * std::exp(-ws * t) +
              std::exp(-expfactor * t) *
                 (expfactor * (std::cos(arg * t) - sinfactor * std::sin(arg * t)) +
                  arg * (std::sin(arg * t) + sinfactor * std::cos(arg * t))));
}

// Charge function
Double_t TrackModel2::Charge(Double_t *x, Double_t *par)
{
   Double_t t = x[0];

   if (t <= 0)
      return 0;
   // Factors
   double sigma =
      std::sqrt(2 * t / fRC + fDt * fDt * fdrift); // includes transverse diffusion
   double coeff1 = std::sqrt(2 * (1 + fm * fm) / M_PI) * sigma;
   double denom = 2 * (1 + fm * fm) * sigma * sigma;
   double sqrt2_m2_sigma = std::sqrt(2 * (1 + fm * fm)) * sigma;

   // Terms
   double expRC1 = std::exp(-std::pow(-fylow + fxright * fm + fq, 2) / denom);
   double expRC2 = std::exp(-std::pow(-fylow + fxleft * fm + fq, 2) / denom);
   double expRC3 = std::exp(-std::pow(-fyhigh + fxleft * fm + fq, 2) / denom);
   double expRC4 = std::exp(-std::pow(-fyhigh + fxright * fm + fq, 2) / denom);

   double erfRC1 =
      (fylow - fxleft * fm - fq) * std::erf((-fylow + fxleft * fm + fq) / sqrt2_m2_sigma);
   double erfRC2 = (fyhigh - fxleft * fm - fq) *
                   std::erf((-fyhigh + fxleft * fm + fq) / sqrt2_m2_sigma);
   double erfRC3 = (fylow - fxright * fm - fq) *
                   std::erf((-fylow + fxright * fm + fq) / sqrt2_m2_sigma);
   double erfRC4 = (fyhigh - fxright * fm - fq) *
                   std::erf((-fyhigh + fxright * fm + fq) / sqrt2_m2_sigma);

   // Return result
   return linearCharge * std::sqrt(1 + fm * fm) / (2 * fm) *
          (coeff1 * (expRC1 - expRC2 + expRC3 - expRC4) + erfRC1 - erfRC2 - erfRC3 +
           erfRC4);
}