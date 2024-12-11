#include "TrackModel.h"

#include <TROOT.h>
#include <TF1Convolution.h>

// Constructor
TrackModel::TrackModel(){
    TF1 tf1_ETF(";time (ns); ADC counts/ns",   this, &TrackModel::ETF,  0, ftmax, 0);
    ptf1_dETFdt =		        new TF1(";time (ns); ADC counts/ns",   this, &TrackModel::dETFdt,  0, ftmax, 0);
    ETFmax =                    tf1_ETF.GetMaximum();
}


// Destructor
TrackModel::~TrackModel(){
    delete ptf1_dETFdt;
}


// Initialize parameters
void TrackModel::SetParameters( const double &length, const double &impact, 
                                const double &phi, const double &RC, const double &drift, const double &Dt){
    fDt =                       Dt;
    fRC =                       RC;
    fdrift =                    drift;
    flength =                   length;
    fimpact =                   impact;
    fphirad =                   phi/180*M_PI;
    // Safeguards for numerical instabilities
    if(phi < 1e-5) fphirad =    1e-5/180*M_PI;
    if(phi > 90-1e-5) fphirad = (90-1e-5)/180*M_PI;
    // Set parameters
    fm =                        std::tan(fphirad);
    fq =                        (cos(fphirad)*fxpad-sin(fphirad)*fxpad+fimpact)/cos(fphirad);
}

// Initialize real signal
void TrackModel::SetVariables(const double &signal){
    fsignalreal =               signal;
}


void TrackModel::ComputeAmplitudeLoss(  const double &length, const double &impact, 
                                        const double &phi, const double &RC, const double &drift, const double &Dt){
    SetParameters(length, impact, phi, RC, drift, Dt);
    ComputeAmplitudeLoss();
}
void TrackModel::ComputeAmplitudeLoss() {
    // Initialize functions
    TF1 tf1_Charge(";time (ns); Charge (fC)", this, &TrackModel::Charge, 0, ftmax, 0);
    TF1Convolution tf1_Convolution(&tf1_Charge, ptf1_dETFdt, -ftmax, ftmax, true);
    tf1_Convolution.SetNofPointsFFT(fnconvolpoints);
    TF1 tf1_Signal("", tf1_Convolution, 0, ftmax, tf1_Convolution.GetNpar());
    ptf1_Signal = new TF1(tf1_Signal);

    fchargemodel = flength * ETFmax;
    fsignalmodel = tf1_Signal.GetMaximum();
    famplitudeloss = fsignalmodel / fchargemodel;
}



void TrackModel::ComputeRealCharge(const double &ADCmax, const double &length, const double &impact, 
                                 const double &phi, const double &RC, const double &drift, const double &Dt){
    SetParameters(length, impact, phi, RC, drift, Dt);
    SetVariables(ADCmax);
    ComputeRealCharge();
}
void TrackModel::ComputeRealCharge(){
    ComputeAmplitudeLoss();

    // Compute charge
    /* The real charge deposited over the pad is linearly proportional to the length in the pad.
    A normalization has to be chosen; since the data signal depends on the electronics transfer 
    function, its max is a convenient choice */
    fchargereal =               fsignalreal/famplitudeloss;
}


// Eletronic transfer function
Double_t TrackModel::ETF(Double_t *x, Double_t *par){
    Double_t t = x[0];

	if(t <= 0) return 0;
	else return fnormelec*(std::exp(-ws*t) + std::exp(-expfactor*t) * (sinfactor*std::sin(arg*t) - std::cos(arg*t)));
}

Double_t TrackModel::dETFdt(Double_t *x, Double_t *par){
    Double_t t = x[0];

    if(t <= 0) return 0;
    else return fnormelec*(-ws*std::exp(-ws*t) + std::exp(-expfactor*t) * (expfactor*(std::cos(arg*t) - sinfactor*std::sin(arg*t)) + arg*(std::sin(arg*t) + sinfactor*std::cos(arg*t))));
}


// Charge function
Double_t TrackModel::Charge(Double_t *x, Double_t *par){
    Double_t t = x[0];
    
    if(t <= 0) return 0;
    // Factors
    double sigma = std::sqrt(2 * t / fRC + fDt * fDt * fdrift);  // includes transverse diffusion
    double coeff1 = std::sqrt(2 * (1 + fm * fm) / M_PI) * sigma;
    double denom = 2 * (1 + fm * fm) * sigma * sigma;
    double sqrt2_m2_sigma = std::sqrt(2 * (1 + fm * fm)) * sigma;

    // Terms
    double expRC1 = std::exp(-std::pow(-fylow  + fxright * fm + fq, 2) / denom);
    double expRC2 = std::exp(-std::pow(-fylow  + fxleft  * fm + fq, 2) / denom);
    double expRC3 = std::exp(-std::pow(-fyhigh + fxleft  * fm + fq, 2) / denom);
    double expRC4 = std::exp(-std::pow(-fyhigh + fxright * fm + fq, 2) / denom);

    double erfRC1 = (fylow  - fxleft  * fm - fq) * std::erf((-fylow  + fxleft  * fm + fq) / sqrt2_m2_sigma);
    double erfRC2 = (fyhigh - fxleft  * fm - fq) * std::erf((-fyhigh + fxleft  * fm + fq) / sqrt2_m2_sigma);
    double erfRC3 = (fylow  - fxright * fm - fq) * std::erf((-fylow  + fxright * fm + fq) / sqrt2_m2_sigma);
    double erfRC4 = (fyhigh - fxright * fm - fq) * std::erf((-fyhigh + fxright * fm + fq) / sqrt2_m2_sigma);

    // Return result
    return std::sqrt(1 + fm * fm) / (2 * fm) * (coeff1 * (expRC1 - expRC2 + expRC3 - expRC4) + erfRC1 - erfRC2 - erfRC3 + erfRC4);
}
