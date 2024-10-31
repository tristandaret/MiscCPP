#ifndef TRACK_MODEL_H
#define TRACK_MODEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <TF1.h>
#include <TF1Convolution.h>

class TrackModel
{
    public:
    TrackModel();
    virtual ~TrackModel();

    // Recompute the deposited charge based on model
    // distances in mm, angles in degrees
    double GetRealCharge(const double &ADCmax, const double &length, const double &impact, 
                         const double &phi, const double &RC, const double &drift);

    double time =                        0;
    
    private:

    // Convolution
    TF1 *ptf1_Signal;
    TF1Convolution *ptf1_Convolution;
    const double fnconvolpoints =        1e3;
    const double ftmax =                 1e3; // Maximum time

    // Eletronics transfer function
    Double_t ETF    (Double_t *x, Double_t *par);
    Double_t dETFdt (Double_t *x, Double_t *par);
    TF1 *ptf1_ETF;
    TF1 *ptf1_dETFdt;
    double fpeakingTime =                412.;
    const double Q =                     2./3.;
    const double ws =                    2/fpeakingTime;
    const double expfactor =             ws/(2*Q);
    const double arg =                   ws/2*std::sqrt(4-1/std::pow(Q,2));
    const double sinfactor =             std::sqrt((2*Q-1)/(2*Q+1));
    const double fnormelec =             4096./120./0.16233963; // 4096 max ADC for 120 fC | 0.16... max of ETF

    // Charge function
    Double_t Charge (Double_t *x, Double_t *par);
    TF1 *ptf1_Charge;
    // Track inputs
    double flength =                     0;
    double fimpact =                     0;
    double fphirad =                     0;
    double fdrift =                      0;
    double fRC =                         0;
    double fm =                          0; // slope
    double fq =                          0; // intercept
    // Computation intermediate variables
    double fsignalreal =                 0;
    double fchargemodel =                0;
    double fsignalmodel =                0;
    double famplitudeloss =              0;
    double fchargereal =                 0;
    // Pad parameters
    const double fxpad =                 0.; // pad x position
    const double fypad =                 0.; // pad y position
    const double fxleft =                -11.28/2;
    const double fxright =               11.28/2;
    const double fylow =                 -10.19/2;
    const double fyhigh =                10.19/2;
    // Transverse diffusion coefficent
    const double fDt =                   310./pow(10,3.5); // conversion to sqrt(mm)
};

#endif