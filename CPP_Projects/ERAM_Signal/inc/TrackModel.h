#ifndef TRACK_MODEL_H
#define TRACK_MODEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <TMath.h>
#include <TF1.h>
#include <TF1Convolution.h>

void DrawTrackModel();

class TrackModel
{
    public:
    TrackModel(const float &impact, const float &phi, const float &RC, const float &drift, 
                const float &xpad, const float &ypad);
    virtual ~TrackModel();

    // Signal functions
    Double_t dETFdt (Double_t *x, Double_t *par);
    Double_t Charge (Double_t *x, Double_t *par);

    TF1 *ptf1_dETFdt;
    TF1 *ptf1_Charge;
    TF1 *ptf1_Signal;

    TF1Convolution *ptf1_Convolution;

    void SetAmax(const float &Amax){fAmax = Amax;}
    void Sett_Amax(const float &t_Amax){ft_Amax = t_Amax;}

    float GetAmax(){return fAmax;}
    float Gett_Amax(){return ft_Amax;}

    private:

    float fAmax =                        0; // Maximum amplitude
    float ft_Amax =                      0; // Time of maximum amplitude

    float ftmax =                       5e3; // Maximum time

    // Eletronics transfer function
    float fpeakingTime =                412;
    float ws;
    float arg;
    float expfactor;

    constexpr static float Q =          2./3.;
    const float sinfactor =             std::sqrt((2*Q-1)/(2*Q+1));
    const float flambdaG =              1.61e-4*1493*10;

    constexpr static float fSETFMAX =   0.162345;

    // Charge function
    float fimpact = 0, fphirad = 0, fdrift = 0; // inputs
    float fm = 0, fq = 0; // slope and intercept of the track
    float fxpad = 0, fypad = 0, fRC = 0; // pad parameters
    float fxleft = 0, fxright = 0, fylow = 0, fyhigh = 0; // pad borders
    constexpr static float fSXWIDTH =   11.28;
    constexpr static float fSYWIDTH =   10.19;
    float erfx, erfy;
    // Transverse diffusion coefficent
    const float fDt =                   310/pow(10,3.5); // conversion to sqrt(mm)

    // Convolution
    float fnconvolpoints =               1e4+1;
};

#endif